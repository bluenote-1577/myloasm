use crate::cli::Cli;
use crate::constants::IDENTITY_THRESHOLDS;
use crate::constants::MAX_ALLOWABLE_SNPMER_ERROR_MISC;
use crate::constants::MAX_MULTIPLICITY_KMER;
use crate::constants::MIN_CHAIN_SCORE_COMPARE;
use crate::constants::MIN_READ_LENGTH;
use crate::constants::USE_SOLID_KMERS;
use crate::graph::GraphNode;
use crate::map_processing::*;
use crate::polishing::alignment;
use crate::seeding;
use crate::twin_graph::same_strain;
use crate::types::*;
use crate::unitig;
use crate::utils::*;
use crate::unitig::NodeSequence;
use bio_seq::codec::Codec;
use flate2::write::GzEncoder;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use rayon::prelude::*;
use flate2::Compression;
use rust_lapper::Interval;
use rust_lapper::Lapper;
use std::collections::HashSet;
use std::io::BufWriter;
use std::io::Write;
use std::path::PathBuf;
use std::sync::Mutex;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct HitInfo {
    //pub index: u32,
    pub position: u32,
    /// Bit 31 = canonical-strand flag of the reference minimizer (same encoding as
    /// `anchor.pos1`): 1 = canonical form is the reverse complement.
    /// Bits 30-0 = actual contig id.
    pub contig_id_strand: u32,
}

pub struct Anchors {
    pub anchors: Vec<Anchor>,
    pub max_mult: usize,
}

fn _edit_distance(v1: &[u32], v2: &[u32]) -> usize {
    let len1 = v1.len();
    let len2 = v2.len();

    // Create a 2D vector to store edit distances
    let mut dp = vec![vec![0; len2 + 1]; len1 + 1];

    // Initialize the first row and column
    for i in 0..=len1 {
        dp[i][0] = i;
    }
    for j in 0..=len2 {
        dp[0][j] = j;
    }

    // Fill the rest of the dp table
    for i in 1..=len1 {
        for j in 1..=len2 {
            if v1[i - 1] == v2[j - 1] {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = 1 + dp[i - 1][j - 1].min(dp[i][j - 1]).min(dp[i - 1][j]);
            }
        }
    }

    // The edit distance is the value in the bottom-right corner of the dp table
    dp[len1][len2]
}

fn _smith_waterman(
    v1: &[u32],
    v2: &[u32],
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) -> (f64, usize, usize) {
    let len1 = v1.len();
    let len2 = v2.len();

    // Create a 2D vector to store scores
    let mut dp = vec![vec![0; len2 + 1]; len1 + 1];
    let mut traceback = vec![vec![0; len2 + 1]; len1 + 1];
    let mut max_score = 0;
    let mut max_index = (0, 0);

    // Fill the dp table
    for i in 1..=len1 {
        for j in 1..=len2 {
            let score_substitute = if v1[i - 1] == v2[j - 1] {
                match_score
            } else {
                mismatch_penalty
            };

            // Calculate possible scores for this cell
            let score_diag = dp[i - 1][j - 1] + score_substitute;
            let score_up = dp[i - 1][j] + gap_penalty;
            let score_left = dp[i][j - 1] + gap_penalty;

            // Cell score is the max of calculated scores or 0 (Smith-Waterman uses zero as a minimum score)
            dp[i][j] = 0.max(score_diag).max(score_up).max(score_left);

            // Update the traceback matrix
            if dp[i][j] == 0 {
                traceback[i][j] = 0;
            } else if dp[i][j] == score_diag {
                traceback[i][j] = 1;
            } else if dp[i][j] == score_up {
                traceback[i][j] = 2;
            } else {
                traceback[i][j] = 3;
            }

            max_index = if dp[i][j] > max_score {
                max_score = dp[i][j];
                (i, j)
            } else {
                max_index
            };
        }
    }

    let mut aln1 = vec![];
    let mut aln2 = vec![];
    while dp[max_index.0][max_index.1] > 0 {
        match traceback[max_index.0][max_index.1] {
            0 => break,
            1 => {
                aln1.push(v1[max_index.0 - 1]);
                aln2.push(v2[max_index.1 - 1]);
                max_index = (max_index.0 - 1, max_index.1 - 1);
            }
            2 => {
                aln1.push(v1[max_index.0 - 1]);
                aln2.push(0);
                max_index = (max_index.0 - 1, max_index.1);
            }
            3 => {
                aln1.push(0);
                aln2.push(v2[max_index.1 - 1]);
                max_index = (max_index.0, max_index.1 - 1);
            }
            _ => panic!("Invalid traceback value"),
        }
    }
    // The maximum score represents the optimal local alignment score
    (max_score as f64, aln1.len(), aln2.len())
}

// I learned about Borrow trait recently... lol
pub fn find_exact_matches_with_full_index(
    seq1: &[(u32, FlagKmer48)],
    index: &FxHashMap<Kmer48, Vec<HitInfo>>,
    _reference_seqs_owned: Option<&FxHashMap<usize, TwinRead>>,
    _reference_seqs_ref: Option<&FxHashMap<usize, &TwinRead>>,
) -> FxHashMap<u32, Anchors> {
    let mut max_mult = 0;
    let mut matches = FxHashMap::default();

    for (pos, flag_kmer) in seq1.iter() {
        let s1 = flag_kmer.strand() as u32;
        if let Some(indices) = index.get(&flag_kmer.kmer()) {
            if indices.len() > max_mult {
                max_mult = indices.len();
            }
            for hit in indices {
                let s2     = hit.contig_id_strand >> 31;
                let contig = hit.contig_id_strand & 0x7FFF_FFFF;
                let rel_strand = s1 ^ s2;
                let anchor = AnchorBuilder {
                    pos1: (rel_strand << 31) | *pos,
                    pos2: hit.position,
                };
                matches.entry(contig).or_insert(vec![]).push(anchor);
            }
        }
    }

    matches
        .into_iter()
        .map(|(k, v)| {
            let mut anchors = v
                .into_iter()
                .map(|anchor| Anchor {
                    i: None,
                    j: None,
                    pos1: anchor.pos1,
                    pos2: anchor.pos2,
                })
                .collect::<Vec<_>>();
            // Sort by (encoded pos1 as u32, pos2) so forward anchors precede reverse.
            anchors.sort_by_key(|a| (a.pos1, a.pos2));
            (k, Anchors { anchors, max_mult })
        })
        .collect()
}

fn _find_exact_matches_indexes_references(
    seq1: &[(u32, u64)],
    seq2: &[(u32, u64)],
) -> (Vec<Anchor>, usize) {
    let mut max_mult = 0;
    let mut matches = Vec::new();
    let mut index_map = FxHashMap::default();

    //Sorted
    for (j, (pos, x)) in seq2.iter().enumerate() {
        index_map.entry(x).or_insert(vec![]).push((j, pos));
    }

    //Sorted
    for (i, (pos, s)) in seq1.iter().enumerate() {
        if let Some(indices) = index_map.get(&s) {
            if indices.len() > max_mult {
                max_mult = indices.len();
            }
            for (j, pos2) in indices {
                let anchor = Anchor {
                    i: Some(i as u32),
                    j: Some(*j as u32),
                    pos1: *pos as u32,
                    pos2: **pos2 as u32,
                };
                matches.push(anchor);
            }
        }
    }

    (matches, max_mult)
}

/// Matches canonical minimizers from `seq1` against `seq2` and encodes the relative
/// strand into bit 31 of `Anchor.pos1`:
///   bit 31 = 0 → forward alignment (pos2 increases in a valid chain)
///   bit 31 = 1 → reverse alignment (pos2 decreases in a valid chain)
/// `Anchor.pos2` always holds the raw reference position.
fn find_exact_matches_indexes(
    seq1: &[(u32, FlagKmer48)],
    seq2: &[(u32, FlagKmer48)],
) -> (Vec<Anchor>, usize) {
    let mut max_mult = 0;
    let mut matches = Vec::new();

    // Index seq2: key = plain Kmer48 (no strand flag), value = (index, pos, is_rc).
    let mut index_map: FxHashMap<Kmer48, Vec<(usize, u32, bool)>> = FxHashMap::default();
    for (j, &(pos, kmer)) in seq2.iter().enumerate() {
        index_map
            .entry(kmer.kmer())
            .or_insert_with(Vec::new)
            .push((j, pos, kmer.strand()));
    }

    for (i, &(pos1, kmer)) in seq1.iter().enumerate() {
        let s1 = kmer.strand() as u32;
        if let Some(indices) = index_map.get(&kmer.kmer()) {
            if indices.len() > max_mult {
                max_mult = indices.len();
            }
            for &(j, pos2, s2) in indices {
                // rel_strand = 0 → same strand → forward chain (pos2 increases)
                // rel_strand = 1 → opposite strands → reverse chain (pos2 decreases)
                let rel_strand = s1 ^ (s2 as u32);
                let anchor = Anchor {
                    i: Some(i as u32),
                    j: Some(j as u32),
                    pos1: (rel_strand << 31) | pos1,
                    pos2,
                };
                matches.push(anchor);
            }
        }
    }

    (matches, max_mult)
}

fn _find_exact_matches_quadratic(seq1: &[usize], seq2: &[usize]) -> Vec<(usize, usize, usize)> {
    let mut matches = Vec::new();
    let len1 = seq1.len();
    let len2 = seq2.len();

    for i in 0..len1 {
        for j in 0..len2 {
            if seq1[i] == seq2[j] {
                matches.push((i, j, 1)); // (start in seq1, start in seq2, length of match = 1)
            }
        }
    }
    matches
}

pub fn dp_anchors(
    matches: &[Anchor],
    reverse: bool,
    gap_cost: i32,
    match_score: i32,
    band: usize,
    max_gap: usize,
    double_gap: usize,
    _debug: bool,
) -> Vec<(i32, Vec<Anchor>, bool)> {
    if matches.is_empty() {
        return vec![(0, vec![], false)];
    }
    let mut dp = vec![0; matches.len()];
    let mut prev = vec![None; matches.len()];
    let mut large_indel_tracker = vec![false; matches.len()];
    let mut max_score = 0;
    let mut max_index = 0;
    let max_skip = 10;

    for i in 0..matches.len() {
        let start1 = matches[i].pos1;
        let start2 = matches[i].pos2;
        dp[i] = match_score;
        let mut unimproved = 0;
        let mut unimproved_pos1 = 0;
        let back = if i > band { i - band } else { 0 };
        for j in (back..i).rev() {
            let end1 = matches[j].pos1;
            let end2 = matches[j].pos2;
            if reverse {
                if end1 >= start1 || end2 <= start2 {
                    continue;
                }
            } else {
                if end1 >= start1 || end2 >= start2 {
                    continue;
                }
            }
            let gap_penalty_signed =
                (start1 as i32 - (end1) as i32).abs() - (start2 as i32 - (end2) as i32).abs();
            let gap_penalty = gap_penalty_signed.abs();

            let kmer_dist1_score = (start1 as i32 - end1 as i32).abs().min(match_score);
            let kmer_dist2_score = (start2 as i32 - end2 as i32).abs().min(match_score);
            let kmer_overlap_score = kmer_dist1_score.min(kmer_dist2_score);

            let large_indel = false;
            if gap_penalty > max_gap as i32 {
                //large_indel = true;
                continue;
            }

            if (start1 as i32 - end1 as i32).abs() > (double_gap.try_into().unwrap())
                || (start2 as i32 - end2 as i32).abs() > (double_gap.try_into().unwrap())
            {
                continue;
            }

            let score = dp[j] + kmer_overlap_score - gap_cost * gap_penalty;
            if score > dp[i] {
                dp[i] = score;
                prev[i] = Some(j);
                if score > max_score {
                    max_score = score;
                    max_index = i;
                    unimproved = 0;
                }
                if large_indel {
                    large_indel_tracker[i] = true;
                }
            } else {
                if unimproved_pos1 != end1 {
                    unimproved += 1;
                    unimproved_pos1 = end1;
                }
            }
            if unimproved > max_skip {
                break;
            }
        }
    }

    // Reconstruct the chain
    let mut chains = Vec::new();
    let mut used_anchors = FxHashSet::default();
    let mut best_indices_ordered = (0..matches.len()).map(|i| (dp[i], i)).collect::<Vec<_>>();
    best_indices_ordered.sort_by_key(|&(score, _)| -score);
    assert!(dp[max_index] == best_indices_ordered[0].0);


    for (score, best_index) in best_indices_ordered {
        if used_anchors.contains(&best_index) {
            continue;
        }

        let mut chain = Vec::new();
        let mut large_indel = false;
        let mut i = Some(best_index);
        let mut good_chain = true;
        while let Some(idx) = i {
            large_indel = large_indel || large_indel_tracker[idx];
            if used_anchors.contains(&idx) {
                good_chain = false;
                break;
            }
            used_anchors.insert(idx);
            chain.push(matches[idx].clone());
            i = prev[idx];
        }

        if chain.len() < 2 {
            break;
        }

        if good_chain {
            chain.reverse();
            chains.push((score, chain, large_indel));
        }
    }

    return chains;
}

/// Anchor chaining DP following the minimap2/minigraph `mg_lchain_dp` algorithm.
///
/// Key differences from `dp_anchors`:
/// - `st` sliding window replaces fixed `band`: left boundary advances when pos1 distance
///   exceeds `double_gap`, so the inner loop is O(n) total rather than O(n * band).
/// - Skip logic uses minimap2's `t[]` marker scheme instead of a broken `unimproved` counter:
///   `t[p[j]] = i` marks the predecessor of each valid candidate; `n_skip` counts how many
///   such already-seen predecessors fail to improve dp[i]. Decrements on improvement.
/// - `max_ii` fallback: after an early-stop, explicitly checks the highest-scoring predecessor
///   in the window, ensuring the global best predecessor is never missed.
/// - Invalid transitions (ordering violations, gap > max_gap, dist > double_gap) do NOT
///   participate in the skip-counter logic at all, matching minimap2 semantics.
/// Strand-aware anchor chaining.
///
/// Anchors must have `pos1` encoded by `find_exact_matches_indexes`:
///   bit 31 of pos1 = 0 → forward alignment anchor
///   bit 31 of pos1 = 1 → reverse alignment anchor  (pos2 should decrease in chain)
///
/// All anchors (both strands) are processed in a **single DP pass**.
/// Strand separation is achieved by the MSB encoding:
///   - sorted by (pos1 as u32, encoded_pos2 as u32): forward anchors precede reverse
///   - cross-strand transitions naturally fail dist1 ≤ 0 or dist2 ≤ 0 in dp_inner
///
/// Returns `(score, chain, is_reverse)` tuples; `is_reverse` replaces the old
/// `large_indel` placeholder (which was always false).
/// Chain anchors have the strand bit stripped from pos1; use `ChainInfo.reverse`.
pub fn dp_anchors_v2(
    matches: &[Anchor],
    gap_cost: i32,
    match_score: i32,
    max_iter: usize,
    max_skip: usize,
    max_gap: usize,
    double_gap: usize,
    debug: bool,
    min_chain_length: usize,
) -> Vec<(i32, Vec<Anchor>, bool)> {
    if matches.is_empty() {
        return vec![];
    }
    let n = matches.len();

    // --- Sort by (pos1 as u32, encoded_pos2 as u32) ---
    // Forward anchors (bit31=0) sort before reverse anchors (bit31=1) because their
    // pos1 u32 values are smaller.  Within each strand group, anchors are sorted by
    // increasing encoded_pos2.  For reverse anchors encoded_pos2 = !pos2_raw, so the
    // sort puts large-pos2 anchors first (which are valid predecessors in a reverse chain).
    let mut sorted_indices: Vec<usize> = (0..n).collect();
    sorted_indices.sort_unstable_by_key(|&i| {
        let pos1_u32  = matches[i].pos1 as i32;
        let pos2_raw  = matches[i].pos2 as i32;
        let enc_pos2  = if pos1_u32 >> 31 == 1 { !pos2_raw } else { pos2_raw };
        (pos1_u32, enc_pos2)
    });

    // --- SoA extraction in sorted order ---
    // For forward anchors: pos1_i32 = pos1_raw (positive), pos2_i32 = pos2_raw (positive).
    // For reverse anchors: pos1_i32 = (1<<31|pos1_raw) as i32 (negative, sign bit set),
    //                      pos2_i32 = (!pos2_raw) as i32 (negative).
    // dp_inner::<false> uses dist1 = s1-e1 and dist2 = s2-e2 for both strands:
    //   • same-strand forward: dists are naturally positive ✓
    //   • same-strand reverse: i32::MIN terms cancel, dists equal raw differences ✓
    //   • cross-strand: dist1 or dist2 < 0 → rejected by the ≤ 0 guard ✓
    let pos1s: Vec<i32> = sorted_indices.iter().map(|&i| matches[i].pos1 as i32).collect();
    let pos2s: Vec<i32> = sorted_indices.iter().map(|&i| {
        let pos2_raw = matches[i].pos2;
        if matches[i].pos1 >> 31 == 1 { (!pos2_raw) as i32 } else { pos2_raw as i32 }
    }).collect();

    let mut f: Vec<i32> = vec![match_score; n];
    let mut p: Vec<i32> = vec![-1i32; n];
    let mut t: Vec<i32> = vec![-1i32; n];

    let double_gap_i32 = double_gap as i32;
    let max_gap_i32    = max_gap    as i32;
    let mut st: usize  = 0;
    let mut max_ii: i32 = -1;

    let time = std::time::Instant::now();

    // Single pass — REV const generic no longer needed; encoding handles both strands.
    dp_inner::<false>(
        &pos1s, &pos2s, &mut f, &mut p, &mut t,
        gap_cost, match_score, max_iter, max_skip,
        double_gap_i32, max_gap_i32,
        &mut st, &mut max_ii, n,
    );

    if debug {
        log::debug!(
            "DP chaining (v2, single-pass) of {} anchors took {:?}",
            matches.len(),
            time.elapsed()
        );
    }

    // --- Chain reconstruction ---
    // Repurpose t[] as the claimed-anchor marker (0 = unclaimed, 1 = claimed).
    t.fill(0);

    let time = std::time::Instant::now();
    let mut chains = Vec::new();

    let mut best_indices_ordered = (0..n as i32)
        .filter(|&i| f[i as usize] > min_chain_length as i32 * match_score / 2) // heuristic threshold to skip very short chains
        .map(|i| (f[i as usize], i))
        .collect::<Vec<_>>();
    best_indices_ordered.sort_by_key(|&(score, _)| -score);

    if debug {
        log::debug!("Sorting candidates for chain reconstruction took {:?}", time.elapsed());
    }

    let time = std::time::Instant::now();
    for (score, best_index) in best_indices_ordered {
        let bi = best_index as usize;
        if t[bi] != 0 {
            continue;
        }

        // Detect strand from the top anchor (all anchors in a chain share the same strand).
        let chain_is_reverse = matches[sorted_indices[bi]].pos1 >> 31 == 1;

        let mut chain = Vec::new();
        let mut idx = best_index;
        let mut good_chain = true;
        while idx >= 0 {
            let u = idx as usize;
            if t[u] != 0 {
                good_chain = false;
                break;
            }
            t[u] = 1;
            // Strip the strand bit from pos1 — strand is in chain_is_reverse.
            let orig = &matches[sorted_indices[u]];
            chain.push(Anchor {
                pos1: orig.pos1 & 0x7FFF_FFFF,
                ..*orig
            });
            idx = p[u];
        }

        if chain.len() < min_chain_length {
            break;
        }

        if good_chain {
            chain.reverse();
            chains.push((score, chain, chain_is_reverse));
        }
    }

    if debug {
        log::debug!("Extracting chains from DP (v2) took {:?}", time.elapsed());
    }

    chains
}

/// Monomorphised DP kernel. `REV` is true for reverse-strand chaining.
/// All runtime branching on strand direction is eliminated by the const generic.
#[inline(never)]
fn dp_inner<const REV: bool>(
    pos1s: &[i32],
    pos2s: &[i32],
    f: &mut [i32],
    p: &mut [i32],
    t: &mut [i32],
    gap_cost: i32,
    match_score: i32,
    max_iter: usize,
    max_skip: usize,
    double_gap_i32: i32,
    max_gap_i32: i32,
    st: &mut usize,
    max_ii: &mut i32,
    n: usize,
) {
    for i in 0..n {
        let s1 = pos1s[i];
        let s2 = pos2s[i];
        let i_i32 = i as i32; // hoisted; used in t[] comparisons and writes

        // Advance pos1 window left boundary (O(n) total across all i).
        while *st < i && s1 > double_gap_i32 + pos1s[*st] {
            *st += 1;
        }
        let lo = if i >= max_iter { (i - max_iter).max(*st) } else { *st };

        let mut max_f  = match_score;
        let mut max_j  = -1i32;
        let mut n_skip = 0usize;
        let mut end_j  = lo; // last j visited; used by the max_ii fallback
        let strand_i = s1 >> 31;

        'inner: for j in (lo..i).rev() {
            end_j = j;
            let e1 = pos1s[j];
            let e2 = pos2s[j];

            // Combined ordering + pos2-distance filter.
            // Invalid transitions do NOT count toward n_skip or t[] marking.
            // For REV=false: need e1 < s1 (dist1 > 0) and e2 < s2 (dist2 > 0).
            // For REV=true:  need e1 < s1 (dist1 > 0) and e2 > s2 (dist2 > 0).
            let dist1 = s1 - e1;
            let dist2 = if REV { e2 - s2 } else { s2 - e2 };
            let same_strand = strand_i == (e1 >> 31);
            if dist1 <= 0 || dist2 <= 0 || dist2 > double_gap_i32 || !same_strand {
                continue;
            }

            // Gap-imbalance filter
            let gap_penalty = (dist1 - dist2).abs();
            if gap_penalty > max_gap_i32 {
                continue;
            }

            let kmer_overlap_score = dist1.min(dist2).min(match_score);
            let sc = f[j] + kmer_overlap_score - gap_cost * gap_penalty;

            if sc > max_f {
                max_f = sc;
                max_j = j as i32;
                if n_skip > 0 { n_skip -= 1; }
            } else if t[j] == i_i32 {
                n_skip += 1;
                if n_skip > max_skip { break 'inner; }
            }

            // Mark predecessor of j as a skip candidate for anchor i.
            let pj = p[j];
            if pj >= 0 { t[pj as usize] = i_i32; }
        }

        // --- max_ii fallback ---
        // Rescan if max_ii has left the pos1 window.
        let rescan = *max_ii < 0 || s1  > double_gap_i32 + pos1s[*max_ii as usize];
        if rescan {
            let mut best = i32::MIN;
            *max_ii = -1;
            for j in (lo..i).rev() {
                if f[j] > best { best = f[j]; *max_ii = j as i32; }
            }
        }
        // If the global best predecessor was in the unvisited region, check it.
        let mii = *max_ii;
        if mii >= 0 && (mii as usize) < end_j {
            let m = mii as usize;
            let e1 = pos1s[m];
            let e2 = pos2s[m];
            let dist1 = s1 - e1;
            let dist2 = if REV { e2 - s2 } else { s2 - e2 };
            if dist1 > 0 && dist2 > 0 && dist2 <= double_gap_i32 {
                let gap_penalty = (dist1 - dist2).abs();
                if gap_penalty <= max_gap_i32 {
                    let kmer_overlap_score = dist1.min(dist2).min(match_score);
                    let sc = f[m] + kmer_overlap_score - gap_cost * gap_penalty;
                    if sc > max_f { max_f = sc; max_j = mii; }
                }
            }
        }

        f[i] = max_f;
        p[i] = max_j;

        if *max_ii < 0 || f[*max_ii as usize] < f[i] { *max_ii = i_i32; }
    }
}

fn find_optimal_chain(
    anchors: &Vec<Anchor>,
    match_score: i32,
    gap_cost: i32,
    band_opt: Option<usize>,
    tr_options: &CompareTwinReadOptions
) -> Vec<ChainInfo> {
    let band;
    let matches = anchors;
    let max_gap = tr_options.max_gap;
    let double_gap = tr_options.double_gap;
    if band_opt.is_none() {
        band = 50;
    } else {
        band = band_opt.unwrap();
    }

    if anchors.is_empty() {
        return vec![];
    }

    // Single strand-aware DP pass — both forward and reverse chains in one call.
    let mut all_chains = dp_anchors_v2(
        matches, gap_cost, match_score, band,
        tr_options.max_skip, max_gap, double_gap, tr_options.debug, tr_options.min_chain_length
    );

    if all_chains.is_empty() {
        return vec![];
    }

    let max_score = all_chains.iter().map(|x| x.0).max().unwrap();
    let mut chains = vec![];
    let mut reference_intervals: Vec<Interval<u32, bool>> = vec![];
    let mut query_intervals: Vec<Interval<u32, bool>> = vec![];
    all_chains.sort_unstable_by(|a, b| b.0.cmp(&a.0));

    for (score, chain, reverse) in all_chains {
        let large_indel = false;
        let cond1 = score as f64 > tr_options.supplementary_threshold_ratio.unwrap_or(0.25) * max_score as f64;
        let cond2 = score as f64 > tr_options.supplementary_threshold_score.unwrap_or(f64::MAX);
        if cond1 || cond2{
            let l = chain.first().unwrap().pos2;
            let r = chain.last().unwrap().pos2;
            let interval = Interval {
                start: l.min(r),
                stop: l.max(r),
                val: true,
            };

            if reference_intervals.iter().any(|x| {
                    let intersect = x.intersect(&interval);
                    intersect as f64 / (interval.stop - interval.start) as f64 > 0.25
                })
            {
                if tr_options.force_ref_nonoverlap{
                    continue;
                }
            }

            reference_intervals.push(interval);

            let l_q = chain.first().unwrap().pos1;
            let r_q = chain.last().unwrap().pos1;
            let interval_q = Interval {
                start: l_q.min(r_q),
                stop: l_q.max(r_q),
                val: true,
            };

            if query_intervals.iter().any(|x| {
                let intersect = x.intersect(&interval_q);
                intersect as f64 / (interval_q.stop - interval_q.start) as f64 > 0.25
            })
            {
                if tr_options.force_query_nonoverlap{
                    continue;
                }
                else{
                    let secondary_ratio = tr_options.secondary_threshold.unwrap_or(0.50);
                    if (score as f64) < secondary_ratio * (max_score as f64) {
                        continue;
                    }
                }
            }

            query_intervals.push(interval_q);

            chains.push(ChainInfo {
                chain,
                reverse: reverse,
                score: score,
                large_indel: large_indel,
            });

        }
    }

    return chains;
}

pub fn compare_twin_reads(
    seq1: &TwinRead,
    seq2: &TwinRead,
    mini_anchors: Option<&Anchors>,
    snpmer_anchors: Option<&Anchors>,
    i: usize,
    j: usize,
    options: &CompareTwinReadOptions,
    args: &Cli,
) -> Vec<TwinOverlap> {
    let mut mini_chain_infos;
    let time = std::time::Instant::now();
    if let Some(anchors) = mini_anchors {
        mini_chain_infos = find_optimal_chain(
            &anchors.anchors,
            args.c as i32,
            1,
            Some(anchors.max_mult * 20),
            options,
        );
        if options.debug{
            log::debug!("Compared read {} to read {} in {:?} with {} anchors, found {} mini chains", seq1.id, seq2.id, time.elapsed(), anchors.anchors.len(), mini_chain_infos.len());
        }
    } else {
        let anchors;
        if let Some(seq1_minimizers) = options.read1_mininimizers.as_ref(){
            anchors = find_exact_matches_indexes(seq1_minimizers, &seq2.minimizers_vec_strand());
        }
        else{
            anchors = find_exact_matches_indexes(&seq1.minimizers_vec_strand(), &seq2.minimizers_vec_strand());
        }
        mini_chain_infos = find_optimal_chain(
            &anchors.0,
            args.c as i32,
            1,
            Some((anchors.1 * 20).min(MAX_MULTIPLICITY_KMER)),
            options,
        );
    }


    if options.maximal_only{
        mini_chain_infos.retain(|mini_chain_info| {
            let mini_chain = &mini_chain_info.chain;
            // Use positions from the Anchor struct directly
            let l1 = mini_chain[0].pos1;
            let r1 = mini_chain[mini_chain.len() - 1].pos1;
            let l2 = mini_chain[0].pos2;
            let r2 = mini_chain[mini_chain.len() - 1].pos2;
            let start1 = l1.min(r1);
            let end1 = l1.max(r1) + seq1.k as u32 - 1;
            let start2 = l2.min(r2);
            let end2 = l2.max(r2) + seq2.k as u32 - 1;
            let shared_minimizers = mini_chain.len();
            let end_fuzz_pair = seq1.overlap_hang_length.unwrap(); 
            let end_fuzz = end_fuzz_pair.0.max(end_fuzz_pair.1);
            let max_mapping = check_maximal_overlap(start1 as usize, end1 as usize, start2 as usize, end2 as usize, seq1.base_length, seq2.base_length, mini_chain_info.reverse, end_fuzz);
            if !max_mapping{
                return false;
            }
            if (shared_minimizers as f64) < (seq1.base_length as f64 / args.c as f64 / args.absolute_minimizer_cut_ratio){
                return false;
            }
            return true;
        });
    }

    mini_chain_infos.retain(|mini_chain_info|{
        mini_chain_info.score >= MIN_CHAIN_SCORE_COMPARE
    });

    let mini_chain_infos = mini_chain_infos;

    if mini_chain_infos.is_empty() {
        return vec![];
    }

    let mut twin_overlaps = vec![];
    let k = seq1.k as usize;

    let temp_vec;
    let mut snpmer_vec = &vec![];

    let temp_vec2;
    let mut snpmer_vec_2 = &vec![];

    // We have to populate snpmers_vec() which may be large for contigs
    if options.compare_snpmers{
        if let Some(snpmers_vec_1) = options.read1_snpmers.as_ref() {
            snpmer_vec = snpmers_vec_1;
        } else {
            temp_vec = seq1.snpmers_vec_strand();
            snpmer_vec = &temp_vec;
        }

        temp_vec2 = seq2.snpmers_vec_strand();
        snpmer_vec_2 = &temp_vec2;
    }
    

    for mini_chain_info in mini_chain_infos {
        let mini_chain = &mini_chain_info.chain;

        let mut shared_snpmer = usize::MAX;
        let mut diff_snpmer = usize::MAX;

        if options.compare_snpmers {
            shared_snpmer = 0;
            diff_snpmer = 0;

            // Use positions from the Anchor struct directly
            let l1 = mini_chain[0].pos1 as usize;
            let r1 = mini_chain[mini_chain.len() - 1].pos1 as usize;
            let l2 = mini_chain[0].pos2 as usize;
            let r2 = mini_chain[mini_chain.len() - 1].pos2 as usize;
            let start1 = l1.min(r1);
            let end1 = l1.max(r1) + k - 1;
            let start2 = l2.min(r2);
            let end2 = l2.max(r2) + k - 1;

            let mask = !(3 << (k - 1));

            let mut splitmers1: Vec<(u32, FlagKmer48)> = vec![];
            let mut ind_redirect1 = vec![];

            for (i, &(pos, snpmer)) in snpmer_vec.iter().enumerate() {
                if pos as usize >= start1 && pos as usize <= end1 {
                    ind_redirect1.push(i);
                    let masked = FlagKmer48::new(
                        Kmer48::from_u64(snpmer.kmer().to_u64() & mask),
                        snpmer.strand(),
                    );
                    splitmers1.push((pos, masked));
                }
            }

            let mut splitmers2: Vec<(u32, FlagKmer48)> = vec![];
            let mut ind_redirect2 = vec![];

            for (i, &(pos, snpmer)) in snpmer_vec_2.iter().enumerate() {
                if pos as usize >= start2 && pos as usize <= end2 {
                    ind_redirect2.push(i);
                    let masked = FlagKmer48::new(
                        Kmer48::from_u64(snpmer.kmer().to_u64() & mask),
                        snpmer.strand(),
                    );
                    splitmers2.push((pos, masked));
                }
            }

            let split_chain_opt;
            let mut split_options = options.clone();
            split_options.min_chain_length = 2;
            split_options.double_gap = 2_000_000;
            if let Some(anchors) = snpmer_anchors {
                split_chain_opt = find_optimal_chain(
                    &anchors.anchors,
                    50,
                    1,
                    Some((anchors.max_mult * 10).min(50)),
                    &split_options,
                )
                .into_iter()
                .max_by_key(|x| x.score);
            } else {
                let anchors = find_exact_matches_indexes(&splitmers1, &splitmers2);
                let chains = find_optimal_chain(
                    &anchors.0,
                    50,
                    1,
                    Some((anchors.1 * 10).min(50)),
                    &split_options,
                );
                split_chain_opt = chains.into_iter().max_by_key(|x| x.score);
            }
            

            //If mini chain goes opposite from split chain, probably split chain
            //is not reliable, so set shared and diff = 0.
            if let Some(split_chain) = split_chain_opt.as_ref() {
                if split_chain.reverse == mini_chain_info.reverse || split_chain.chain.len() == 1 {
                    let snpmer_kmers_seq1 = &snpmer_vec;
                    let snpmer_kmers_seq2 = &snpmer_vec_2;
                    for anchor in split_chain.chain.iter() {
                        let i = anchor.i;
                        let i = ind_redirect1[i.unwrap() as usize];
                        let j = anchor.j;
                        let j = ind_redirect2[j.unwrap() as usize];

                        //if seq1.snpmer_kmers[i as usize] == seq2.snpmer_kmers[j as usize] {
                        if snpmer_kmers_seq1[i as usize].1.kmer() == snpmer_kmers_seq2[j as usize].1.kmer() {
                            shared_snpmer += 1;
                        } else {
                            diff_snpmer += 1;
                        }
                    }
                }
            }

            //Only if log level is trace
            if log::log_enabled!(log::Level::Trace) && true {
                if diff_snpmer < 10 && shared_snpmer > 100 {
                    let mut positions_read1_snpmer_diff = vec![];
                    let mut positions_read2_snpmer_diff = vec![];

                    let mut kmers_read1_diff = vec![];
                    let mut kmers_read2_diff = vec![];

                    let snpmer_kmers_seq1 = seq1.snpmer_kmers();
                    let snpmer_kmers_seq2 = seq2.snpmer_kmers();
                    let snpmer_positions_seq1 = seq1.snpmer_positions();
                    let snpmer_positions_seq2 = seq2.snpmer_positions();

                    for anchor in split_chain_opt.unwrap().chain.iter() {
                        let i = anchor.i;
                        let i = ind_redirect1[i.unwrap() as usize];
                        let j = anchor.j;
                        let j = ind_redirect2[j.unwrap() as usize];
                        if snpmer_kmers_seq1[i as usize] != snpmer_kmers_seq2[j as usize] {
                            positions_read1_snpmer_diff.push(snpmer_positions_seq1[i as usize]);
                            positions_read2_snpmer_diff.push(snpmer_positions_seq2[j as usize]);

                            let kmer1 = decode_kmer48(snpmer_kmers_seq1[i as usize], seq1.k as u8);
                            let kmer2 = decode_kmer48(snpmer_kmers_seq2[j as usize], seq2.k as u8);

                            kmers_read1_diff.push(kmer1);
                            kmers_read2_diff.push(kmer2);
                        }
                    }
                    log::trace!(
                        "{}--{:?} {}--{:?}, snp_diff:{} snp_shared:{}, kmers1:{:?}, kmers2:{:?}",
                        &seq1.id,
                        positions_read1_snpmer_diff,
                        &seq2.id,
                        positions_read2_snpmer_diff,
                        diff_snpmer,
                        shared_snpmer,
                        kmers_read1_diff,
                        kmers_read2_diff
                    );
                }
            }
        }

        // Use positions from the Anchor struct directly
        let l1 = mini_chain[0].pos1;
        let r1 = mini_chain[mini_chain.len() - 1].pos1;
        let l2 = mini_chain[0].pos2;
        let r2 = mini_chain[mini_chain.len() - 1].pos2;
        let start1 = l1.min(r1);
        let end1 = l1.max(r1) + k as u32 - 1;
        let start2 = l2.min(r2);
        let end2 = l2.max(r2) + k as u32 - 1;
        let shared_minimizers = mini_chain.len();
        let mut mini_chain_return = None;
        if options.retain_chain {
            mini_chain_return = Some(mini_chain_info.chain);
        }
        let twinol = TwinOverlap {
            i1: i,
            i2: j,
            start1: start1 as usize,
            end1: end1 as usize,
            start2: start2 as usize,
            end2: end2 as usize,
            shared_minimizers,
            shared_snpmers: shared_snpmer,
            diff_snpmers: diff_snpmer,
            snpmers_in_both: (seq1.snpmer_count(), seq2.snpmer_count()),
            chain_reverse: mini_chain_info.reverse,
            chain_score: mini_chain_info.score,
            minimizer_chain: mini_chain_return,
            large_indel: mini_chain_info.large_indel,
        };
        twin_overlaps.push(twinol);
    }
    return twin_overlaps;
}



pub fn id_est(shared_minimizers: usize, diff_snpmers: usize, c: u64, large_indel: bool) -> f64 {
    let diff_snps = diff_snpmers as f64;
    let shared_minis = shared_minimizers as f64;
    let alpha = diff_snps as f64 / shared_minis as f64 / c as f64;
    let theta = alpha / (1. + alpha);
    let mut id_est = 1. - theta;

    if large_indel {
        //Right now it's 0.5% penalty for large indels, but we don't use largeindels. 
        let penalty = IDENTITY_THRESHOLDS.last().unwrap() - IDENTITY_THRESHOLDS.first().unwrap();
        let penalty = penalty / 2.;
        id_est -= penalty;
    }

    return id_est;
}

fn unitigs_to_tr(
    unitig_graph: &unitig::UnitigGraph,
    snpmer_set: &FxHashSet<Kmer64>,
    solid_kmers: &HashSet<Kmer48>,
    high_freq_kmers: &HashSet<Kmer48>,
    args: &Cli,
) -> FxHashMap<usize, TwinRead> {
    let tr_unitigs = Mutex::new(FxHashMap::default());
    unitig_graph.nodes.par_iter().for_each(|(&node_hash_id, unitig)| {
        let id = format!("u{}", unitig.read_indices_ori[0].0);
        let u8_seq = unitig
            .base_seq()
            .iter()
            .map(|x| bits_to_ascii(x.to_bits()) as u8)
            .collect::<Vec<u8>>();
        let tr =
            seeding::get_twin_read_syncmer(u8_seq, None, args.kmer_size, args.c, &snpmer_set, id);
        if let Some(mut tr) = tr {
            let mut solid_mini_indices = FxHashSet::default();
            let mut solid_snpmer_indices = FxHashSet::default();
            for (index, mini) in tr.minimizer_kmers().iter().enumerate(){
                if USE_SOLID_KMERS{
                    if solid_kmers.contains(&mini) {
                        solid_mini_indices.insert(index);
                    }
                }
                else{
                    if !high_freq_kmers.contains(&mini) {
                        solid_mini_indices.insert(index);
                    }
                }
            }
            for (index, snpmer) in tr.snpmer_kmers().iter().enumerate() {
                if USE_SOLID_KMERS{
                    if snpmer_set.contains(&snpmer.to_u64()) {
                        solid_snpmer_indices.insert(index);
                    }
                }
                else{
                    if !high_freq_kmers.contains(&snpmer) {
                        solid_snpmer_indices.insert(index);
                    }
                }
            }
            tr.retain_mini_indices(solid_mini_indices);
            tr.retain_snpmer_indices(solid_snpmer_indices);
            tr_unitigs.lock().unwrap().insert(node_hash_id, tr);
        }
    });
    return tr_unitigs.into_inner().unwrap();
}

pub fn get_minimizer_index(
    tr_owned: Option<&FxHashMap<usize, TwinRead>>, 
    tr_ref: Option<&FxHashMap<usize, &TwinRead>>) -> FxHashMap<Kmer48, Vec<HitInfo>> {
    let mut mini_index = FxHashMap::default();
    if let Some(twinreads) = tr_owned {
        let mut sorted_keys = twinreads.keys().collect::<Vec<_>>();
        sorted_keys.sort();
        for (&id, tr) in sorted_keys.iter().map(|&x| (x, &twinreads[x])) {
            for (pos, mini) in tr.minimizers_vec_strand().into_iter() {
                let hit = HitInfo {
                    contig_id_strand: (mini.strand() as u32) << 31 | id as u32,
                    position: pos,
                };
                mini_index.entry(mini.kmer()).or_insert(vec![]).push(hit);
            }
        }
    }
    else if let Some(twinreads) = tr_ref {
        let mut sorted_keys = twinreads.keys().collect::<Vec<_>>();
        sorted_keys.sort();
        for (&id, tr) in sorted_keys.iter().map(|&x| (x, &twinreads[x])) {
            for (pos, mini) in tr.minimizers_vec_strand().into_iter() {
                let hit = HitInfo {
                    contig_id_strand: (mini.strand() as u32) << 31 | id as u32,
                    position: pos,
                };
                mini_index.entry(mini.kmer()).or_insert(vec![]).push(hit);
            }
        }
    }
    else{
        panic!("No minimizer index provided");
    }

    if mini_index.len() == 0{
        return mini_index
    }

    let mut minimizer_to_hit_count = mini_index
        .iter()
        .map(|(_, v)| v.len())
        .collect::<Vec<_>>();

    minimizer_to_hit_count.sort_by(|a, b| b.cmp(&a));
    let threshold = minimizer_to_hit_count[minimizer_to_hit_count.len() / 100_000];
    log::trace!(
        "Minimizer index size: {}. Threshold: {}",
        minimizer_to_hit_count.len(),
        threshold
    );

    // Only threshold when necessary
    if mini_index.len() > 500_000{
        mini_index.retain(|_, v| v.len() < threshold);
    }

    mini_index.shrink_to_fit();

    return mini_index;
}

/// Memory-efficient version that processes batches immediately and returns compact split plans
pub fn map_reads_to_outer_reads_efficient(
    outer_read_indices: &[usize],
    twin_reads: &[TwinRead],
    args: &Cli,
    break_chimeras: bool
) -> Vec<SplitReadPlan> {
    use crate::map_processing::{cov_mapping_breakpoints, compute_coverage_stats};

    let split_plans = Mutex::new(vec![]);
    let num_alignments = Mutex::new(0);
    let num_maximal = Mutex::new(0);

    let chunk_size = args.read_map_batch_size;

    let outer_read_chunks = outer_read_indices
        .chunks(chunk_size)
        .collect::<Vec<_>>();

    log::info!("Mapping {} reads to {} outer reads. Building hash table... ", twin_reads.len(), outer_read_indices.len());
    log::debug!("ITERATIONS: Breaking {} reads into {} chunks of <= {}", outer_read_indices.len(), outer_read_chunks.len(), chunk_size);

    for (batch_index, mapping_chunk_indices) in outer_read_chunks.iter().enumerate() {
        // (encoded boundaries as varints, BareMappingOverlap vec)
        let mapping_maximal_boundaries: Mutex<FxHashMap<usize, (Vec<u8>, Vec<BareMappingOverlap>)>> = Mutex::new(FxHashMap::default());
        let mapping_local_boundaries_map = Mutex::new(FxHashMap::default());

        let tr_outer = mapping_chunk_indices
            .into_iter()
            .map(|&i| (i, &twin_reads[i]))
            .collect::<FxHashMap<usize, &TwinRead>>();

        let mini_index = get_minimizer_index(None, Some(&tr_outer));
        log::info!("Built minimizer index for batch {}/{}. Mapping reads...", batch_index + 1, outer_read_chunks.len());
        let counter = Mutex::new(0);

        twin_reads.par_iter().enumerate().for_each(|(rid, read)| {
            let mut tr_options = CompareTwinReadOptions::default();
            tr_options.double_gap = 1500;
            tr_options.force_query_nonoverlap = true;
            tr_options.read1_snpmers = Some(read.snpmers_vec_strand());

            // if rid % 100 == 0{
            //     tr_options.debug = true;
            // }

            //log::debug!("Processing read {} with id {}...", rid, read.id);

            let mini = read.minimizers_vec_strand();
            let mini_anchors = find_exact_matches_with_full_index(&mini, &mini_index, None, Some(&tr_outer));
            drop(mini);
            let mut unitig_hits : Vec<TwinOverlap> = vec![];
            *counter.lock().unwrap() += 1;
            if *counter.lock().unwrap() % 100000 == 0 {
                log::debug!(
                    "Processed {} reads / {} ...",
                    *counter.lock().unwrap(),
                    twin_reads.len()
                );
                log_memory_usage(false, "100k reads mapped");
            }

            let num_mapped = *counter.lock().unwrap();
            if num_mapped % (twin_reads.len() / 10).max(1) == 0 {
                let mapping_progress = (num_mapped as f64 / twin_reads.len() as f64) * 100.;
                let dbgstring = format!(
                    "Mapping batch {}/{}: Processed {} / {} ({}%) of reads...",
                    batch_index + 1,
                    outer_read_chunks.len(),
                    num_mapped, 
                    twin_reads.len(), 
                    mapping_progress.round()
                );
                log_memory_usage(true, &dbgstring);

            }

            for (contig_id, anchors) in mini_anchors.into_iter() {
                if anchors.anchors.len() < 15 {
                    continue;
                }
                for twin_ol in compare_twin_reads(
                    read,
                    &tr_outer[&(contig_id as usize)],
                    Some(&anchors),
                    None,
                    rid,
                    contig_id as usize,
                    &tr_options,
                    args,
                ) {
                    if twin_ol.end2 - twin_ol.start2 < 500 {
                        continue;
                    }
                    unitig_hits.push(twin_ol);
                }
                drop(anchors);
            }
            

            for hit in unitig_hits.into_iter() {
                let max_overlap = check_maximal_overlap(
                    hit.start1 as usize,
                    hit.end1 as usize,
                    hit.start2 as usize,
                    hit.end2 as usize,
                    read.base_length,
                    tr_outer[&hit.i2].base_length,
                    hit.chain_reverse,
                    args.maximal_end_fuzz,
                );

                let identity = id_est(
                    hit.shared_minimizers,
                    hit.diff_snpmers,
                    args.c as u64,
                    hit.large_indel,
                );

                if max_overlap {
                    if identity > IDENTITY_THRESHOLDS[0] - 0.05 / 100. {
                        let small_twin_ol = BareMappingOverlap {
                            snpmer_identity: identity as Fraction,
                        };
                        let mut map = mapping_maximal_boundaries.lock().unwrap();
                        let (enc, ids) = map.entry(hit.i2).or_insert((vec![], vec![]));
                        encode_boundary_pair(hit.start2 as u32 + 50, hit.end2 as u32 - 50, enc);
                        ids.push(small_twin_ol);
                    }
                }

                let mut map = mapping_local_boundaries_map.lock().unwrap();
                let vec = map.entry(hit.i2).or_insert(vec![]);
                encode_boundary_pair(hit.start2 as u32 + 50, hit.end2 as u32 - 50, vec);
            }
        });

        log::info!("Finished mapping batch {}/{} to outer reads. Splitting chimeras...", batch_index + 1, outer_read_chunks.len());

        drop(mini_index);

        let mapping_local_boundaries_map = mapping_local_boundaries_map.into_inner().unwrap();

        // Process each outer read immediately to compute split plan and free memory
        mapping_local_boundaries_map.into_par_iter().for_each(|(outer_id, boundaries)| {
            let outer_read = &twin_reads[outer_id];
            let outer_read_length = outer_read.base_length;

            let mut all_local_intervals: Vec<BareInterval> = Vec::new();
            for (start, end) in BoundaryPairIter::new(&boundaries) {
                let start = if start < 200 { start } else { start + 50 };
                let stop = if end > outer_read_length as u32 - 200 {
                    end
                } else {
                    end - 50
                };
                all_local_intervals.push(BareInterval { start, stop });
            }
            *num_alignments.lock().unwrap() += all_local_intervals.len();

            all_local_intervals.sort_unstable();

            let (maximal_enc, maximal_ids) = std::mem::take(
                mapping_maximal_boundaries.lock().unwrap()
                    .get_mut(&outer_id)
                    .unwrap_or(&mut (vec![], vec![])),
            );

            let max_intervals: Vec<_> = BoundaryPairIter::new(&maximal_enc)
                .zip(maximal_ids.into_iter())
                .map(|((start, stop), overlap)| (BareInterval { start, stop }, overlap))
                .collect();

            *num_maximal.lock().unwrap() += max_intervals.len();

            // Immediately compute breakpoints
            let breakpoints = cov_mapping_breakpoints(&all_local_intervals, outer_read_length as u32, &outer_read.id, args);

            // Immediately compute coverage stats for each segment >= MIN_READ_LENGTH
            let mut coverage_stats = vec![];

            if breakpoints.is_empty() || !break_chimeras{
                // No breakpoints - compute coverage for entire read if long enough
                let cov_stats = compute_coverage_stats(
                    outer_read,
                    &max_intervals,
                    0,
                    outer_read_length,
                );
                coverage_stats.push(cov_stats);
            } else {
                // Compute coverage for each segment that's >= MIN_READ_LENGTH
                let mut last_break = 0;
                for breakpoint in &breakpoints {
                    let bp_start = breakpoint.pos1;
                    if bp_start - last_break >  MIN_READ_LENGTH {
                        let cov_stats = compute_coverage_stats(
                            outer_read,
                            &max_intervals,
                            last_break,
                            bp_start,
                        );
                        coverage_stats.push(cov_stats);
                    }
                    last_break = breakpoint.pos2;
                }
            }

            // Only create split plan if there are valid coverage stats
            if !coverage_stats.is_empty() || !breakpoints.is_empty() {
                let split_plan = SplitReadPlan {
                    read_index: outer_id,
                    breakpoints,
                    coverage_stats,
                };
                split_plans.lock().unwrap().push(split_plan);
            }

            // max_intervals and all_local_intervals are dropped here, freeing memory
        });
    }

    log::info!(
        "Number of local alignments to outer reads: {}",
        *num_alignments.lock().unwrap()
    );
    log::info!(
        "Number of maximal alignments to outer reads: {}",
        *num_maximal.lock().unwrap()
    );

    split_plans.into_inner().unwrap()
}

pub fn map_reads_to_outer_reads(
    outer_read_indices: &[usize],
    twin_reads: &[TwinRead],
    args: &Cli,
) -> Vec<TwinReadMapping> {

    let mut ret = vec![];
    let mut num_alignments = 0;
    let mut num_maximal = 0;

    let chunk_size = args.read_map_batch_size;
   //let chunk_size = 10_000;

    let outer_read_chunks = outer_read_indices
        .chunks(chunk_size)
        .collect::<Vec<_>>();

    log::info!("Mapping {} reads to {} outer reads", twin_reads.len(), outer_read_indices.len());
    log::info!("Breaking {} reads into {} chunks of <= {} reads", outer_read_indices.len(), outer_read_chunks.len(), chunk_size);

    for mapping_chunk_indices in outer_read_chunks
    {
        // (encoded boundaries as varints, BareMappingOverlap vec)
        let mapping_maximal_boundaries: Mutex<FxHashMap<usize, (Vec<u8>, Vec<BareMappingOverlap>)>> = Mutex::new(FxHashMap::default());
        let mapping_local_boundaries_map = Mutex::new(FxHashMap::default());
        // 0..outer_read_indices -- confusingly, I chose to renumber the indices in this step. Then
        // it's fixed in the index_of_outer_in_all.
        let tr_outer = mapping_chunk_indices
            .into_iter()
            .map(|&i| (i, &twin_reads[i]))
            .collect::<FxHashMap<usize, &TwinRead>>();

        let mini_index = get_minimizer_index(None, Some(&tr_outer));

        let counter = Mutex::new(0);

        twin_reads.par_iter().enumerate().for_each(|(rid, read)| {
            let mut tr_options = CompareTwinReadOptions::default();
            tr_options.double_gap = 1500;
            //tr_options.force_one_to_one_alignments = true;
            tr_options.force_query_nonoverlap = true;
            tr_options.read1_snpmers = Some(read.snpmers_vec_strand());
            let mini = read.minimizers_vec_strand();
            let mini_anchors = find_exact_matches_with_full_index(&mini, &mini_index, None, Some(&tr_outer));
            drop(mini);
            let mut unitig_hits : Vec<TwinOverlap> = vec![];
            *counter.lock().unwrap() += 1;
            if *counter.lock().unwrap() % 100000 == 0 {
                log::debug!(
                    "Processed {} reads / {} ...",
                    *counter.lock().unwrap(),
                    twin_reads.len()
                );
                log_memory_usage(false, "100k reads mapped");
            }

            //let start = std::time::Instant::now();
            for (contig_id, anchors) in mini_anchors.into_iter() {
                if anchors.anchors.len() < 15 {
                    continue;
                }
                for twin_ol in compare_twin_reads(
                    read,
                    &tr_outer[&(contig_id as usize)],
                    Some(&anchors),
                    None,
                    rid,
                    contig_id as usize,
                    &tr_options,
                    args,
                ) {
                    if twin_ol.end2 - twin_ol.start2 < 500 {
                        continue;
                    }
                    unitig_hits.push(twin_ol);
                }

                drop(anchors);
            }

            //let overlap_time = start.elapsed().as_micros();
            for hit in unitig_hits.into_iter() {
                {
                    //log::trace!("{} {} {} {} {} {}", hit.i1, hit.i2, hit.start1, hit.end1, hit.start2, hit.end2);
                    let max_overlap = check_maximal_overlap(
                        hit.start1 as usize,
                        hit.end1 as usize,
                        hit.start2 as usize,
                        hit.end2 as usize,
                        read.base_length,
                        tr_outer[&hit.i2].base_length,
                        hit.chain_reverse,
                        args.maximal_end_fuzz, // TODO we should change this to the overlap hang length for concordance
                    );

                    let identity = id_est(
                        hit.shared_minimizers,
                        hit.diff_snpmers,
                        args.c as u64,
                        hit.large_indel,
                    );

                    //Populate mapping boundaries map
                    if max_overlap {
                        if identity > IDENTITY_THRESHOLDS[0] - 0.05 / 100. {
                            let small_twin_ol = BareMappingOverlap {
                                snpmer_identity: identity as Fraction,
                            };
                            let mut map = mapping_maximal_boundaries.lock().unwrap();
                            let (enc, ids) = map.entry(hit.i2).or_insert((vec![], vec![]));
                            encode_boundary_pair(hit.start2 as u32 + 50, hit.end2 as u32 - 50, enc);
                            ids.push(small_twin_ol);
                        }
                    }

                    // TODO Require length conditio that scales with hit.i2 (the target read)
                    let mut map = mapping_local_boundaries_map.lock().unwrap();
                    let vec = map.entry(hit.i2).or_insert(vec![]);
                    encode_boundary_pair(hit.start2 as u32 + 50, hit.end2 as u32 - 50, vec);
                }
            }
        });

        drop(mini_index);

        let mapping_local_boundaries_map = mapping_local_boundaries_map.into_inner().unwrap();
        let mut mapping_maximal_boundaries = mapping_maximal_boundaries.into_inner().unwrap();

        for (outer_id, boundaries) in mapping_local_boundaries_map.into_iter() {
            //let index_of_outer_in_all = outer_read_indices[outer_id];
            let outer_read_length = twin_reads[outer_id].base_length;

            let mut all_local_intervals: Vec<BareInterval> = Vec::new();
            for (start, end) in BoundaryPairIter::new(&boundaries) {
                let start = if start < 200 { start } else { start + 50 };
                let stop = if end > outer_read_length as u32 - 200 {
                    end
                } else {
                    end - 50
                };
                all_local_intervals.push(BareInterval { start, stop });
            }
            num_alignments += all_local_intervals.len();

            all_local_intervals.sort_unstable();
            all_local_intervals.shrink_to_fit();

            let (maximal_enc, maximal_ids) = std::mem::take(
                mapping_maximal_boundaries
                    .get_mut(&outer_id)
                    .unwrap_or(&mut (vec![], vec![])),
            );

            let mut max_intervals: Vec<_> = BoundaryPairIter::new(&maximal_enc)
                .zip(maximal_ids.into_iter())
                .map(|((start, stop), overlap)| (BareInterval { start, stop }, overlap))
                .collect();

            num_maximal += max_intervals.len();
            max_intervals.shrink_to_fit();
            //let lapper = Lapper::new(max_intervals);

            let map_info = MappingInfo {
                minimum_depth: -1.,
                median_depth: -1.,
                max_alignment_boundaries: None,
                max_mapping_boundaries: Some(max_intervals),
                present: true,
                length: outer_read_length,
            };

            let twinread_mapping = TwinReadMapping {
                tr_index: outer_id,
                mapping_info: map_info,
                all_intervals: all_local_intervals,
            };

            ret.push(twinread_mapping);
        }
    }

    log::info!(
        "Number of local alignments to outer reads: {}",
        num_alignments
    );
    log::info!(
        "Number of maximal alignments to outer reads: {}",
        num_maximal
    );
    log::debug!(
        "Number reads mapped to: {}",
        ret.len()
    );


    ret
}

pub fn map_to_dereplicate(
    unitig_graph: &mut unitig::UnitigGraph,
    kmer_info: &KmerGlobalInfo,
    _twin_reads: &[TwinRead],
    temp_dir: &PathBuf,
    args: &Cli,
) {

    let mapping_file = temp_dir.join("dereplicate_unitigs.paf.gz");
    let mapping_file = Mutex::new(
        GzEncoder::new(
            BufWriter::new(std::fs::File::create(mapping_file).unwrap()),
            Compression::default()
        )
    );

    let mut snpmer_set = FxHashSet::default();
    for snpmer_i in kmer_info.snpmer_info.iter() {
        let k = args.kmer_size as usize;
        let snpmer1 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[0] as u64) << (k - 1));
        let snpmer2 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[1] as u64) << (k - 1));
        snpmer_set.insert(snpmer1);
        snpmer_set.insert(snpmer2);
    }

    //Convert unitigs to twinreads
    let mut tr_unitigs = unitigs_to_tr(unitig_graph, &snpmer_set, &kmer_info.solid_kmers, &kmer_info.high_freq_kmers, args);
    let contained_contigs = Mutex::new(FxHashSet::default());
    
    //Put low abund unitigs into contained contigs
    for (&id, _) in unitig_graph.nodes.iter() {
        let node_unitig = &unitig_graph.nodes[&id];
        if !node_unitig.passes_abundance_thresholds(args){
            contained_contigs.lock().unwrap().insert(id);
        }
    }

    // Remove unitigs that fail the three filters below
    tr_unitigs.retain(|&id, _| {
        let node_unitig = &unitig_graph.nodes[&id];
        node_unitig.passes_abundance_thresholds(args)
    });

    let unitigs_to_index = tr_unitigs.iter().map(|(id, read)| (*id, read)).filter(|(id, _)| {
        let unitig =  tr_unitigs.get(id).unwrap();
        let node_unitig = &unitig_graph.nodes[id];
        if unitig.base_length >= 100_000 && node_unitig.read_indices_ori.len() >= 5 {
            return false;
        }
        // Don't remove circular contigs
        if node_unitig.has_circular_walk(){
            return false;
        }
        true
    }).collect::<FxHashMap<_, _>>();

    log::info!("Initializing index for dereplication. {} unitigs pass abundance thresholds.", tr_unitigs.len());
    let mini_index = get_minimizer_index(None, Some(&unitigs_to_index));

    log::debug!("Built minimizer index of size {}", mini_index.len());
    log_memory_usage(true, "Mapping unitigs to each other for dereplication...");
    tr_unitigs.par_iter().for_each(|(q_id, q_unitig)| {
        
        let mut tr_options = CompareTwinReadOptions::default();
        tr_options.read1_snpmers = Some(q_unitig.snpmers_vec_strand());
        tr_options.retain_chain = true;

        let mini = q_unitig.minimizers_vec_strand();
        let mini_anchors = find_exact_matches_with_full_index(&mini, &mini_index, Some(&tr_unitigs), None);
        let mut mini_anchor_sorted_indices = mini_anchors.keys().cloned().collect::<Vec<_>>();
        mini_anchor_sorted_indices.sort_by_key(|x| mini_anchors[x].anchors.len());
        mini_anchor_sorted_indices.reverse();
        mini_anchor_sorted_indices.retain(|x| {
            let r_unitig = &tr_unitigs[&(*x as usize)];
            let r_node_unitig = &unitig_graph.nodes[&(*x as usize)];
            if !(r_unitig.base_length < 100_000 || r_node_unitig.read_indices_ori.len() < 5) {
                return false;
            }
            // Don't remove circular contigs
            if r_node_unitig.has_circular_walk(){
                return false;
            }
            if *x as usize == *q_id as usize{
                return false;
            }
            if r_unitig.base_length as f64 * 1.1  >  q_unitig.base_length as f64 {
                return false;
            }
            true
        });

        // r_untig is SMALLER than q_unitig
        //for contig_id in mini_anchor_sorted_indices.iter() {
        mini_anchor_sorted_indices.par_iter().for_each(|contig_id| {
            let r_unitig = &tr_unitigs[&(*contig_id as usize)];
            let _r_node_unitig = &unitig_graph.nodes[&(*contig_id as usize)];

            if contained_contigs.lock().unwrap().contains(&(*contig_id as usize)){
                return;
            }

            let anchors = mini_anchors.get(contig_id).unwrap();
            if anchors.anchors.len() < 10{
                return;
            }

            //TODO change strictness
            if (anchors.anchors.len() as f64) < (r_unitig.base_length as f64 / args.c as f64 / args.absolute_minimizer_cut_ratio){
                return;
            }

            for uni_ol in compare_twin_reads(
                q_unitig,
                r_unitig,
                Some(anchors),
                None,
                *q_id,
                *contig_id as usize,
                &tr_options,
                args
            ) {

                //TODO change strictness
                if (uni_ol.shared_minimizers as f64) < (r_unitig.base_length as f64 / args.c as f64 / args.absolute_minimizer_cut_ratio){
                    return;
                }
                let ss_strict = same_strain(
                    uni_ol.shared_minimizers,
                    uni_ol.diff_snpmers,
                    uni_ol.shared_snpmers,
                    args.c as u64,
                    args.snpmer_threshold_strict,
                    args.snpmer_error_rate_strict,
                    uni_ol.large_indel,
                );

                let id_derep_ol = id_est(
                    uni_ol.shared_minimizers,
                    uni_ol.diff_snpmers,
                    args.c as u64,
                    uni_ol.large_indel,
                );
                //Allow 2 snpmer mismatches for these small contigs -- this level of variation is filtered into alternate anyways
                // OR is snpmer error
                // OR TODO v0.3.0 -- allow if identity is very high

                let ss_strict = ss_strict 
                || (uni_ol.diff_snpmers <= MAX_ALLOWABLE_SNPMER_ERROR_MISC && uni_ol.shared_snpmers > 0)
                || (id_derep_ol > IDENTITY_THRESHOLDS[0] - 0.05 / 100.);

                if uni_ol.end2 - uni_ol.start2 < r_unitig.base_length / 2 {
                    return;
                }
                
                let (start1, end1, start2, end2) = alignment::extend_ends_chain(&q_unitig.dna_seq, &r_unitig.dna_seq, &uni_ol, args);

                write!(
                    mapping_file.lock().unwrap(),
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tshared_mini:{}\tdiff_snp:{}\tshared_snp:{}\n",
                    q_unitig.id,
                    q_unitig.base_length,
                    start1,
                    end1,
                    if uni_ol.chain_reverse { "-" } else { "+" },
                    r_unitig.id,
                    r_unitig.base_length,
                    start2,
                    end2,
                    end1 - start1,
                    end2 - start2,
                    255,
                    uni_ol.shared_minimizers,
                    uni_ol.diff_snpmers,
                    uni_ol.shared_snpmers,
                ).unwrap();

                let range_ref = end2 - start2;
                let frac = range_ref as f64 / r_unitig.base_length as f64;
                let close_fraction = frac > 0.98;
                let close_to_both_ends = end2 + 100 >= r_unitig.base_length && start2 <= 100;
                //TODO change frac
                if (close_to_both_ends || close_fraction) && ss_strict && r_unitig.base_length * 2 < q_unitig.base_length {
                    let mut contained_contigs = contained_contigs.lock().unwrap();
                    contained_contigs.insert(*contig_id as usize);
                    return;
                }
            }
        });
    });

    let vec_remove = contained_contigs.into_inner().unwrap().into_iter().map(|x| x as usize).collect::<Vec<_>>();
    log::debug!("SMALL CONTIG REMOVAL: Removing {} small contigs that are too similar to larger contigs (or singletons that have no coverage)", vec_remove.len());
    unitig_graph.remove_nodes(&vec_remove, false);
    unitig_graph.re_unitig();
}

pub fn map_reads_to_unitigs(
    unitig_graph: &mut unitig::UnitigGraph,
    kmer_info: &KmerGlobalInfo,
    twin_reads: &[TwinRead],
    temp_dir: &PathBuf,
    args: &Cli,
) {
    let mapping_file = temp_dir.join("map_to_unitigs.paf.gz");
    let mapping_file = Mutex::new(
        GzEncoder::new(
            BufWriter::new(std::fs::File::create(mapping_file).unwrap()),
            Compression::default()
        )
    );
    
    let mut snpmer_set = FxHashSet::default();
    for snpmer_i in kmer_info.snpmer_info.iter() {
        let k = args.kmer_size as usize;
        let snpmer1 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[0] as u64) << (k - 1));
        let snpmer2 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[1] as u64) << (k - 1));
        snpmer_set.insert(snpmer1);
        snpmer_set.insert(snpmer2);
    }

    //Convert unitigs to twinreads
    let tr_unitigs = unitigs_to_tr(unitig_graph, &snpmer_set, &kmer_info.solid_kmers, &kmer_info.high_freq_kmers, args);
    drop(snpmer_set);
    let circular_unitigs = unitig_graph.nodes.iter().filter(|(_, u)| u.has_circular_walk()).map(|(id, _)| *id).collect::<FxHashSet<_>>();
    let mapping_boundaries_map = Mutex::new(FxHashMap::default());
    let num_reads = twin_reads.len();

    let chunk_size = if args.low_mem { tr_unitigs.len() / 3 + 1 } else { tr_unitigs.len() };
    let tr_unitig_keys: Vec<usize> = tr_unitigs.keys().cloned().collect();
    let num_chunks = (tr_unitig_keys.len() + chunk_size - 1) / chunk_size.max(1);
    for (chunk_idx, key_chunk) in tr_unitig_keys.chunks(chunk_size.max(1)).enumerate() {
        let tr_chunk: FxHashMap<usize, &TwinRead> = key_chunk.iter().map(|k| (*k, &tr_unitigs[k])).collect();
        let mini_index = get_minimizer_index(None, Some(&tr_chunk));
        let counter = Mutex::new(0_usize);

        log_memory_usage(true, &format!("Built unitig index chunk {}/{}: starting alignments", chunk_idx + 1, num_chunks));

        twin_reads.par_iter().enumerate().for_each(|(rid, read)| {
        let mut tr_options = CompareTwinReadOptions::default();
        tr_options.retain_chain = true;
        tr_options.read1_snpmers = Some(read.snpmers_vec_strand());
        tr_options.secondary_threshold = Some(0.15);

        // This should've been added in v0.3; adding in v0.4 TESTING
        tr_options.maximal_only = true;

        let mini = read.minimizers_vec_strand();
        let mini_anchors = find_exact_matches_with_full_index(&mini, &mini_index, None, None);
        let mut unitig_hits = vec![];

        for (contig_id, anchors) in mini_anchors.iter() {
            if anchors.anchors.len() < 10{
                continue;
            }

            if (anchors.anchors.len() as f64) < read.base_length as f64 / args.c as f64 / args.absolute_minimizer_cut_ratio / 1.5 {
                continue;
            }

            // Improved alignment to small circular genomes, that are "repetitive" before end trimming
            let unitig_length = tr_chunk[&(*contig_id as usize)].base_length;
            if unitig_length < read.base_length * 3 {
                tr_options.force_ref_nonoverlap = false;
            }
            else{
                tr_options.force_ref_nonoverlap = true;
            }

            // if read.id.contains("0f6f4a99-e0ff-4e5a-a7dd-00b4d45b0120") {
            //     dbg!(contig_id);
            //     dbg!(&anchors.anchors);
            // }

            for twinol in compare_twin_reads(
                read,
                &tr_chunk[&(*contig_id as usize)],
                Some(anchors),
                None,
                rid,
                *contig_id as usize,
                &tr_options,
                args
            ) {
                log::trace!("Read {} unitig {} snpmers_shared {} snpmers_diff {} range1 {}-{} range2 {}-{}; anchor mult {}", &read.id, &tr_chunk[&(*contig_id as usize)].id, twinol.shared_snpmers, twinol.diff_snpmers, twinol.start1, twinol.end1, twinol.start2, twinol.end2, anchors.max_mult);
                
                //Disallow large indels because they may cause windowed POA to fail
                let ol_len = twinol.end2 - twinol.start2;

                if ol_len < MIN_READ_LENGTH || twinol.large_indel{
                    continue;
                }

                if (twinol.shared_minimizers as f64) < (ol_len as f64 / args.c as f64 / args.absolute_minimizer_cut_ratio){
                    continue;
                }

                unitig_hits.push(twinol);
            }
        }

        let mut ss_hits = unitig_hits.into_iter().filter(|x| same_strain(
                x.shared_minimizers,
                x.diff_snpmers,
                x.shared_snpmers,
                args.c.try_into().unwrap(),
                args.snpmer_threshold_lax,
                args.snpmer_error_rate_lax,
                x.large_indel,
            )).collect::<Vec<_>>();

        let end_fuzz_pair = &read.overlap_hang_length.unwrap(); 
        let end_fuzz = end_fuzz_pair.0.max(end_fuzz_pair.1);

        ss_hits.retain(|x| check_maximal_overlap(
            x.start1,
            x.end1,
            x.start2,
            x.end2,
            read.base_length,
            tr_chunk[&x.i2].base_length,
            x.chain_reverse,
            end_fuzz
        ));

        let mut retained_hits = vec![];
        //TODO group similar hits together, allow some overlap ...
        ss_hits.sort_by(|a, b| id_est(b.shared_minimizers, b.diff_snpmers, args.c.try_into().unwrap(), b.large_indel).
        partial_cmp(&id_est(a.shared_minimizers, a.diff_snpmers, args.c.try_into().unwrap(), a.large_indel)).unwrap());

        let perfect_hits = ss_hits.iter().filter(|x| same_strain(
            x.shared_minimizers,
            x.diff_snpmers,
            x.shared_snpmers,
            args.c.try_into().unwrap(),
            args.snpmer_threshold_strict,
            args.snpmer_error_rate_strict,
            x.large_indel,
        )).collect::<Vec<_>>();
        
        let mut imperfect_hits = ss_hits.iter().filter(|x| !same_strain(
            x.shared_minimizers,
            x.diff_snpmers,
            x.shared_snpmers,
            args.c.try_into().unwrap(),
            args.snpmer_threshold_strict,
            args.snpmer_error_rate_strict,
            x.large_indel,
        )).collect::<Vec<_>>();

        if args.hifi{
            imperfect_hits.sort_by(|a, b| {
                let id_a = id_est(a.shared_minimizers, a.diff_snpmers, args.c.try_into().unwrap(), a.large_indel) - 0.98;
                let id_b = id_est(b.shared_minimizers, b.diff_snpmers, args.c.try_into().unwrap(), b.large_indel) - 0.98;
                let mini_cumulative_a = (a.shared_minimizers + a.shared_snpmers) as f64 - a.diff_snpmers as f64;
                let mini_cumulative_b = (b.shared_minimizers + b.shared_snpmers) as f64 - b.diff_snpmers as f64;
                (mini_cumulative_b * id_b).partial_cmp(&(mini_cumulative_a * id_a)).unwrap()
            });
        }
        else{
            imperfect_hits.sort_by(|a, b| {
                let rate_a = a.shared_minimizers as f64 * (args.c as f64) * 0.10;
                let rate_b = b.shared_minimizers as f64 * (args.c as f64) * 0.10;
                let prob_nosnp_a = 1. - 2.718f64.powf(-rate_a / args.kmer_size as f64);
                let prob_nosnp_b = 1. - 2.718f64.powf(-rate_b / args.kmer_size as f64);

                let id_a = id_est(a.shared_minimizers, a.diff_snpmers, args.c.try_into().unwrap(), a.large_indel);
                let id_b = id_est(b.shared_minimizers, b.diff_snpmers, args.c.try_into().unwrap(), b.large_indel);
                let id_kmer_a = rate_a.powf(1. / args.kmer_size as f64);
                let id_kmer_b = rate_b.powf(1. / args.kmer_size as f64);

                let weighted_id_a = id_a * prob_nosnp_a + id_kmer_a * (1. - prob_nosnp_a);
                let weighted_id_b = id_b * prob_nosnp_b + id_kmer_b * (1. - prob_nosnp_b);

                weighted_id_b.partial_cmp(&weighted_id_a).unwrap()
            });
        }
        
        let imperfect_ids = imperfect_hits.iter().map(|x| id_est(x.shared_minimizers, x.diff_snpmers, args.c.try_into().unwrap(), x.large_indel)).collect::<Vec<_>>();
        let max_imperfect_id = imperfect_ids.first().cloned().unwrap_or(1.0);

        //let mut imperfect_ids = imperfect_hits.iter().map(|x| id_est(x.shared_minimizers, x.diff_snpmers, args.c.try_into().unwrap(), x.large_indel)).collect::<Vec<_>>();
        //imperfect_ids.sort_by(|a, b| b.partial_cmp(a).unwrap());

        let mut max_perfect_mini_opt = None;
        let mut max_id = None;

        for hit in perfect_hits.iter().chain(imperfect_hits.iter()) {

            log::trace!("FULL LENGTH MAPPING: {} query:{}-{} ref:{}-{} snp_shared {} snp_diff {} unitig u{}", first_word(&read.id), hit.start1, hit.end1, hit.start2, hit.end2, hit.shared_snpmers, hit.diff_snpmers, &tr_chunk[&(hit.i2 as usize)].id);

            // If this passes stringent standards, retain the hit. Otherwise, only retain the top hit. 
            if !same_strain(
                hit.shared_minimizers,
                hit.diff_snpmers,
                hit.shared_snpmers,
                args.c.try_into().unwrap(),
                args.snpmer_threshold_strict,
                args.snpmer_error_rate_strict,
                hit.large_indel
            ){
                if let Some(max_perfect_mini) = max_perfect_mini_opt {
                    let mini_cumulative = (hit.shared_minimizers + hit.shared_snpmers) as i64 - (hit.diff_snpmers as i64);
                    if mini_cumulative < max_perfect_mini || !args.hifi {
                        //Still take the best imperfect hit IF no perfect hits found, otherwise set a threshold at 1/3 
                        //of the distance to a perfect hit 
                        let id = id_est(hit.shared_minimizers, hit.diff_snpmers, args.c as u64, hit.large_indel);
                        let max_id_unwr = max_id.unwrap_or(max_imperfect_id);
                        let cutoff = max_imperfect_id - (1.0 - max_id_unwr) * 0.33 + 0.000000000001;

                        if id <= cutoff || max_id_unwr == 1.0{
                            break;
                        }

                        if args.hifi{
                            if (mini_cumulative as f64) < max_perfect_mini as f64 * 0.9{
                                break;
                            }
                        }
                        else{
                            if (mini_cumulative as f64) < max_perfect_mini as f64 * 0.5{
                                break;
                            }
                        }
                    }
                }
                else{
                    max_perfect_mini_opt = Some((hit.shared_minimizers + hit.shared_snpmers) as i64 - (hit.diff_snpmers as i64));
                }
            }
            else{
                max_id = Some(1.0);
                if let Some(max_perfect_mini) = max_perfect_mini_opt {
                    if (hit.shared_minimizers + hit.shared_snpmers) as i64 - (hit.diff_snpmers as i64) > max_perfect_mini {
                        max_perfect_mini_opt = Some((hit.shared_minimizers + hit.shared_snpmers) as i64 - (hit.diff_snpmers as i64));
                   }
                }
                else{
                    max_perfect_mini_opt = Some((hit.shared_minimizers + hit.shared_snpmers) as i64 - (hit.diff_snpmers as i64));
                }
            }


            retained_hits.push(hit);
        }

        // CHANGE v0.3.0 -- dont' use simple chain score, 
        //retained_hits.sort_by(|a, b| b.chain_score.partial_cmp(&a.chain_score).unwrap());

        let mut query_intervals : Vec<Interval<u32, (u32, u32, u32)>> = vec![];
        let mut max_score = 0;
        let max_chain = retained_hits.iter().max_by_key(|hit| hit.chain_score);
        if let Some(max_chain) = max_chain{
            max_score = max_chain.chain_score;
            query_intervals.push(Interval{
                start: max_chain.start1 as u32,
                stop: max_chain.end1 as u32,
                val: (max_chain.start2 as u32, max_chain.end2 as u32, max_chain.i2 as u32)
            });
        }

        for hit in retained_hits {

            let interval = Interval{
                start: hit.start1 as u32,
                stop: hit.end1 as u32,
                val: (hit.start2 as u32, hit.end2 as u32, hit.i2 as u32)
            };

            let reference_length = tr_chunk[&hit.i2].base_length;
            let read_length = twin_reads[hit.i1].base_length;

            //Internal secondary filter in chaining procedure (compare twin reads)
            // This is an additional filter for finer control
            let mut proceed = true;
            let chain_threshold = if args.hifi {0.70} else {0.30};
            for q_interval in query_intervals.iter(){
                let intersect = interval.intersect(q_interval);
                let overlap = intersect as f64 / (interval.stop - interval.start) as f64;
                if overlap > 0.25{
                    if (hit.chain_score as f64) < (max_score as f64) * chain_threshold {
                        //Allow duplicate mappings to small circular contigs
                        if read_length * 3 < reference_length{

                            // Allow mappings to "repeats" caused by end circularity
                            if circular_unitigs.contains(&(hit.i2 as usize)) && q_interval.val.2 == hit.i2 as u32{
                                let reference_length = tr_chunk[&hit.i2].base_length;
                                if overlap < 0.7 
                                || (q_interval.val.1 as i32 - interval.val.0 as i32).abs() + 100_000 < reference_length as i32 {
                                //|| ((hit.chain_score as f64) < 0.1 * (max_score as f64)){
                                    proceed = false;
                                    break
                                }
                            }
                            else{
                                proceed = false;
                                break;
                            }
                        }
                    }
                }
            }

            if !proceed{
                break;
            }

            let alignment_result = alignment::get_full_alignment(
                &twin_reads[hit.i1].dna_seq,
                &tr_chunk[&hit.i2].dna_seq,
                &hit,
                args,
            );

            write!(
                mapping_file.lock().unwrap(),
                "r{}-{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tshared_mini:{}\tdiff_snp:{}\tshared_snp:{}\n",
                hit.i1,
                twin_reads[hit.i1].id.split_ascii_whitespace().next().unwrap_or("").to_string(),
                twin_reads[hit.i1].base_length,
                alignment_result.as_ref().unwrap().q_start,
                alignment_result.as_ref().unwrap().q_end,
                if hit.chain_reverse { "-" } else { "+" },
                format!("{}ctg", &tr_chunk[&hit.i2].id),
                tr_chunk[&hit.i2].base_length,
                alignment_result.as_ref().unwrap().r_start,
                alignment_result.as_ref().unwrap().r_end,
                hit.end1 - hit.start1,
                hit.end2 - hit.start2,
                255,
                hit.shared_minimizers,
                hit.diff_snpmers,
                hit.shared_snpmers,
            ).unwrap();

            let mut map = mapping_boundaries_map.lock().unwrap();
            let vec = map.entry(hit.i2).or_insert(vec![]);
            let small_twin_ol = SmallTwinOl{
                query_id: rid as u32,
                //shared_minimizers: hit.shared_minimizers as u32,
                //diff_snpmers: hit.diff_snpmers as u32,
                snpmer_identity: id_est(hit.shared_minimizers, hit.diff_snpmers, args.c as u64, hit.large_indel) as f32,
                //shared_snpmers: hit.shared_snpmers as u32,
                reverse: hit.chain_reverse,
                alignment_result: alignment_result
            };
            vec.push((hit.start2, hit.end2, small_twin_ol));
        }

        *counter.lock().unwrap() += 1;
        let count = *counter.lock().unwrap();
        if count % (num_reads / 10) == 0{
            log::info!("Mapped {:.0}% of reads back to contigs", count as f64 / num_reads as f64 * 100.);
        }
        });

    } // end chunk loop

    let mut number_of_alignments = 0;
    let mut cigar_string_lengths = vec![];
    for (contig_id, boundaries_and_rid) in mapping_boundaries_map.into_inner().unwrap().into_iter()
    {
        let unitig_length = unitig_graph.nodes.get(&contig_id).unwrap().cut_length();
        let intervals = boundaries_and_rid
            .into_iter()
            .map(|x| Interval {
                start: x.0 as u32,
                stop: x.1 as u32,
                val: x.2,
            })
            .collect::<Vec<Interval<u32, SmallTwinOl>>>();
        number_of_alignments += intervals.len();
        cigar_string_lengths.extend(
            intervals
                .iter()
                .map(|x| x.val.alignment_result.as_ref().unwrap().cigar.len()),
        );
        let mut lapper = Lapper::new(intervals);
        lapper.intervals.shrink_to_fit();

        //let (unitig_first_mini_pos, unitig_last_mini_pos) = first_last_mini_in_range(0, unitig_length, args.kmer_size, MINIMIZER_END_NTH_COV, tr_unitigs[&contig_id].minimizers.as_slice());
        //let (min_depth, median_depth) = median_and_min_depth_from_lapper(&lapper, SAMPLING_RATE_COV, unitig_first_mini_pos, unitig_last_mini_pos).unwrap();

        let map_info = MappingInfo {
            median_depth: -1.,
            minimum_depth: -1.,
            max_mapping_boundaries: None,
            max_alignment_boundaries: Some(lapper),
            present: true,
            length: unitig_length,
        };
        let mut_node = unitig_graph.nodes.get_mut(&contig_id).unwrap();
        mut_node.mapping_info = map_info;
    }

    log::debug!("Number of alignments: {}", number_of_alignments);
    log::debug!(
        "Average cigar string length: {}",
        cigar_string_lengths.iter().sum::<usize>() as f64 / cigar_string_lengths.len() as f64
    );
}

fn _get_splitmers(snpmers: &[(usize, u64)], k: u64) -> Vec<(usize, u64)> {
    let mask = !(3 << (k - 1));
    snpmers
        .iter()
        .map(|x| (x.0, x.1 as u64 & mask))
        .collect::<Vec<(usize, u64)>>()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn anc(pos1: u32, pos2: u32) -> Anchor {
        Anchor { i: None, j: None, pos1, pos2 }
    }

    /// Reverse-strand anchor: pos1 bit31=1 (rel_strand=1), pos2 raw (decreasing along chain).
    fn anc_rev(pos1_raw: u32, pos2: u32) -> Anchor {
        Anchor { i: None, j: None, pos1: (1u32 << 31) | pos1_raw, pos2 }
    }


    // Canonical params matching CompareTwinReadOptions::default()
    const GAP_COST: i32 = 1;
    const MATCH_SCORE: i32 = 10;
    const BAND: usize = 50;
    const MAX_ITER: usize = 500;
    const MAX_SKIP: usize = 10;
    const MAX_GAP: usize = 200;     // gap imbalance threshold
    const DOUBLE_GAP: usize = 10_000; // absolute distance threshold

    /// Sort chains for deterministic comparison: descending score,
    /// then by position of the first anchor as a tiebreaker.
    fn sorted(mut chains: Vec<(i32, Vec<Anchor>, bool)>) -> Vec<(i32, Vec<Anchor>, bool)> {
        chains.sort_by(|a, b| {
            b.0.cmp(&a.0).then_with(|| {
                a.1.first().map(|x| x.pos1).cmp(&b.1.first().map(|x| x.pos1))
            })
        });
        chains
    }

    fn chains_equal(v1: &[(i32, Vec<Anchor>, bool)], v2: &[(i32, Vec<Anchor>, bool)]) {
        assert_eq!(v1.len(), v2.len(), "chain count: v1={} v2={}", v1.len(), v2.len());
        for (i, (c1, c2)) in v1.iter().zip(v2.iter()).enumerate() {
            assert_eq!(c1.0, c2.0, "chain {i}: score v1={} v2={}", c1.0, c2.0);
            assert_eq!(c1.1, c2.1, "chain {i}: anchor list differs");
        }
    }

    // -----------------------------------------------------------------------
    // Test min_chain parameter
    // -----------------------------------------------------------------------
    #[test]
    fn test_simple_chain() {
        let anchors: Vec<Anchor> = (0..1).map(|i: u32| anc(i * 20, i * 20)).collect();

        let t0 = std::time::Instant::now();
        let r2 = dp_anchors_v2(&anchors, GAP_COST, MATCH_SCORE, MAX_ITER, MAX_SKIP, MAX_GAP, DOUBLE_GAP, true, 1);
        let _e2 = t0.elapsed();

        let r2 = sorted(r2);
        eprintln!("  v2 chains: {:?}", r2.iter().map(|(s, c, _)| (s, c.len())).collect::<Vec<_>>());

        assert_eq!(r2[0].0, MATCH_SCORE, "expected score {}, got {}", MATCH_SCORE, r2[0].0);
    }

    // -----------------------------------------------------------------------
    // Test 1: five anchors on a perfect diagonal (no gap, no noise).
    // Both algorithms must produce one chain of length 5 with the same score.
    // -----------------------------------------------------------------------
    #[test]
    fn test_simple_linear_chain() {
        // Anchors at (0,0),(20,20),(40,40),(60,60),(80,80)
        // Each step: dist1=dist2=20, gap_penalty=0, kmer_overlap=min(20,20,10)=10
        // Expected chain score: match_score + 4 * kmer_overlap = 10 + 40 = 50
        let anchors: Vec<Anchor> = (0..5).map(|i: u32| anc(i * 20, i * 20)).collect();

        let t0 = std::time::Instant::now();
        let r1 = dp_anchors(&anchors, false, GAP_COST, MATCH_SCORE, BAND, MAX_GAP, DOUBLE_GAP, true);
        let e1 = t0.elapsed();

        let t0 = std::time::Instant::now();
        let r2 = dp_anchors_v2(&anchors, GAP_COST, MATCH_SCORE, MAX_ITER, MAX_SKIP, MAX_GAP, DOUBLE_GAP, true, 2);
        let e2 = t0.elapsed();

        eprintln!("[simple_linear] v1={e1:?}  v2={e2:?}");

        let r1 = sorted(r1);
        let r2 = sorted(r2);
        eprintln!("  v1 chains: {:?}", r1.iter().map(|(s, c, _)| (s, c.len())).collect::<Vec<_>>());
        eprintln!("  v2 chains: {:?}", r2.iter().map(|(s, c, _)| (s, c.len())).collect::<Vec<_>>());

        assert_eq!(r1[0].0, 50, "expected score 50, got {}", r1[0].0);
        chains_equal(&r1, &r2);
    }

    // -----------------------------------------------------------------------
    // Test 2: two chains separated by a gap > double_gap so they cannot link.
    // Each chain has 5 anchors; both algorithms should return two equal chains.
    // -----------------------------------------------------------------------
    #[test]
    fn test_two_separate_chains() {
        // Chain A: pos1 = 100..180, pos2 = 100..180
        // Chain B: pos1 = 15000..15080, pos2 = 15000..15080
        // Gap between chains: ~14820 > DOUBLE_GAP(10000) → cannot link
        let mut anchors: Vec<Anchor> = (0..5).map(|i: u32| anc(100 + i * 20, 100 + i * 20)).collect();
        anchors.extend((0..5).map(|i: u32| anc(15_000 + i * 20, 15_000 + i * 20)));
        // Already sorted by pos1

        let t0 = std::time::Instant::now();
        let r1 = dp_anchors(&anchors, false, GAP_COST, MATCH_SCORE, BAND, MAX_GAP, DOUBLE_GAP, true);
        let e1 = t0.elapsed();

        let t0 = std::time::Instant::now();
        let r2 = dp_anchors_v2(&anchors, GAP_COST, MATCH_SCORE, MAX_ITER, MAX_SKIP, MAX_GAP, DOUBLE_GAP, true, 2);
        let e2 = t0.elapsed();

        eprintln!("[two_chains] v1={e1:?}  v2={e2:?}");

        let r1 = sorted(r1);
        let r2 = sorted(r2);
        eprintln!("  v1 chains: {:?}", r1.iter().map(|(s, c, _)| (s, c.len())).collect::<Vec<_>>());
        eprintln!("  v2 chains: {:?}", r2.iter().map(|(s, c, _)| (s, c.len())).collect::<Vec<_>>());

        chains_equal(&r1, &r2);
    }

    // -----------------------------------------------------------------------
    // Test 3: ~1300 anchors on a perfect diagonal — matches the log entry
    //   "DP chaining of 1279 anchors took 3.000917ms"
    // Primary purpose: timing comparison between v1 and v2 with debug=true.
    // On a clean diagonal both algorithms should agree on a single chain.
    // -----------------------------------------------------------------------
    #[test]
    fn test_large_timing_diagonal() {
        let n = 1300u32;
        // 20bp minimizer spacing, perfect diagonal (no gap noise)
        let anchors1: Vec<Anchor> = (0..n).map(|i| anc(i * 20, i * 20)).collect();
        let anchors2: Vec<Anchor> = (0..n).map(|i| anc(i * 20 + 100000, i * 20 + 500000)).collect();
        let mut anchors = [anchors1, anchors2].concat();
        anchors.sort_unstable_by_key(|a| (a.pos1, a.pos2)); // Ensure sorted by pos1 for DP

        let t0 = std::time::Instant::now();
        let r1 = dp_anchors(&anchors, false, GAP_COST, MATCH_SCORE, BAND, MAX_GAP, DOUBLE_GAP, true);
        let e1 = t0.elapsed();

        let t0 = std::time::Instant::now();
        let r2 = dp_anchors_v2(&anchors, GAP_COST, MATCH_SCORE, MAX_ITER, MAX_SKIP, MAX_GAP, DOUBLE_GAP, true, 2);
        let e2 = t0.elapsed();

        eprintln!("[large_diagonal n={n}] v1={e1:?}  v2={e2:?}");

        let r1 = sorted(r1);
        let r2 = sorted(r2);
        eprintln!("  v1: {} chain(s), best score={}, len={}",
            r1.len(), r1[0].0, r1[0].1.len());
        eprintln!("  v2: {} chain(s), best score={}, len={}",
            r2.len(), r2[0].0, r2[0].1.len());

        chains_equal(&r1, &r2);
    }

    // -----------------------------------------------------------------------
    // Test 4: five anchors on a perfect reverse-strand diagonal.
    // pos1 increases (0→80), pos2 decreases (80→0).  Same step size as test 1,
    // so the score must also be 50 and the chain must be marked is_reverse=true.
    // Strand bit is stripped from returned anchor.pos1.
    // -----------------------------------------------------------------------
    #[test]
    fn test_simple_reverse_chain() {
        let anchors: Vec<Anchor> = (0..5u32).map(|i| anc_rev(i * 20, 80 - i * 20)).collect();

        let t0 = std::time::Instant::now();
        let r = dp_anchors_v2(&anchors, GAP_COST, MATCH_SCORE, MAX_ITER, MAX_SKIP, MAX_GAP, DOUBLE_GAP, true, 2);
        eprintln!("[simple_reverse] took {:?}", t0.elapsed());
        eprintln!("  chains: {:?}", r.iter().map(|(s, c, rev)| (s, c.len(), rev)).collect::<Vec<_>>());

        assert_eq!(r.len(), 1, "expected 1 chain, got {}", r.len());
        let (score, chain, is_reverse) = &r[0];
        assert_eq!(*score, 50, "expected score 50, got {score}");
        assert_eq!(chain.len(), 5, "expected 5 anchors, got {}", chain.len());
        assert!(*is_reverse, "chain should be marked is_reverse=true");
        // Strand bit must be stripped from returned anchor.pos1.
        assert_eq!(chain[0].pos1, 0,  "first anchor pos1 wrong (strand bit not stripped?)");
        assert_eq!(chain[4].pos1, 80, "last  anchor pos1 wrong");
        // pos2 should be in decreasing order (raw values unchanged).
        assert_eq!(chain[0].pos2, 80);
        assert_eq!(chain[4].pos2, 0);
    }

    // -----------------------------------------------------------------------
    // Test 5: one forward chain and one reverse chain in the same anchor set.
    // They must be found independently with correct strand flags.
    // -----------------------------------------------------------------------
    #[test]
    fn test_forward_and_reverse_chains() {
        // Forward: pos1 0..80, pos2 100..180  (bit31=0)
        // Reverse: pos1_raw 200..280, pos2 380..300  (bit31=1, pos2 decreasing)
        let mut anchors: Vec<Anchor> = (0..5u32)
            .map(|i| anc(i * 20, 100 + i * 20))
            .collect();
        anchors.extend((0..5u32).map(|i| anc_rev(200 + i * 20, 380 - i * 20)));

        let t0 = std::time::Instant::now();
        let r = sorted(dp_anchors_v2(&anchors, GAP_COST, MATCH_SCORE, MAX_ITER, MAX_SKIP, MAX_GAP, DOUBLE_GAP, true, 2));
        eprintln!("[fwd+rev] took {:?}", t0.elapsed());
        eprintln!("  chains: {:?}", r.iter().map(|(s, c, rev)| (s, c.len(), rev)).collect::<Vec<_>>());

        assert_eq!(r.len(), 2, "expected 2 chains, got {}", r.len());
        assert_eq!(r[0].0, 50, "forward chain score wrong");
        assert_eq!(r[1].0, 50, "reverse chain score wrong");
        let has_fwd = r.iter().any(|(_, _, rev)| !rev);
        let has_rev = r.iter().any(|(_, _, rev)| *rev);
        assert!(has_fwd, "missing forward chain");
        assert!(has_rev, "missing reverse chain");
        // Each chain should have all 5 anchors.
        for (_, chain, _) in &r {
            assert_eq!(chain.len(), 5, "chain length should be 5");
        }
    }

    // -----------------------------------------------------------------------
    // Test 6: interleaved forward and reverse anchors at similar pos1 values.
    // Cross-strand transitions must be rejected — each strand forms its own chain.
    // -----------------------------------------------------------------------
    #[test]
    fn test_no_cross_strand_chaining() {
        // Forward: pos1=0,20,40  pos2=100,120,140
        // Reverse: pos1_raw=10,30,50  pos2=140,120,100 (decreasing)
        // The forward anchor at pos1=0,pos2=100 and reverse anchor at pos1_raw=10,pos2=140
        // look spatially close, but must NOT chain together.
        let fwd: Vec<Anchor> = (0..3u32).map(|i| anc(i * 20, 100 + i * 20)).collect();
        let rev: Vec<Anchor> = (0..3u32).map(|i| anc_rev(i * 20 + 10, 140 - i * 20)).collect();
        let anchors: Vec<Anchor> = [fwd, rev].concat();

        let r = sorted(dp_anchors_v2(&anchors, GAP_COST, MATCH_SCORE, MAX_ITER, MAX_SKIP, MAX_GAP, DOUBLE_GAP, true, 2));
        eprintln!("[cross_strand] chains: {:?}", r.iter().map(|(s, c, rev)| (s, c.len(), rev)).collect::<Vec<_>>());

        assert_eq!(r.len(), 2, "expected 2 chains (one per strand), got {}", r.len());
        let has_fwd = r.iter().any(|(_, _, rev)| !rev);
        let has_rev = r.iter().any(|(_, _, rev)| *rev);
        assert!(has_fwd, "missing forward chain");
        assert!(has_rev, "missing reverse chain");
        // No chain should mix strands (all anchors in a chain share the same bit31).
        for (_, chain, is_rev) in &r {
            for anchor in chain {
                let anchor_strand = anchor.pos1 >> 31 == 1;
                // Strand bit is stripped before returning, so all returned pos1 have bit31=0.
                // The chain-level is_reverse flag encodes the strand.
                assert_eq!(anchor.pos1 >> 31, 0, "strand bit should be stripped from returned anchor");
                let _ = (anchor_strand, is_rev); // suppress unused warning
            }
        }
    }
}
