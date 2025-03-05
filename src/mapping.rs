use crate::cli::Cli;
use crate::constants::IDENTITY_THRESHOLDS;
use crate::constants::ID_THRESHOLD_ITERS;
use rust_lapper::Interval;
use crate::constants::MAX_GAP_CHAINING;
use crate::constants::MIN_CHAIN_SCORE_COMPARE;
use crate::constants::MIN_READ_LENGTH;
use std::io::Write;
use std::io::BufWriter;
use std::path::Path;
use crate::seeding;
use crate::twin_graph::same_strain;
use std::collections::HashSet;
use crate::types::*;
use crate::unitig;
use crate::unitig::NodeSequence;
use crate::polishing::alignment;
use bio_seq::codec::Codec;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use rayon::prelude::*;
use rust_lapper::Lapper;
use std::sync::Mutex;
use crate::map_processing::*;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct HitInfo {
    pub index: u32,
    pub contig_id: u32,
    pub pos: u32,
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

pub fn find_exact_matches_with_full_index(
    seq1: &[(u32, u64)],
    index: &FxHashMap<u64, Vec<HitInfo>>,
) -> FxHashMap<u32, Anchors> {
    let mut max_mult = 0;
    let mut matches = FxHashMap::default();

    for (i, (pos, s)) in seq1.iter().enumerate() {
        if let Some(indices) = index.get(s) {
            if indices.len() > max_mult {
                max_mult = indices.len();
            }
            for hit in indices {
                let contig = hit.contig_id;
                let anchor = Anchor {
                    i: i as u32, 
                    j: hit.index as u32,
                    pos1: *pos as u32,
                    pos2: hit.pos as u32,
                };

                matches.entry(contig).or_insert(vec![]).push(anchor);
            }
        }
    }

    matches
        .into_iter()
        .map(|(k, v)| {
            (
                k,
                Anchors {
                    anchors: v,
                    max_mult,
                },
            )
        })
        .collect()
}

fn find_exact_matches_indexes_references(
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
                    i: i as u32,
                    j: *j as u32,
                    pos1: *pos as u32,
                    pos2: **pos2 as u32,
                };
                matches.push(anchor);
            }
        }
    }

    (matches, max_mult)
}

fn find_exact_matches_indexes<T>(
    seq1: T,
    seq2: T,
) -> (Vec<Anchor>, usize) 
    where T: Iterator<Item = (u32, u64)>
{
    let mut max_mult = 0;
    let mut matches = Vec::new();
    let mut index_map = FxHashMap::default();

    //Sorted
    for (j, (pos, x)) in seq2.enumerate() {
        index_map.entry(x).or_insert(vec![]).push((j, pos));
    }

    //Sorted
    for (i, (pos, s)) in seq1.enumerate() {
        if let Some(indices) = index_map.get(&s) {
            if indices.len() > max_mult {
                max_mult = indices.len();
            }
            for (j, pos2) in indices {
                let anchor = Anchor {
                    i: i as u32,
                    j: *j as u32,
                    pos1: pos as u32,
                    pos2: *pos2 as u32,
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
                (start1 as i32 - (end1) as i32).abs() 
                - (start2 as i32 - (end2) as i32).abs();
            let gap_penalty = gap_penalty_signed.abs();

            let large_indel = false;
            if gap_penalty > MAX_GAP_CHAINING as i32 {
                //large_indel = true;
                continue;
            }

            let score = dp[j] + match_score - gap_cost * gap_penalty;
            if score > dp[i] {
                dp[i] = score;
                prev[i] = Some(j);
                if score > max_score {
                    max_score = score;
                    max_index = i;
                }
                if large_indel{
                    large_indel_tracker[i] = true;
                }
            }
            else{
                unimproved += 1;
            }
            if unimproved > max_skip{
                break;
            }
        }
    }

    // Reconstruct the chain
    let mut chains = Vec::new();
    let mut used_anchors = FxHashSet::default();
    let mut reference_ranges : Vec<Interval<u32,bool>> = Vec::new();
    let mut best_indices_ordered = (0..matches.len())
        .map(|i| (dp[i], i))
        .collect::<Vec<_>>();
    best_indices_ordered.sort_unstable_by_key(|&(score, _)| -score);
    assert!(dp[max_index] == best_indices_ordered[0].0);

    for (score, best_index) in best_indices_ordered{
        if used_anchors.contains(&best_index) {
            continue;
        }
        if score < max_score * 1 / 4 {
            break;
        }

        let mut chain = Vec::new();
        let mut large_indel = false;
        let mut i = Some(best_index);
        while let Some(idx) = i {
            large_indel = large_indel || large_indel_tracker[idx];
            used_anchors.insert(idx);
            chain.push(matches[idx].clone());
            i = prev[idx];
        }

        if chain.len() < 2 {
            break;
        }

        chain.reverse();

        let l = chain.first().unwrap().pos2;
        let r = chain.last().unwrap().pos2;
        let interval = Interval{start: l.min(r), stop: l.max(r), val: true};
        if reference_ranges.iter().any(|x| {
            let intersect = x.intersect(&interval);
            intersect as f64 / (interval.stop - interval.start) as f64 > 0.25
        }) {
            continue;
        }
        chains.push((score, chain, large_indel));
        reference_ranges.push(interval);
    }
    
    return chains
}

fn find_optimal_chain(
    anchors: &Vec<Anchor>,
    match_score: i32,
    gap_cost: i32,
    band_opt: Option<usize>,
) -> Vec<ChainInfo> {
    let band;
    let matches = anchors;
    if band_opt.is_none() {
        band = 50;
    } else {
        band = band_opt.unwrap();
    }

    if anchors.is_empty() {
        return vec![];
    }

    let mut scores_and_chains_f = vec![];

    #[cfg(any(target_arch = "x86_64"))]
    {
        let vals = dp_anchors(matches, false, gap_cost, match_score, band);
        scores_and_chains_f.extend(vals);
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let vals = dp_anchors(matches, false, gap_cost, match_score, band);
        scores_and_chains_f.extend(vals);
    }

    let mut scores_and_chains_r = vec![];
    #[cfg(any(target_arch = "x86_64"))]
    {
        let vals = dp_anchors(matches, true, gap_cost, match_score, band);
        scores_and_chains_r.extend(vals);
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let vals = dp_anchors(matches, true, gap_cost, match_score, band);
        scores_and_chains_r.extend(vals);
    }

    if scores_and_chains_f.is_empty() && scores_and_chains_r.is_empty() {
        return vec![];
    }

    let max_score = scores_and_chains_f.iter().chain(scores_and_chains_r.iter()).map(|x| x.0).max().unwrap();
    let mut chains = vec![];
    let mut reference_intervals : Vec<Interval<u32, bool>> = vec![];
    let mut both_chains = scores_and_chains_f.into_iter().map(|x| (x, false)).chain(scores_and_chains_r.into_iter().map(|x| (x, true))).collect::<Vec<_>>();
    both_chains.sort_by(|a, b| b.0.0.cmp(&a.0.0));

    for ((score, chain, large_indel), reverse) in both_chains{
        if score as f64 > 0.75 * max_score as f64 {
            let l = chain.first().unwrap().pos2;
            let r = chain.last().unwrap().pos2;
            let interval = Interval{start: l.min(r), stop: l.max(r), val: true};

            if reference_intervals.iter().any(|x| {
                let intersect = x.intersect(&interval);
                intersect as f64 / (interval.stop - interval.start) as f64 > 0.25
            }) {
                continue;
            }
            chains.push(ChainInfo {
                chain,
                reverse: reverse,
                score: score,
                large_indel: large_indel,
            });
            reference_intervals.push(interval);


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
    compare_snpmers: bool,
    retain_chain: bool,
    args: &Cli,
) -> Vec<TwinOverlap> {
    let mini_chain_infos;
    if let Some(anchors) = mini_anchors {
        mini_chain_infos = find_optimal_chain(
            &anchors.anchors,
            args.c as i32,
            1,
            Some((anchors.max_mult * 10).min(50)),
        );
    } else {
        let anchors = find_exact_matches_indexes(seq1.minimizers(), seq2.minimizers());
        mini_chain_infos = find_optimal_chain(&anchors.0, args.c as i32, 1, Some((anchors.1 * 10).min(50)));
    }
    let mut twin_overlaps = vec![];
    let k = seq1.k as usize;
    for mini_chain_info in mini_chain_infos{
        let mini_chain = &mini_chain_info.chain;
        if mini_chain_info.score < MIN_CHAIN_SCORE_COMPARE {
            continue;
        }

        let mut shared_snpmer = usize::MAX;
        let mut diff_snpmer = usize::MAX;

        let mut intersect_split = 0;
        let mut intersection_snp = 0;

        if compare_snpmers{
            shared_snpmer = 0;
            diff_snpmer = 0;

            let l1 = seq1.minimizer_positions[mini_chain[0].i as usize] as usize;
            let r1 = seq1.minimizer_positions[mini_chain[mini_chain.len() - 1].i as usize] as usize;
            let l2 = seq2.minimizer_positions[mini_chain[0].j as usize] as usize;
            let r2 = seq2.minimizer_positions[mini_chain[mini_chain.len() - 1].j as usize] as usize;
            let start1 = l1.min(r1);
            let end1 = l1.max(r1) + k - 1;
            let start2 = l2.min(r2);
            let end2 = l2.max(r2) + k - 1;

            let mask = !(3 << (k - 1));
            
            let mut splitmers1 = vec![];
            let mut ind_redirect1 = vec![];

            for (i, (pos, snpmer)) in seq1.snpmers().enumerate() {
                if pos as usize >= start1 && pos as usize <= end1 {
                    ind_redirect1.push(i);
                    splitmers1.push((pos, snpmer & mask));
                }
            }

            let mut splitmers2 = vec![];
            let mut ind_redirect2 = vec![];

            for (i, (pos, snpmer)) in seq2.snpmers().enumerate() {
                if pos as usize >= start2 && pos as usize <= end2 {
                    ind_redirect2.push(i);
                    splitmers2.push((pos, snpmer & mask));
                }
            }
            
            let split_chain_opt;
            if let Some(anchors) = snpmer_anchors {
                split_chain_opt = find_optimal_chain(&anchors.anchors, 50, 1, Some((anchors.max_mult * 10).min(50))).into_iter().max_by_key(|x| x.score);
            } else {
                let anchors = find_exact_matches_indexes_references(&splitmers1, &splitmers2);
                let chains = find_optimal_chain(&anchors.0, 50, 1, Some((anchors.1 * 10).min(50)));
                split_chain_opt = chains.into_iter().max_by_key(|x| x.score);
            }

            //If mini chain goes opposite from split chain, probably split chain
            //is not reliable, so set shared and diff = 0. 
            if let Some(split_chain) = split_chain_opt.as_ref(){
                if split_chain.reverse == mini_chain_info.reverse || split_chain.chain.len() == 1 {
                    for anchor in split_chain.chain.iter() {
                        let i = anchor.i;
                        let i = ind_redirect1[i as usize];
                        let j = anchor.j;
                        let j = ind_redirect2[j as usize];

                        if seq1.snpmer_kmers[i as usize] == seq2.snpmer_kmers[j as usize] {
                            shared_snpmer += 1;
                        } else {
                            diff_snpmer += 1;
                        }
                    }
                }
            }

            //Only if log level is trace
            if log::log_enabled!(log::Level::Trace) && false{
                if diff_snpmer < 10 && shared_snpmer > 100{
                    let mut positions_read1_snpmer_diff = vec![];
                    let mut positions_read2_snpmer_diff = vec![];

                    let mut kmers_read1_diff = vec![];
                    let mut kmers_read2_diff = vec![];

                    for anchor in split_chain_opt.unwrap().chain.iter() {
                        let i = anchor.i;
                        let j = anchor.j;
                        if seq1.snpmer_kmers[i as usize] != seq2.snpmer_kmers[j as usize] {
                            positions_read1_snpmer_diff.push(seq1.snpmer_positions[i as usize]);
                            positions_read2_snpmer_diff.push(seq2.snpmer_positions[j as usize]);

                            let kmer1 = decode_kmer(seq1.snpmer_kmers[i as usize], seq1.k as u8);
                            let kmer2 = decode_kmer(seq2.snpmer_kmers[j as usize], seq2.k as u8);

                            kmers_read1_diff.push(kmer1);
                            kmers_read2_diff.push(kmer2);
                        }
                    }
                    log::trace!("{}--{:?} {}--{:?}, snp_diff:{} snp_shared:{}, kmers1:{:?}, kmers2:{:?}", &seq1.id, positions_read1_snpmer_diff, &seq2.id, positions_read2_snpmer_diff, diff_snpmer, shared_snpmer, kmers_read1_diff, kmers_read2_diff);
                }
            }

            intersect_split = splitmers1
                .iter()
                .map(|x| x.1)
                .collect::<FxHashSet<_>>()
                .intersection(&splitmers2.iter().map(|x| x.1).collect::<FxHashSet<_>>())
                .count();
            intersection_snp = seq1
                .snpmer_kmers
                .iter()
                .collect::<FxHashSet<_>>()
                .intersection(&seq2.snpmer_kmers.iter().collect::<FxHashSet<_>>())
                .count();
        }

        let l1 = seq1.minimizer_positions[mini_chain[0].i as usize];
        let r1 = seq1.minimizer_positions[mini_chain[mini_chain.len() - 1].i as usize];
        let l2 = seq2.minimizer_positions[mini_chain[0].j as usize];
        let r2 = seq2.minimizer_positions[mini_chain[mini_chain.len() - 1].j as usize];
        let start1 = l1.min(r1);
        let end1 = l1.max(r1);
        let start2 = l2.min(r2);
        let end2 = l2.max(r2);
        let shared_minimizers = mini_chain.len();
        let mut mini_chain_return = None;
        if retain_chain{
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
            snpmers_in_both: (seq1.snpmer_kmers.len(), seq2.snpmer_kmers.len()),
            chain_reverse: mini_chain_info.reverse,
            intersect: (intersect_split, intersection_snp),
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
        //Right now it's 0.5% penalty for large indels.
        let penalty = IDENTITY_THRESHOLDS.last().unwrap() - IDENTITY_THRESHOLDS.first().unwrap();
        let penalty = penalty / 2.;
        id_est -= penalty;
    }

    return id_est;
}

fn unitigs_to_tr(
    unitig_graph: &unitig::UnitigGraph,
    snpmer_set: &FxHashSet<u64>,
    solid_kmers: &HashSet<u64>,
    args: &Cli,
) -> FxHashMap<usize, TwinRead> {
    let mut tr_unitigs = FxHashMap::default();
    for (&node_hash_id, unitig) in unitig_graph.nodes.iter() {
        let id = format!("u{}", unitig.read_indices_ori[0].0);
        let u8_seq = unitig
            .base_seq()
            .iter()
            .map(|x| bits_to_ascii(x.to_bits()) as u8)
            .collect::<Vec<u8>>();
        let tr = seeding::get_twin_read_syncmer(u8_seq, None, args.kmer_size, args.c, &snpmer_set, id);
        if let Some(mut tr) = tr {
            let mut solid_mini_positions = FxHashSet::default();
            let mut solid_snpmer_positions = FxHashSet::default();
            for (pos,mini) in tr.minimizers(){
                if solid_kmers.contains(&mini){
                    solid_mini_positions.insert(pos as usize);
                }
            }
            for (pos, snpmer) in tr.snpmers(){
                if snpmer_set.contains(&snpmer){
                    solid_snpmer_positions.insert(pos as usize);
                }
            }
            tr.retain_mini_positions(solid_mini_positions);
            tr.retain_snpmer_positions(solid_snpmer_positions);
            tr_unitigs.insert(node_hash_id, tr);
        }
    }
    return tr_unitigs;
}



pub fn get_minimizer_index(twinreads: &FxHashMap<usize, TwinRead>) -> FxHashMap<u64, Vec<HitInfo>> {
    let mut mini_index = FxHashMap::default();
    for (&id, tr) in twinreads.iter() {
        for (i, (pos, mini)) in tr.minimizers().enumerate() {
            let hit = HitInfo {
                index: i as u32,
                contig_id: id as u32,
                pos: pos as u32,
            };
            mini_index.entry(mini).or_insert(vec![]).push(hit);
        }
    }
    mini_index
}

pub fn get_minimizer_index_ref(
    twinreads: &FxHashMap<usize, &TwinRead>,
) -> FxHashMap<u64, Vec<HitInfo>> {
    let mut mini_index = FxHashMap::default();
    for (&id, tr) in twinreads.iter() {
        for (i, (pos, mini)) in tr.minimizers().enumerate() {
            let hit = HitInfo {
                index: i as u32,
                contig_id: id as u32,
                pos: pos as u32,
            };
            mini_index.entry(mini).or_insert(vec![]).push(hit);
        }
    }
    mini_index
}

pub fn map_reads_to_outer_reads(
    outer_read_indices: &[usize],
    twin_reads: &[TwinRead],
    args: &Cli,
) -> Vec<TwinReadMapping> {
    let mut ret = vec![];
    // 0..outer_read_indices -- confusingly, I chose to renumber the indices in this step. Then 
    // it's fixed in the index_of_outer_in_all. 
    let tr_outer = outer_read_indices
        .iter()
        .enumerate()
        .map(|(e, &i)| (e, &twin_reads[i]))
        .collect::<FxHashMap<usize, &TwinRead>>();

    let mini_index = get_minimizer_index_ref(&tr_outer);
    let mapping_maximal_boundaries_map = Mutex::new(FxHashMap::default());
    let mapping_local_boundaries_map = Mutex::new(FxHashMap::default());
    //let kmer_count_map = Mutex::new(FxHashMap::default());
    
    let counter = Mutex::new(0);

    twin_reads.par_iter().enumerate().for_each(|(rid, read)| {
        let mini = read.minimizers_vec();
        //let start = std::time::Instant::now();
        let mini_anchors = find_exact_matches_with_full_index(&mini, &mini_index);
        //let anchor_finding_time = start.elapsed().as_micros();
        let mut unitig_hits = vec![];
        *counter.lock().unwrap() += 1;
        if *counter.lock().unwrap() % 10000 == 0 {
            log::info!("Processed {} reads / {} ...", *counter.lock().unwrap(), twin_reads.len());
        }

        //let start = std::time::Instant::now();
        for (contig_id, anchors) in mini_anchors.iter() {
            if anchors.anchors.len() < 10 {
                continue;
            }
            for twin_ol in compare_twin_reads(
                read,
                &tr_outer[&(*contig_id as usize)],
                Some(anchors),
                None,
                rid,
                *contig_id as usize,
                true,
                false,
                args
            ) {
                if twin_ol.end2 - twin_ol.start2 < 500 {
                    continue;
                }
                unitig_hits.push(twin_ol);
            }
        }
        //let overlap_time = start.elapsed().as_micros();
        for hit in unitig_hits.into_iter() {
            {
                let max_overlap = check_maximal_overlap(
                    hit.start1 as usize,
                    hit.end1 as usize,
                    hit.start2 as usize,
                    hit.end2 as usize,
                    read.base_length,
                    tr_outer[&hit.i2].base_length,
                    hit.chain_reverse,
                );

                let identity = id_est(hit.shared_minimizers, hit.diff_snpmers, args.c as u64, hit.large_indel);

                //Used for EC
                let _strain_specific_thresh_lax = same_strain(
                    hit.shared_minimizers,
                    hit.diff_snpmers,
                    hit.shared_snpmers,
                    args.c as u64,
                    args.snpmer_threshold_lax,
                    args.snpmer_error_rate_lax,
                    hit.large_indel,
                );

                //Populate mapping boundaries map
                if max_overlap {
                    let small_twin_ol = BareMappingOverlap{
                        snpmer_identity: identity as f32,
                    };
                    let mut map = mapping_maximal_boundaries_map.lock().unwrap();
                    let vec = map.entry(hit.i2).or_insert(vec![]);
                    vec.push((hit.start2 as u32 + 50, hit.end2 as u32 - 50, small_twin_ol));

                } 

                let mut map = mapping_local_boundaries_map.lock().unwrap();
                let vec = map.entry(hit.i2).or_insert(vec![]);
                vec.push((hit.start2 as u32 + 50, hit.end2 as u32 - 50));
            }
        }
        //println!("Anchor finding time: {}us", anchor_finding_time);
        //println!("Overlap time: {}us", overlap_time);
    });

    drop(mini_index);

    //let kmer_count_map = kmer_count_map.into_inner().unwrap();
    let mut num_alignments = 0;
    let mut num_maximal = 0;

    let mapping_local_boundaries_map = mapping_local_boundaries_map.into_inner().unwrap();
    let mut mapping_maximal_boundaries_map = mapping_maximal_boundaries_map.into_inner().unwrap();

    for (outer_id, boundaries) in mapping_local_boundaries_map.into_iter() {

        let index_of_outer_in_all = outer_read_indices[outer_id];
        let outer_read_length = twin_reads[index_of_outer_in_all].base_length;
        num_alignments += boundaries.len();

        let mut all_local_intervals = boundaries.into_iter().map(|x| BareInterval{start: x.0, stop: x.1}).collect::<Vec<_>>();
        all_local_intervals.sort_unstable();

        let maximal_boundaries = std::mem::take(mapping_maximal_boundaries_map.get_mut(&outer_id).unwrap_or(&mut vec![]));

        let max_intervals = maximal_boundaries
            .into_iter()
            .map(|x| Interval {
                start: x.0 as u32,
                stop: x.1 as u32,
                val: x.2,
            })
            .collect::<Vec<Interval<u32, BareMappingOverlap>>>();

        num_maximal += max_intervals.len();
        let lapper = Lapper::new(max_intervals);
        
        let map_info = MappingInfo {
            minimum_depth: -1.,
            median_depth: -1.,
            max_alignment_boundaries: None,
            max_mapping_boundaries: Some(lapper),
            present: true,
            length: outer_read_length,
        };

        let twinread_mapping = TwinReadMapping {
            tr_index: outer_read_indices[outer_id],
            mapping_info: map_info,
            all_intervals: all_local_intervals,
        };

        ret.push(twinread_mapping);
    }
    log::info!("Number of local alignments to outer reads: {}", num_alignments);
    log::info!("Number of maximal alignments to outer reads: {}", num_maximal);

    ret
}

pub fn map_reads_to_unitigs(
    unitig_graph: &mut unitig::UnitigGraph,
    kmer_info: &KmerGlobalInfo,
    twin_reads: &[TwinRead],
    args: &Cli,
) {
    let mapping_file = Path::new(args.output_dir.as_str()).join("map_to_unitigs.paf");
    let mapping_file = Mutex::new(BufWriter::new(std::fs::File::create(mapping_file).unwrap()));
    let mut snpmer_set = FxHashSet::default();
    for snpmer_i in kmer_info.snpmer_info.iter() {
        let k = snpmer_i.k as usize;
        let snpmer1 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[0] as u64) << (k - 1));
        let snpmer2 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[1] as u64) << (k - 1));
        snpmer_set.insert(snpmer1);
        snpmer_set.insert(snpmer2);
    }

    //Convert unitigs to twinreads
    let tr_unitigs = unitigs_to_tr(unitig_graph, &snpmer_set, &kmer_info.solid_kmers, args);
    let mini_index = get_minimizer_index(&tr_unitigs);
    let mapping_boundaries_map = Mutex::new(FxHashMap::default());

    twin_reads.par_iter().enumerate().for_each(|(rid, read)| {
        let mini = read.minimizers_vec();
        let mini_anchors = find_exact_matches_with_full_index(&mini, &mini_index);
        let mut unitig_hits = vec![];

        for (contig_id, anchors) in mini_anchors.iter() {
            if anchors.anchors.len() < 10{
                continue;
            }
            for twinol in compare_twin_reads(
                read,
                &tr_unitigs[&(*contig_id as usize)],
                Some(anchors),
                None,
                rid,
                *contig_id as usize,
                true,
                true,
                args
            ) {
                //Disallow large indels because they may cause windowed POA to fail
                if twinol.end2 - twinol.start2 < MIN_READ_LENGTH || twinol.large_indel{
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

        let mut retained_hits = vec![];
        ss_hits.sort_by(|a, b| id_est(b.shared_minimizers, b.diff_snpmers, args.c.try_into().unwrap(), b.large_indel).
        partial_cmp(&id_est(a.shared_minimizers, a.diff_snpmers, args.c.try_into().unwrap(), a.large_indel)).unwrap());

        let mut already_best_hit = false;
        for hit in ss_hits{
            let read_length = twin_reads[hit.i1].base_length;
            let unitig_length = tr_unitigs[&hit.i2].base_length;
            let retained = check_maximal_overlap(hit.start1, hit.end1, hit.start2, hit.end2, read_length, unitig_length, hit.chain_reverse);

            if !retained{
                continue;
            }

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
                if already_best_hit{
                    break;
                }
                else{
                    already_best_hit = true;
                }
            }

            retained_hits.push(hit);
        }

        for hit in retained_hits {
            let alignment_result = alignment::get_full_alignment(
                &twin_reads[hit.i1].dna_seq,
                &tr_unitigs[&hit.i2].dna_seq,
                &hit,
                args,
            );

            write!(
                mapping_file.lock().unwrap(),
                "r{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tshared_mini:{}\tdiff_snp:{}\tshared_snp:{}\n",
                hit.i1,
                twin_reads[hit.i1].base_length,
                alignment_result.as_ref().unwrap().q_start,
                alignment_result.as_ref().unwrap().q_end,
                if hit.chain_reverse { "-" } else { "+" },
                tr_unitigs[&hit.i2].id,
                tr_unitigs[&hit.i2].base_length,
                alignment_result.as_ref().unwrap().r_start,
                alignment_result.as_ref().unwrap().r_end,
                hit.end2 - hit.start2,
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
    });

    drop(mini_index);

    let mut number_of_alignments = 0;
    let mut cigar_string_lengths = vec![];
    for (contig_id, boundaries_and_rid) in mapping_boundaries_map.into_inner().unwrap().into_iter() {
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
        cigar_string_lengths.extend(intervals.iter().map(|x| x.val.alignment_result.as_ref().unwrap().cigar.len()));
        let mut lapper = Lapper::new(intervals);
        lapper.intervals.shrink_to_fit();

        //let (unitig_first_mini_pos, unitig_last_mini_pos) = first_last_mini_in_range(0, unitig_length, args.kmer_size, MINIMIZER_END_NTH_COV, tr_unitigs[&contig_id].minimizers.as_slice());
        //let (min_depth, median_depth) = median_and_min_depth_from_lapper(&lapper, SAMPLING_RATE_COV, unitig_first_mini_pos, unitig_last_mini_pos).unwrap();

        // println!("Unitig id kmer counts: u{}", &unitig_graph.nodes.get(&contig_id).unwrap().node_id);
        //let median_kmer_count = sliding_window_kmer_coverages(&kmer_vec, 10, args.kmer_size as usize);
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
    log::debug!("Average cigar string length: {}", cigar_string_lengths.iter().sum::<usize>() as f64 / cigar_string_lengths.len() as f64);
}

fn _get_splitmers(snpmers: &[(usize, u64)], k: u64) -> Vec<(usize, u64)> {
    let mask = !(3 << (k - 1));
    snpmers
        .iter()
        .map(|x| (x.0, x.1 as u64 & mask))
        .collect::<Vec<(usize, u64)>>()
}


