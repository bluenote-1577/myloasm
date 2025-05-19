use crate::cli::Cli;
use crate::constants::IDENTITY_THRESHOLDS;
use crate::constants::ID_THRESHOLD_ITERS;
use crate::constants::MAX_ALLOWABLE_SNPMER_ERROR_MISC;
use crate::constants::MAX_MULTIPLICITY_KMER;
use crate::constants::MIN_CHAIN_SCORE_COMPARE;
use crate::constants::MIN_READ_LENGTH;
use crate::constants::USE_SOLID_KMERS;
use crate::graph::GraphNode;
use crate::map_processing::*;
use crate::twin_graph;
use crate::polishing::alignment;
use crate::seeding;
use crate::twin_graph::same_strain;
use crate::types::*;
use crate::unitig;
use crate::utils::log_memory_usage;
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
    pub index: u32,
    pub contig_id: u32,
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
    seq1: &[(u32, Kmer48)],
    index: &FxHashMap<Kmer48, Vec<HitInfo>>,
    reference_seqs_owned: Option<&FxHashMap<usize, TwinRead>>,
    reference_seqs_ref: Option<&FxHashMap<usize, &TwinRead>>,
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
                let anchor = AnchorBuilder {
                    i: i as u32,
                    j: hit.index as u32,
                    pos1: *pos as u32,
                };

                matches.entry(contig).or_insert(vec![]).push(anchor);
            }
        }
    }

    matches
        .into_iter()
        .map(|(k, v)| {
            let reference_kmer_positions; 
            if let Some(reference_seqs_owned) = reference_seqs_owned {
                reference_kmer_positions = &reference_seqs_owned[&(k as usize)].minimizer_positions;
            } else if let Some(reference_seqs_ref) = reference_seqs_ref {
                reference_kmer_positions = &reference_seqs_ref[&(k as usize)].minimizer_positions;
            } else {
                panic!("No reference sequences provided");
            }
            let mut anchors = v
                .into_iter()
                .map(|anchor| {
                    let pos2 = reference_kmer_positions[anchor.j as usize];
                    Anchor {
                        i: anchor.i,
                        j: anchor.j,
                        pos1: anchor.pos1,
                        pos2: pos2 as u32,
                    }
                })
                .collect::<Vec<_>>();
            anchors.sort_by(|a, b| (a.pos1, a.pos2).cmp(&(b.pos1, b.pos2)));
            return (k,
            Anchors {
                anchors,
                max_mult,
            });
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

fn find_exact_matches_indexes(seq1: &Vec<(u32,Kmer48)> , seq2: &Vec<(u32, Kmer48)>) -> (Vec<Anchor>, usize)
where
{
    let mut max_mult = 0;
    let mut matches = Vec::new();
    let mut index_map = FxHashMap::default();

    //Sorted
    for (j, &(pos, x)) in seq2.iter().enumerate() {
        index_map.entry(x).or_insert(vec![]).push((j, pos));
    }

    //Sorted
    for (i, &(pos, s)) in seq1.iter().enumerate() {
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
    max_gap: usize,
    double_gap: usize,
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

            if (start1 as i32 - end1 as i32).abs() > (double_gap as i32)
                || (start2 as i32 - end2 as i32).abs() > (double_gap as i32)
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
                }
                if large_indel {
                    large_indel_tracker[i] = true;
                }
            } else {
                unimproved += 1;
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

    let mut scores_and_chains_f = vec![];

    #[cfg(any(target_arch = "x86_64"))]
    {
        let vals = dp_anchors(matches, false, gap_cost, match_score, band, max_gap, double_gap);
        scores_and_chains_f.extend(vals);
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let vals = dp_anchors(matches, false, gap_cost, match_score, band, max_gap, double_gap);
        scores_and_chains_f.extend(vals);
    }

    let mut scores_and_chains_r = vec![];
    #[cfg(any(target_arch = "x86_64"))]
    {
        let vals = dp_anchors(matches, true, gap_cost, match_score, band, max_gap, double_gap);
        scores_and_chains_r.extend(vals);
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let vals = dp_anchors(matches, true, gap_cost, match_score, band, max_gap, double_gap);
        scores_and_chains_r.extend(vals);
    }

    if scores_and_chains_f.is_empty() && scores_and_chains_r.is_empty() {
        return vec![];
    }

    let max_score = scores_and_chains_f
        .iter()
        .chain(scores_and_chains_r.iter())
        .map(|x| x.0)
        .max()
        .unwrap();
    let mut chains = vec![];
    let mut reference_intervals: Vec<Interval<u32, bool>> = vec![];
    let mut query_intervals: Vec<Interval<u32, bool>> = vec![];
    let mut both_chains = scores_and_chains_f
        .into_iter()
        .map(|x| (x, false))
        .chain(scores_and_chains_r.into_iter().map(|x| (x, true)))
        .collect::<Vec<_>>();
    both_chains.sort_by(|a, b| b.0 .0.cmp(&a.0 .0));

    for ((score, chain, large_indel), reverse) in both_chains {
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
                continue;
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
                    if tr_options.force_one_to_one_alignments{
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
    let mini_chain_infos;
    if let Some(anchors) = mini_anchors {
        mini_chain_infos = find_optimal_chain(
            &anchors.anchors,
            args.c as i32,
            1,
            Some((anchors.max_mult * 10).min(MAX_MULTIPLICITY_KMER)),
            options,
        );
    } else {
        let anchors;
        if let Some(seq1_minimizers) = options.read1_mininimizers.as_ref(){
            anchors = find_exact_matches_indexes(seq1_minimizers, &seq2.minimizers_vec());
        }
        else{
            anchors = find_exact_matches_indexes(&seq1.minimizers_vec(), &seq2.minimizers_vec());
        }
        mini_chain_infos = find_optimal_chain(
            &anchors.0,
            args.c as i32,
            1,
            Some((anchors.1 * 10).min(MAX_MULTIPLICITY_KMER)),
            options,
        );
    }
    let mut twin_overlaps = vec![];
    let k = seq1.k as usize;

    let temp_vec;
    let mut snpmer_vec = &vec![];

    let temp_vec2;
    let mut snpmer_vec_2 = &vec![];
    if options.compare_snpmers{
        if let Some(snpmers_vec_1) = options.read1_snpmers.as_ref() {
            snpmer_vec = snpmers_vec_1;
        } else {
            temp_vec = seq1.snpmers_vec();
            snpmer_vec = &temp_vec;
        }

        temp_vec2 = seq2.snpmers_vec();
        snpmer_vec_2 = &temp_vec2;
    }

    for mini_chain_info in mini_chain_infos {
        let mini_chain = &mini_chain_info.chain;
        if mini_chain_info.score < MIN_CHAIN_SCORE_COMPARE {
            continue;
        }

        let mut shared_snpmer = usize::MAX;
        let mut diff_snpmer = usize::MAX;

        if options.compare_snpmers {
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

            for (i, &(pos, snpmer)) in snpmer_vec.iter().enumerate() {
                if pos as usize >= start1 && pos as usize <= end1 {
                    ind_redirect1.push(i);
                    splitmers1.push((pos, snpmer.to_u64() & mask));
                }
            }

            let mut splitmers2 = vec![];
            let mut ind_redirect2 = vec![];

            for (i, &(pos, snpmer)) in snpmer_vec_2.iter().enumerate() {
                if pos as usize >= start2 && pos as usize <= end2 {
                    ind_redirect2.push(i);
                    splitmers2.push((pos, snpmer.to_u64() & mask));
                }
            }

            let split_chain_opt;
            if let Some(anchors) = snpmer_anchors {
                split_chain_opt = find_optimal_chain(
                    &anchors.anchors,
                    50,
                    1,
                    Some((anchors.max_mult * 10).min(50)),
                    options,
                )
                .into_iter()
                .max_by_key(|x| x.score);
            } else {
                let anchors = find_exact_matches_indexes_references(&splitmers1, &splitmers2);
                let chains = find_optimal_chain(
                    &anchors.0,
                    50,
                    1,
                    Some((anchors.1 * 10).min(50)),
                    options,
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
                        let i = ind_redirect1[i as usize];
                        let j = anchor.j;
                        let j = ind_redirect2[j as usize];

                        //if seq1.snpmer_kmers[i as usize] == seq2.snpmer_kmers[j as usize] {
                        if snpmer_kmers_seq1[i as usize].1 == snpmer_kmers_seq2[j as usize].1 {
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

                    for anchor in split_chain_opt.unwrap().chain.iter() {
                        let i = anchor.i;
                        let i = ind_redirect1[i as usize];
                        let j = anchor.j;
                        let j = ind_redirect2[j as usize];
                        if snpmer_kmers_seq1[i as usize] != snpmer_kmers_seq2[j as usize] {
                            positions_read1_snpmer_diff.push(seq1.snpmer_positions[i as usize]);
                            positions_read2_snpmer_diff.push(seq2.snpmer_positions[j as usize]);

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

        let l1 = seq1.minimizer_positions[mini_chain[0].i as usize];
        let r1 = seq1.minimizer_positions[mini_chain[mini_chain.len() - 1].i as usize];
        let l2 = seq2.minimizer_positions[mini_chain[0].j as usize];
        let r2 = seq2.minimizer_positions[mini_chain[mini_chain.len() - 1].j as usize];
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
            snpmers_in_both: (seq1.snpmer_positions.len(), seq2.snpmer_positions.len()),
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
    let mut tr_unitigs = FxHashMap::default();
    for (&node_hash_id, unitig) in unitig_graph.nodes.iter() {
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
            tr_unitigs.insert(node_hash_id, tr);
        }
    }
    return tr_unitigs;
}

pub fn get_minimizer_index(
    tr_owned: Option<&FxHashMap<usize, TwinRead>>, 
    tr_ref: Option<&FxHashMap<usize, &TwinRead>>) -> FxHashMap<Kmer48, Vec<HitInfo>> {
    let mut mini_index = FxHashMap::default();
    if let Some(twinreads) = tr_owned {
        let mut sorted_keys = twinreads.keys().collect::<Vec<_>>();
        sorted_keys.sort();
        for (&id, tr) in sorted_keys.iter().map(|&x| (x, &twinreads[x])) {
            for (i, (_, mini)) in tr.minimizers_vec().into_iter().enumerate() {
                let hit = HitInfo {
                    index: i as u32,
                    contig_id: id as u32,
                };
                mini_index.entry(mini).or_insert(vec![]).push(hit);
            }
        }
    }
    else if let Some(twinreads) = tr_ref {
        let mut sorted_keys = twinreads.keys().collect::<Vec<_>>();
        sorted_keys.sort();
        for (&id, tr) in sorted_keys.iter().map(|&x| (x, &twinreads[x])) {
            for (i, (_, mini)) in tr.minimizers_vec().into_iter().enumerate() {
                let hit = HitInfo {
                    index: i as u32,
                    contig_id: id as u32,
                };
                mini_index.entry(mini).or_insert(vec![]).push(hit);
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

    return mini_index;
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
    log::debug!("ITERATIONS: Breaking {} reads into {} chunks of <= 1 million", outer_read_indices.len(), outer_read_chunks.len());

    for mapping_chunk_indices in outer_read_chunks
    {
        let mapping_maximal_boundaries_map = Mutex::new(FxHashMap::default());
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
            tr_options.double_gap = 750;
            tr_options.force_one_to_one_alignments = true;
            tr_options.read1_snpmers = Some(read.snpmers_vec());
            let mini = read.minimizers_vec();
            //let start = std::time::Instant::now();
            let mini_anchors = find_exact_matches_with_full_index(&mini, &mini_index, None, Some(&tr_outer));
            drop(mini);
            //let anchor_finding_time = start.elapsed().as_micros();
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
                        args.maximal_end_fuzz,
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
                            let mut map = mapping_maximal_boundaries_map.lock().unwrap();
                            let vec = map.entry(hit.i2).or_insert(vec![]);
                            vec.push((hit.start2 as u32 + 50, hit.end2 as u32 - 50, small_twin_ol));
                        }
                    }

                    // TODO Require length conditio that scales with hit.i2 (the target read)
                    let mut map = mapping_local_boundaries_map.lock().unwrap();
                    let vec = map.entry(hit.i2).or_insert(vec![]);
                    vec.push((hit.start2 as u32 + 50, hit.end2 as u32 - 50));
                }
            }
        });

        drop(mini_index);
        
        let mapping_local_boundaries_map = mapping_local_boundaries_map.into_inner().unwrap();
        let mut mapping_maximal_boundaries_map = mapping_maximal_boundaries_map.into_inner().unwrap();

        for (outer_id, boundaries) in mapping_local_boundaries_map.into_iter() {
            //let index_of_outer_in_all = outer_read_indices[outer_id];
            let outer_read_length = twin_reads[outer_id].base_length;
            num_alignments += boundaries.len();

            let mut all_local_intervals = boundaries
                .into_iter()
                .map(|x : (u32,u32)| BareInterval {
                    start: x.0,
                    stop: x.1,
                })
                .collect::<Vec<_>>();

            all_local_intervals.sort_unstable();
            all_local_intervals.shrink_to_fit();

            let maximal_boundaries = std::mem::take(
                mapping_maximal_boundaries_map
                    .get_mut(&outer_id)
                    .unwrap_or(&mut vec![]),
            );

            let mut max_intervals = maximal_boundaries
                .into_iter()
                .map(|x : (u32, u32, BareMappingOverlap)| (BareInterval{
                    start: x.0,
                    stop: x.1,
                }, x.2))
                .collect::<Vec<_>>();

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
        let k = snpmer_i.k as usize;
        let snpmer1 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[0] as u64) << (k - 1));
        let snpmer2 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[1] as u64) << (k - 1));
        snpmer_set.insert(snpmer1);
        snpmer_set.insert(snpmer2);
    }

    //Convert unitigs to twinreads
    let tr_unitigs = unitigs_to_tr(unitig_graph, &snpmer_set, &kmer_info.solid_kmers, &kmer_info.high_freq_kmers, args);
    let mini_index = get_minimizer_index(Some(&tr_unitigs), None);
    let contained_contigs = Mutex::new(FxHashSet::default());

    log::debug!("Built minimizer index of size {}", mini_index.len());
    
    tr_unitigs.par_iter().for_each(|(q_id, q_unitig)| {
        
        let q_node_unitig = &unitig_graph.nodes[q_id];

        //Remove singletons with 0 coverage; probably errors
        if q_node_unitig.read_indices_ori.len() == 1 && 
        (q_node_unitig.min_read_depth_multi.unwrap()[0] <= 1 as f64
        || q_node_unitig.min_read_depth_multi.unwrap()[ID_THRESHOLD_ITERS-1] <= 1 as f64){
            contained_contigs.lock().unwrap().insert(*q_id);
            return;
        }

        let mut tr_options = CompareTwinReadOptions::default();
        tr_options.read1_snpmers = Some(q_unitig.snpmers_vec());
        tr_options.retain_chain = true;

        let mini = q_unitig.minimizers_vec();
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
            if r_unitig.base_length * 2 >  q_unitig.base_length {
                return false;
            }
            true
        });

        // r_untig is SMALLER than q_unitig
        //for contig_id in mini_anchor_sorted_indices.iter() {
        mini_anchor_sorted_indices.par_iter().for_each(|contig_id| {
            let r_unitig = &tr_unitigs[&(*contig_id as usize)];

            if contained_contigs.lock().unwrap().contains(&(*contig_id as usize)){
                return;
            }

            let anchors = mini_anchors.get(contig_id).unwrap();
            if anchors.anchors.len() < 10{
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

                //Allow 2 snpmer mismatches for these small contigs -- this level of variation is filtered into alternate anyways
                // OR is snpmer error
                let ss_strict = ss_strict || (uni_ol.diff_snpmers <= MAX_ALLOWABLE_SNPMER_ERROR_MISC && uni_ol.shared_snpmers > 0);

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
        let k = snpmer_i.k as usize;
        let snpmer1 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[0] as u64) << (k - 1));
        let snpmer2 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[1] as u64) << (k - 1));
        snpmer_set.insert(snpmer1);
        snpmer_set.insert(snpmer2);
    }

    //Convert unitigs to twinreads
    let tr_unitigs = unitigs_to_tr(unitig_graph, &snpmer_set, &kmer_info.solid_kmers, &kmer_info.high_freq_kmers, args);
    let mini_index = get_minimizer_index(Some(&tr_unitigs), None);
    let mapping_boundaries_map = Mutex::new(FxHashMap::default());
    let counter = Mutex::new(0);
    let num_reads = twin_reads.len();

    log::info!("Index built; starting mapping");

    twin_reads.par_iter().enumerate().for_each(|(rid, read)| {
        let mut tr_options = CompareTwinReadOptions::default();
        tr_options.retain_chain = true;
        tr_options.read1_snpmers = Some(read.snpmers_vec());

        let mini = read.minimizers_vec();
        let mini_anchors = find_exact_matches_with_full_index(&mini, &mini_index, Some(&tr_unitigs), None);
        let mut unitig_hits = vec![];

        for (contig_id, anchors) in mini_anchors.iter() {
            if anchors.anchors.len() < 10{
                continue;
            }

            // if read.id.contains("8a4ba6cd-f730-4f92-af9b-5d59da3f29e7") {
            //     dbg!(contig_id);
            //     dbg!(&anchors.anchors);
            // }
            
            for twinol in compare_twin_reads(
                read,
                &tr_unitigs[&(*contig_id as usize)],
                Some(anchors),
                None,
                rid,
                *contig_id as usize,
                &tr_options,
                args
            ) {
                log::trace!("Read {} unitig {} snpmers_shared {} snpmers_diff {} range1 {}-{} range2 {}-{}", &read.id, &tr_unitigs[&(*contig_id as usize)].id, twinol.shared_snpmers, twinol.diff_snpmers, twinol.start1, twinol.end1, twinol.start2, twinol.end2);
                
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

        // imperfect_hits.sort_by(|a, b| {
        //     let id_a = id_est(a.shared_minimizers, a.diff_snpmers, args.c.try_into().unwrap(), a.large_indel);
        //     let id_b = id_est(b.shared_minimizers, b.diff_snpmers, args.c.try_into().unwrap(), b.large_indel);
        //     ((id_b - args.snpmer_threshold_lax/100. - 0.01) * b.shared_minimizers as f64).partial_cmp(&((id_a - &args.snpmer_threshold_lax / 100. - 0.01) * a.shared_minimizers as f64)).unwrap()
        // });

        imperfect_hits.sort_by(|a, b| {
            let id_a = id_est(a.shared_minimizers, a.diff_snpmers, args.c.try_into().unwrap(), a.large_indel) - 0.98;
            let id_b = id_est(b.shared_minimizers, b.diff_snpmers, args.c.try_into().unwrap(), b.large_indel) - 0.98;
            let mini_cumulative_a = (a.shared_minimizers + a.shared_snpmers) as f64 - a.diff_snpmers as f64;
            let mini_cumulative_b = (b.shared_minimizers + b.shared_snpmers) as f64 - b.diff_snpmers as f64;
            (mini_cumulative_b * id_b).partial_cmp(&(mini_cumulative_a * id_a)).unwrap()
        });

        let mut imperfect_ids = imperfect_hits.iter().map(|x| id_est(x.shared_minimizers, x.diff_snpmers, args.c.try_into().unwrap(), x.large_indel)).collect::<Vec<_>>();
        imperfect_ids.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let mut max_perfect_mini_opt = None;
        let mut max_id = None;

        for hit in perfect_hits.iter().chain(imperfect_hits.iter()) {

            let read_length = twin_reads[hit.i1].base_length;
            let unitig_length = tr_unitigs[&hit.i2].base_length;
            let end_fuzz_pair = twin_graph::overlap_hang_length(&read, args);
            let end_fuzz = end_fuzz_pair.0.max(end_fuzz_pair.1);
            let retained = check_maximal_overlap(hit.start1, hit.end1, hit.start2, hit.end2, read_length, unitig_length, hit.chain_reverse, end_fuzz);

            log::trace!("MAPPING: {} start {} end {} ref {} snp_shared {} snp_diff {}", read.id, hit.start1, hit.end1, hit.start2, hit.shared_snpmers, hit.diff_snpmers);

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
                if let Some(max_perfect_mini) = max_perfect_mini_opt {
                    let mini_cumulative = (hit.shared_minimizers + hit.shared_snpmers) as i64 - (hit.diff_snpmers as i64);
                    if mini_cumulative < max_perfect_mini {

                        //Still take the best imperfect hit IF no perfect hits found, otherwise set a threshold at 1/3 
                        //of the distance to a perfect hit 
                        let id = id_est(hit.shared_minimizers, hit.diff_snpmers, args.c as u64, hit.large_indel);
                        let cutoff = imperfect_ids[0] - (1.0 - max_id.unwrap_or(imperfect_ids[0])) * 0.33 + 0.000000000001;
                        if id < cutoff {
                            break;
                        }
                        if (mini_cumulative as f64) < max_perfect_mini as f64 * 0.9 {
                            break;
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

        retained_hits.sort_by(|a, b| b.chain_score.partial_cmp(&a.chain_score).unwrap());
        let mut query_intervals : Vec<Interval<u32, bool>> = vec![];
        let mut max_score = 0;
        let max_chain = retained_hits.first();
        if let Some(max_chain) = max_chain{
            max_score = max_chain.chain_score;
            query_intervals.push(Interval{
                start: max_chain.start1 as u32,
                stop: max_chain.end1 as u32,
                val: true
            });
        }

        for hit in retained_hits {

            let interval = Interval{
                start: hit.start1 as u32,
                stop: hit.end1 as u32,
                val: true
            };

            //Not repopulating the query intervals for now, could change later
            if query_intervals.iter().any(|x| {
                let intersect = x.intersect(&interval);
                intersect as f64 / (interval.stop - interval.start) as f64 > 0.50
            }) {
                if (hit.chain_score as f64) < (max_score as f64) * 0.70 {
                    break;
                }
            }

            let alignment_result = alignment::get_full_alignment(
                &twin_reads[hit.i1].dna_seq,
                &tr_unitigs[&hit.i2].dna_seq,
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
                tr_unitigs[&hit.i2].id,
                tr_unitigs[&hit.i2].base_length,
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

    drop(mini_index);

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
