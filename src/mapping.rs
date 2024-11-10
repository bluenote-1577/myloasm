use crate::types::*;
use disjoint::DisjointSet;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

fn edit_distance(v1: &[u32], v2: &[u32]) -> usize {
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

pub fn disjoint_distance(
    v1: &[u32],
    v2: &[u32],
    disjoint_set: &DisjointSet,
    used_ids: &FxHashSet<u64>,
    ind1: usize,
    ind2: usize,
) -> TigdexOverlap {
    let v1_as_disjoint = v1
        .iter()
        .map(|&x| (disjoint_set.root_of(x as usize) as u64, x))
        .collect::<Vec<(u64, u32)>>();
    let v2_as_disjoint = v2
        .iter()
        .map(|&x| (disjoint_set.root_of(x as usize) as u64, x))
        .collect::<Vec<(u64, u32)>>();
    let v1_roots = v1_as_disjoint.iter().map(|x| x.0).collect::<Vec<u64>>();
    let v2_roots = v2_as_disjoint.iter().map(|x| x.0).collect::<Vec<u64>>();
    let v1_roots = v1_roots
        .into_iter()
        .enumerate()
        .collect::<Vec<(usize, u64)>>();
    let v2_roots = v2_roots
        .into_iter()
        .enumerate()
        .collect::<Vec<(usize, u64)>>();

    let chain_info = find_optimal_chain(&v1_roots, &v2_roots, 5., 1., None);
    let chain = &chain_info.chain;

    let mut shared_x = 0;
    let mut variable_roots = 0;
    let mut variable_shared_x = 0;

    for anchor in chain.iter() {
        let start1 = anchor.i;
        let start2 = anchor.j;
        if used_ids.contains(&(v1_as_disjoint[start1].1 as u64)) {
            variable_roots += 1;
            if v1_as_disjoint[start1].1 == v2_as_disjoint[start2].1 {
                variable_shared_x += 1;
            }
        }
        if v1_as_disjoint[start1].1 == v2_as_disjoint[start2].1 {
            shared_x += 1;
        }
    }

    if shared_x > v1.len() || shared_x > v2.len() {
        dbg!(shared_x, v1.len(), &chain);
        panic!();
    }

    let l1 = chain[0].pos1;
    let r1 = chain[chain.len() - 1].pos1;
    let l2 = chain[0].pos2;
    let r2 = chain[chain.len() - 1].pos2;
    let start1 = l1.min(r1);
    let end1 = l1.max(r1);
    let start2 = l2.min(r2);
    let end2 = l2.max(r2);
    let tiger = TigdexOverlap {
        tig1: ind1,
        tig2: ind2,
        tig1_start: start1,
        tig1_end: end1,
        tig2_start: start2,
        tig2_end: end2,
        shared_tig: shared_x,
        variable_roots,
        variable_tigs: variable_shared_x,
        chain_reverse: chain_info.reverse,
    };
    return tiger;
}

fn smith_waterman(
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

fn find_exact_matches_indexes(seq1: &[(usize, u64)], seq2: &[(usize, u64)]) -> (Vec<Anchor>, usize) {
    let mut max_mult = 0;
    let mut matches = Vec::new();
    let index_seq = &seq2;
    let query_seq = &seq1;
    let mut index_map = FxHashMap::default();

    //Sorted
    for (j, (pos, x)) in index_seq.iter().enumerate() {
        index_map.entry(x).or_insert(vec![]).push((j, *pos));
    }

    //Sorted
    for (i, (pos, s)) in query_seq.iter().enumerate() {
        if let Some(indices) = index_map.get(s) {
            if indices.len() > 1000 {
                continue;
            }
            if indices.len() > max_mult {
                max_mult = indices.len();
            }
            for (j, pos2) in indices {
                let anchor = Anchor {
                    i,
                    j: *j,
                    pos1: *pos,
                    pos2: *pos2,
                    score: 1.,
                };
                matches.push(anchor);
            }
        }
    }

    (matches, max_mult)
}

fn find_exact_matches_quadratic(seq1: &[usize], seq2: &[usize]) -> Vec<(usize, usize, usize)> {
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

fn dp_anchors(
    matches: &[Anchor],
    reverse: bool,
    gap_cost: f64,
    match_score: f64,
    band: usize,
) -> (f64, Vec<Anchor>) {
    let mut dp = vec![0.; matches.len()];
    let mut prev = vec![None; matches.len()];
    let mut max_score = 0.;
    let mut max_index = 0;

    for i in 0..matches.len() {
        let start1 = matches[i].pos1;
        let start2 = matches[i].pos2;
        let hit_s = matches[i].score;
        dp[i] = match_score * (hit_s as f64);
        let back = if i > band { i - band } else { 0 };
        for j in back..i {
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
            let gap_penalty =
                (start1 as i32 - (end1) as i32).abs() - (start2 as i32 - (end2) as i32).abs();
            let score = dp[j] + match_score * 1. - gap_cost * (gap_penalty.abs() as f64);
            if score > dp[i] {
                dp[i] = score;
                prev[i] = Some(j);
                if score > max_score {
                    max_score = score;
                    max_index = i;
                }
            }
        }
    }

    let mut chain = Vec::new();
    let mut i = Some(max_index);
    while let Some(idx) = i {
        chain.push(matches[idx].clone());
        i = prev[idx];
    }

    chain.reverse();
    (max_score, chain)
}

fn find_optimal_chain(
    seq1: &[(usize, u64)],
    seq2: &[(usize, u64)],
    match_score: f64,
    gap_cost: f64,
    band_opt: Option<usize>,
) -> ChainInfo {
    let (matches, max_mult) = find_exact_matches_indexes(seq1, seq2);
    let band;
    if band_opt.is_none(){
        band = (max_mult * 5).min(50);
    }
    else{
        band = band_opt.unwrap();
    }

    if matches.is_empty() {
        return ChainInfo::default();
    }

    let (score_f, chain_f) = dp_anchors(&matches, false, gap_cost, match_score, band);
    let (score_r, chain_r) = dp_anchors(&matches, true, gap_cost, match_score, band);

    if score_f > score_r {
        ChainInfo {
            chain: chain_f,
            reverse: false,
            score: score_f,
        }
    } else {
        ChainInfo {
            chain: chain_r,
            reverse: true,
            score: score_r,
        }
    }
}

pub fn compare_twin_reads(seq1: & TwinRead, seq2: & TwinRead, i: usize, j: usize) -> Option<TwinOverlap> {
    let mini_chain_info = find_optimal_chain(&seq1.minimizers, &seq2.minimizers, 10., 1., None);
    let mini_chain = &mini_chain_info.chain;
    if mini_chain_info.score < 15. {
        return None;
    }
    let k = seq1.k as u64;
    let mask = !(3 << (k - 1));
    let splitmers1 = seq1
        .snpmers
        .iter()
        .map(|x| (x.0, x.1 as u64 & mask))
        .collect::<Vec<(usize, u64)>>();
    let splitmers2 = seq2
        .snpmers
        .iter()
        .map(|x| (x.0, x.1 as u64 & mask))
        .collect::<Vec<(usize, u64)>>();
    let split_chain = find_optimal_chain(&splitmers1, &splitmers2, 50., 1., None);

    let mut shared_snpmer = 0;
    let mut diff_snpmer = 0;
    for anchor in split_chain.chain.iter() {
        let i = anchor.i;
        let j = anchor.j;
        if seq1.snpmers[i].1 == seq2.snpmers[j].1 {
            shared_snpmer += 1;
        } else {
            diff_snpmer += 1;
        }
    }

    let intersect_split = splitmers1
        .iter()
        .map(|x| x.1)
        .collect::<FxHashSet<_>>()
        .intersection(&splitmers2.iter().map(|x| x.1).collect::<FxHashSet<_>>())
        .count();
    let intersection_snp = seq1
        .snpmers
        .iter()
        .map(|x| x.1)
        .collect::<FxHashSet<_>>()
        .intersection(&seq2.snpmers.iter().map(|x| x.1).collect::<FxHashSet<_>>())
        .count();

    let l1 = seq1.minimizers[mini_chain[0].i].0;
    let r1 = seq1.minimizers[mini_chain[mini_chain.len() - 1].i].0;
    let l2 = seq2.minimizers[mini_chain[0].j].0;
    let r2 = seq2.minimizers[mini_chain[mini_chain.len() - 1].j].0;
    let start1 = l1.min(r1);
    let end1 = l1.max(r1);
    let start2 = l2.min(r2);
    let end2 = l2.max(r2);
    let twinol = TwinOverlap {
        i1: i,
        i2: j,
        start1,
        end1,
        start2,
        end2,
        shared_minimizers: mini_chain.len(),
        shared_snpmers: shared_snpmer,
        diff_snpmers: diff_snpmer,
        snpmers_in_both: (seq1.snpmers.len(), seq2.snpmers.len()),
        chain_reverse: mini_chain_info.reverse,
        intersect: (intersect_split, intersection_snp),
    };
    return Some(twinol);
}



pub fn parse_badread(id: &str) -> Option<(String, String)> {
    let spl = id.split_whitespace().collect::<Vec<&str>>()[1]
        .split(",")
        .collect::<Vec<&str>>();
    if spl.len() < 3 {
        return Some((spl[0].to_string(), "0".to_string()));
    }
    let name = spl[0];
    let range = spl[2];
    return Some((name.to_string(), range.to_string()));
}

pub fn id_est(shared_minimizers: usize, diff_snpmers: usize, c: u64) -> f64 {
    if diff_snpmers == 0 {
        return 1.;
    }
    let diff_snps = diff_snpmers as f64;
    let shared_minis = shared_minimizers as f64;
    let alpha = diff_snps as f64 / shared_minis as f64 / c as f64;
    let theta = alpha / (1. + alpha);
    let id_est = 1. - theta;
    return id_est;
}
