use crate::cli::Cli;
use crate::seeding;
use crate::twin_graph::same_strain;
use crate::types::*;
use crate::unitig;
use bio_seq::codec::Codec;
use fxhash::FxHashMap;
use fxhash::FxHashSet;

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

fn find_exact_matches_with_full_index(
    seq1: &[(usize, u64)],
    index: &FxHashMap<u64, Vec<HitInfo>>,
) -> FxHashMap<usize, Anchors> {
    let mut max_mult = 0;
    let mut matches = FxHashMap::default();

    for (i, (pos, s)) in seq1.iter().enumerate() {
        if let Some(indices) = index.get(s) {
            if indices.len() > 1000 {
                continue;
            }
            if indices.len() > max_mult {
                max_mult = indices.len();
            }
            for hit in indices {
                let contig = hit.contig_id;
                let anchor = Anchor {
                    i,
                    j: hit.index,
                    pos1: *pos,
                    pos2: hit.pos as usize,
                    score: 1.,
                };

                matches.entry(contig).or_insert(vec![]).push(anchor);
            }
        }
    }

    matches.into_iter().map(|(k, v)| (k, Anchors{anchors: v, max_mult})).collect()

}

fn find_exact_matches_indexes(
    seq1: &[(usize, u64)],
    seq2: &[(usize, u64)],
) -> (Vec<Anchor>, usize) {
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

fn chain_bothsides(matches: &[Anchor], gap_cost: f64, match_score: f64, band: usize) -> ChainInfo {
    let (score_f, chain_f) = dp_anchors(matches, false, gap_cost, match_score, band);
    let (score_r, chain_r) = dp_anchors(matches, true, gap_cost, match_score, band);

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

fn find_optimal_chain(
    anchors: &Vec<Anchor>,
    match_score: f64,
    gap_cost: f64,
    band_opt: Option<usize>,
) -> ChainInfo {
    let band;
    if band_opt.is_none() {
        band = 50;
    } else {
        band = band_opt.unwrap();
    }

    if anchors.is_empty() {
        return ChainInfo::default();
    }

    return chain_bothsides(&anchors, gap_cost, match_score, band);
}

pub fn compare_twin_reads(
    seq1: &TwinRead,
    seq2: &TwinRead,
    mini_anchors: Option<&Anchors>,
    snpmer_anchors: Option<&Anchors>,
    i: usize,
    j: usize,
) -> Option<TwinOverlap> {
    let mini_chain_info;
    if let Some(anchors) = mini_anchors {
        mini_chain_info = find_optimal_chain(&anchors.anchors, 10., 1., Some((anchors.max_mult*5).min(50)));
    } else {
        let anchors = find_exact_matches_indexes(&seq1.minimizers, &seq2.minimizers);
        mini_chain_info = find_optimal_chain(&anchors.0, 10., 1., Some((anchors.1*5).min(50)));
    }
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

    let split_chain;
    if let Some(anchors) = snpmer_anchors{
        split_chain = find_optimal_chain(&anchors.anchors, 50., 1., None);
    } else {
        let anchors = find_exact_matches_indexes(&splitmers1, &splitmers2);
        split_chain = find_optimal_chain(&anchors.0, 50., 1., None);
    }

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

fn unitigs_to_tr(unitig_graph: &unitig::UnitigGraph, args: &Cli) -> FxHashMap<usize, TwinRead> {
    let mut tr_unitigs = FxHashMap::default();
    for (&node_id, unitig) in unitig_graph.nodes.iter() {
        let id = format!("u{}", unitig.read_indices_ori[0].0);
        let u8_seq = unitig
            .base_seq
            .as_ref()
            .unwrap()
            .iter()
            .map(|x| bits_to_ascii(x.to_bits()) as u8)
            .collect::<Vec<u8>>();
        let tr = seeding::get_twin_read(u8_seq, args.kmer_size, args.c, &FxHashSet::default(), id);
        if let Some(tr) = tr {
            tr_unitigs.insert(node_id, tr);
        }
    }
    return tr_unitigs;
}

pub fn map_reads_to_unitigs(
    unitig_graph: &mut unitig::UnitigGraph,
    snpmer_info: &[SnpmerInfo],
    twin_reads: &[TwinRead],
    args: &Cli,
) {

    for node in unitig_graph.nodes.values_mut() {
        node.mapping_boundaries = Some(vec![]);
    }

    let mut snpmer_set = FxHashSet::default();
    for snpmer_i in snpmer_info.iter() {
        let k = snpmer_i.k as usize;
        let snpmer1 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[0] as u64) << (k - 1));
        let snpmer2 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[1] as u64) << (k - 1));
        snpmer_set.insert(snpmer1);
        snpmer_set.insert(snpmer2);
    }

    //Convert unitigs to twinreads
    let tr_unitigs = unitigs_to_tr(unitig_graph, args);

    let mut mini_index = FxHashMap::default();
    for (id, tr) in tr_unitigs.iter() {
        for (i, (pos, mini)) in tr.minimizers.iter().enumerate() {
            let hit = HitInfo {
                index: i,
                contig_id: *id,
                pos: *pos,
            };
            mini_index.entry(*mini).or_insert(vec![]).push(hit);
        }
    }
    let mut splitmers_unitigs = vec![];
    for (id, tr) in tr_unitigs.iter() {
        let splitmers = get_splitmers(&tr.snpmers, args.kmer_size as u64);
        splitmers_unitigs.push((*id, splitmers));
    }
    let mut splitmer_index = FxHashMap::default();
    for (id, splitmers) in splitmers_unitigs.iter() {
        for (j, (pos, snpmer)) in splitmers.iter().enumerate() {
            let hit = HitInfo {
                index: j,
                contig_id: *id,
                pos: *pos,
            };
            splitmer_index.entry(*snpmer).or_insert(vec![]).push(hit);
        }
    }

    for (rid, read) in twin_reads.iter().enumerate() {
        let mini = &read.minimizers;
        let mini_anchors = find_exact_matches_with_full_index(mini, &mini_index);
        let split_anchors =
            find_exact_matches_with_full_index(&read.snpmers, &splitmer_index);
        let mut unitig_hits = vec![];

        for (contig_id, anchors) in mini_anchors.iter() {
            if let Some(twinol) = compare_twin_reads(
                read,
                &tr_unitigs[contig_id],
                Some(anchors),
                split_anchors.get(contig_id),
                rid,
                *contig_id,
            ){

                if twinol.end2 - twinol.start2 < 500 {
                    continue;
                }
                unitig_hits.push(twinol);
            }

        }
        for hit in unitig_hits.iter() {
            if same_strain(
                hit.shared_minimizers,
                hit.diff_snpmers,
                hit.shared_snpmers,
                args.c.try_into().unwrap(),
            ) {
                unitig_graph
                    .nodes
                    .get_mut(&hit.i2)
                    .unwrap()
                    .mapping_boundaries
                    .as_mut()
                    .unwrap()
                    .push((hit.start2, hit.end2));
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct HitInfo {
    pub index: usize,
    pub contig_id: usize,
    pub pos: usize,
}

pub struct Anchors{
    anchors: Vec<Anchor>,
    max_mult: usize,
}

fn get_splitmers(snpmers: &[(usize, u64)], k: u64) -> Vec<(usize, u64)> {
    let mask = !(3 << (k - 1));
    snpmers
        .iter()
        .map(|x| (x.0, x.1 as u64 & mask))
        .collect::<Vec<(usize, u64)>>()
}
