use crate::cli::Cli;
use crate::seeding;
use crate::twin_graph::same_strain;
use crate::types::*;
use crate::unitig;
use crate::unitig::NodeSequence;
use bio_seq::codec::Codec;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use rayon::prelude::*;
use rust_lapper::Interval;
use rust_lapper::Lapper;
use std::sync::Mutex;



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
            if indices.len() > max_mult {
                max_mult = indices.len();
            }
            for (j, pos2) in indices {
                let anchor = Anchor {
                    i: i as u32,
                    j: *j as u32,
                    pos1: *pos as u32,
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
) -> Vec<(i32, Vec<Anchor>)> {
    if matches.is_empty() {
        return vec![(0, vec![])];
    }
    let mut dp = vec![0; matches.len()];
    let mut prev = vec![None; matches.len()];
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
            let gap_penalty =
                (start1 as i32 - (end1) as i32).abs() - (start2 as i32 - (end2) as i32).abs();
            let score = dp[j] + match_score - gap_cost * (gap_penalty.abs());
            if score > dp[i] {
                dp[i] = score;
                prev[i] = Some(j);
                if score > max_score {
                    max_score = score;
                    max_index = i;
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

    let mut chain = Vec::new();
    let mut i = Some(max_index);
    while let Some(idx) = i {
        chain.push(matches[idx].clone());
        i = prev[idx];
    }

    chain.reverse();
    vec![(max_score, chain)]
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
        if is_x86_feature_detected!("avx2") {
            use crate::avx2_chaining::dp_anchors_simd;
            unsafe{
                let vals = dp_anchors_simd(matches, false, gap_cost, match_score, band);
                scores_and_chains_f.extend(vals);
            }
        }
        else{
            let vals = dp_anchors(matches, false, gap_cost, match_score, band);
            scores_and_chains_f.extend(vals);
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let vals = dp_anchors(matches, false, gap_cost, match_score, band);
        scores_and_chains_f.extend(vals);
    }

    let mut scores_and_chains_r = vec![];
    #[cfg(any(target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("avx2") {
            use crate::avx2_chaining::dp_anchors_simd;
            unsafe{
                let vals = dp_anchors_simd(matches, true, gap_cost, match_score, band);
                scores_and_chains_r.extend(vals);
            }
        }
        else{
            let vals = dp_anchors(matches, true, gap_cost, match_score, band);
            scores_and_chains_r.extend(vals);
        }
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

    for ((score, chain), reverse) in both_chains{
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
) -> Vec<TwinOverlap> {
    let mini_chain_infos;
    if let Some(anchors) = mini_anchors {
        mini_chain_infos = find_optimal_chain(
            &anchors.anchors,
            10,
            1,
            Some((anchors.max_mult * 10).min(50)),
        );
    } else {
        let anchors = find_exact_matches_indexes(&seq1.minimizers, &seq2.minimizers);
        mini_chain_infos = find_optimal_chain(&anchors.0, 10, 1, Some((anchors.1 * 10).min(50)));
    }
    let mut twin_overlaps = vec![];
    for mini_chain_info in mini_chain_infos{
        let mini_chain = &mini_chain_info.chain;
        if mini_chain_info.score < 15 {
            continue;
        }

        let l1 = seq1.minimizers[mini_chain[0].i as usize].0;
        let r1 = seq1.minimizers[mini_chain[mini_chain.len() - 1].i as usize].0;
        let l2 = seq2.minimizers[mini_chain[0].j as usize].0;
        let r2 = seq2.minimizers[mini_chain[mini_chain.len() - 1].j as usize].0;
        let start1 = l1.min(r1);
        let end1 = l1.max(r1);
        let start2 = l2.min(r2);
        let end2 = l2.max(r2);

        let k = seq1.k as u64;
        let mask = !(3 << (k - 1));
        
        let mut splitmers1 = vec![];
        let mut ind_redirect1 = vec![];

        for (i, (pos, snpmer)) in seq1.snpmers.iter().enumerate() {
            if *pos >= start1 && *pos <= end1 {
                ind_redirect1.push(i);
                splitmers1.push((*pos, *snpmer & mask));
            }
        }

        let mut splitmers2 = vec![];
        let mut ind_redirect2 = vec![];

        for (i, (pos, snpmer)) in seq2.snpmers.iter().enumerate() {
            if *pos >= start2 && *pos <= end2 {
                ind_redirect2.push(i);
                splitmers2.push((*pos, *snpmer & mask));
            }
        }
        
        let split_chain_opt;
        if let Some(anchors) = snpmer_anchors {
            split_chain_opt = find_optimal_chain(&anchors.anchors, 50, 1, Some((anchors.max_mult * 10).min(50))).into_iter().max_by_key(|x| x.score);
        } else {
            let anchors = find_exact_matches_indexes(&splitmers1, &splitmers2);
            let chains = find_optimal_chain(&anchors.0, 50, 1, Some((anchors.1 * 10).min(50)));
            split_chain_opt = chains.into_iter().max_by_key(|x| x.score);
        }

        let mut shared_snpmer = 0;
        let mut diff_snpmer = 0;
        //If mini chain goes opposite from split chain, probably split chain
        //is not reliable, so set shared and diff = 0. 
        if let Some(split_chain) = split_chain_opt.as_ref(){
            if split_chain.reverse == mini_chain_info.reverse || split_chain.chain.len() == 1 {
                for anchor in split_chain.chain.iter() {
                    let i = anchor.i;
                    let i = ind_redirect1[i as usize];
                    let j = anchor.j;
                    let j = ind_redirect2[j as usize];
                    if seq1.snpmers[i as usize].1 == seq2.snpmers[j as usize].1 {
                        shared_snpmer += 1;
                    } else {
                        diff_snpmer += 1;
                    }
                }
            }
        }

        //Only if log level is trace
        if log::log_enabled!(log::Level::Trace) {
            if diff_snpmer < 10 && shared_snpmer > 100{
                let mut positions_read1_snpmer_diff = vec![];
                let mut positions_read2_snpmer_diff = vec![];

                let mut kmers_read1_diff = vec![];
                let mut kmers_read2_diff = vec![];

                for anchor in split_chain_opt.unwrap().chain.iter() {
                    let i = anchor.i;
                    let j = anchor.j;
                    if seq1.snpmers[i as usize].1 != seq2.snpmers[j as usize].1 {
                        positions_read1_snpmer_diff.push(seq1.snpmers[i as usize].0);
                        positions_read2_snpmer_diff.push(seq2.snpmers[j as usize].0);

                        let kmer1 = decode_kmer(seq1.snpmers[i as usize].1, seq1.k as u8);
                        let kmer2 = decode_kmer(seq2.snpmers[j as usize].1, seq2.k as u8);

                        kmers_read1_diff.push(kmer1);
                        kmers_read2_diff.push(kmer2);
                    }
                }
                log::trace!("{}--{:?} {}--{:?}, snp_diff:{} snp_shared:{}, kmers1:{:?}, kmers2:{:?}", &seq1.id, positions_read1_snpmer_diff, &seq2.id, positions_read2_snpmer_diff, diff_snpmer, shared_snpmer, kmers_read1_diff, kmers_read2_diff);
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

        let l1 = seq1.minimizers[mini_chain[0].i as usize].0;
        let r1 = seq1.minimizers[mini_chain[mini_chain.len() - 1].i as usize].0;
        let l2 = seq2.minimizers[mini_chain[0].j as usize].0;
        let r2 = seq2.minimizers[mini_chain[mini_chain.len() - 1].j as usize].0;
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
        twin_overlaps.push(twinol);
    }
    return twin_overlaps;
}

pub fn parse_badread(id: &str) -> Option<(String, String)> {
    return Some(("".to_string(), "".to_string()));
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

fn unitigs_to_tr(
    unitig_graph: &unitig::UnitigGraph,
    snpmer_set: &FxHashSet<u64>,
    solid_kmers: &FxHashSet<u64>,
    args: &Cli,
) -> FxHashMap<usize, TwinRead> {
    let mut tr_unitigs = FxHashMap::default();
    for (&node_id, unitig) in unitig_graph.nodes.iter() {
        let id = format!("u{}", unitig.read_indices_ori[0].0);
        let u8_seq = unitig
            .base_seq()
            .iter()
            .map(|x| bits_to_ascii(x.to_bits()) as u8)
            .collect::<Vec<u8>>();
        let tr = seeding::get_twin_read(u8_seq, None, args.kmer_size, args.c, &snpmer_set, id);
        if let Some(mut tr) = tr {
            let mut solid_minis = vec![];
            for (pos,mini) in tr.minimizers.iter_mut(){
                if solid_kmers.contains(mini){
                    solid_minis.push((*pos,*mini));
                }
            }
            tr.minimizers = solid_minis;
            tr_unitigs.insert(node_id, tr);
        }
    }
    return tr_unitigs;
}

pub struct TwinReadMapping {
    pub tr_index: usize,
    pub mapping_info: MappingInfo,
}

impl NodeMapping for TwinReadMapping {
    fn median_depth(&self) -> f64 {
        self.mapping_info.median_depth
    }
    fn mapping_boundaries(&self) -> &Lapper<u32, bool> {
        &self.mapping_info.mapping_boundaries
    }
    fn mean_depth(&self) -> f64 {
        self.mapping_info.mean_depth
    }
    fn set_mapping_info(&mut self, mapping_info: MappingInfo) {
        self.mapping_info = mapping_info;
    }
    fn mapping_info_present(&self) -> bool {
        self.mapping_info.present
    }
    fn reference_length(&self) -> usize {
        self.mapping_info.length
    }
    fn mapped_indices(&self) -> &Vec<usize> {
        &self.mapping_info.mapped_indices
    }
}

fn get_minimizer_index(twinreads: &FxHashMap<usize, TwinRead>) -> FxHashMap<u64, Vec<HitInfo>> {
    let mut mini_index = FxHashMap::default();
    for (&id, tr) in twinreads.iter() {
        for (i, (pos, mini)) in tr.minimizers.iter().enumerate() {
            let hit = HitInfo {
                index: i,
                contig_id: id,
                pos: *pos,
            };
            mini_index.entry(*mini).or_insert(vec![]).push(hit);
        }
    }
    mini_index
}

fn get_minimizer_index_ref(
    twinreads: &FxHashMap<usize, &TwinRead>,
) -> FxHashMap<u64, Vec<HitInfo>> {
    let mut mini_index = FxHashMap::default();
    for (&id, tr) in twinreads.iter() {
        for (i, (pos, mini)) in tr.minimizers.iter().enumerate() {
            let hit = HitInfo {
                index: i,
                contig_id: id,
                pos: *pos,
            };
            mini_index.entry(*mini).or_insert(vec![]).push(hit);
        }
    }
    mini_index
}

pub fn map_reads_to_outer_reads<'a>(
    outer_read_indices: &[usize],
    twin_reads: &'a [TwinRead],
    _args: &Cli,
) -> Vec<TwinReadMapping> {
    let mut ret = vec![];
    let tr_outer = outer_read_indices
        .iter()
        .enumerate()
        .map(|(e, &i)| (e, &twin_reads[i]))
        .collect::<FxHashMap<usize, &TwinRead>>();

    let mini_index = get_minimizer_index_ref(&tr_outer);
    let mapping_boundaries_map = Mutex::new(FxHashMap::default());

    twin_reads.par_iter().enumerate().for_each(|(rid, read)| {
        let mini = &read.minimizers;
        let mini_anchors = find_exact_matches_with_full_index(mini, &mini_index);
        let mut unitig_hits = vec![];

        for (contig_id, anchors) in mini_anchors.iter() {
            if anchors.anchors.len() < 10 {
                continue;
            }
            for twin_ol in compare_twin_reads(
                read,
                &tr_outer[contig_id],
                Some(anchors),
                None,
                rid,
                *contig_id,
            ) {
                if twin_ol.end2 - twin_ol.start2 < 500 {
                    continue;
                }
                unitig_hits.push(twin_ol);
            }
        }
        for hit in unitig_hits.iter() {
            let mut map = mapping_boundaries_map.lock().unwrap();
            let vec = map.entry(hit.i2).or_insert(vec![]);
            //TODO test
            vec.push((hit.start2 + 50, hit.end2 - 50));
        }
    });

    for (contig_id, boundaries) in mapping_boundaries_map.into_inner().unwrap().into_iter() {
        let index_of_outer_in_all = outer_read_indices[contig_id];
        let outer_read_length = twin_reads[index_of_outer_in_all].base_length;
        let mut map_vec = vec![];
        for mapping in boundaries.iter() {
            let map_len = mapping.1 - mapping.0;
            map_vec.push(map_len);
        }

        let intervals = boundaries
            .into_iter()
            .map(|x| Interval {
                start: x.0 as u32,
                stop: x.1 as u32,
                val: true,
            })
            .collect::<Vec<Interval<u32, bool>>>();
        let lapper = Lapper::new(intervals);

        let mean_depth = map_vec.iter().sum::<usize>() as f64 / outer_read_length as f64;
        let map_info = MappingInfo {
            mapping_boundaries: lapper,
            mean_depth: mean_depth,
            median_depth: mean_depth,
            present: true,
            length: outer_read_length,
            mapped_indices: vec![]
        };
        let twinread_mapping = TwinReadMapping {
            tr_index: outer_read_indices[contig_id],
            mapping_info: map_info,
        };
        ret.push(twinread_mapping);
    }

    ret
}

pub fn map_reads_to_unitigs(
    unitig_graph: &mut unitig::UnitigGraph,
    kmer_info: &KmerGlobalInfo,
    twin_reads: &[TwinRead],
    args: &Cli,
) {
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
        let mini = &read.minimizers;
        let mini_anchors = find_exact_matches_with_full_index(mini, &mini_index);
        let mut unitig_hits = vec![];

        for (contig_id, anchors) in mini_anchors.iter() {
            if anchors.anchors.len() < 10{
                continue;
            }
            for twinol in compare_twin_reads(
                read,
                &tr_unitigs[contig_id],
                Some(anchors),
                None,
                rid,
                *contig_id,
            ) {
                if twinol.end2 - twinol.start2 < 500 {
                    continue;
                }
                unitig_hits.push(twinol);
            }
        }
        let ss_hits = unitig_hits.into_iter().filter(|x| same_strain(
            x.shared_minimizers,
            x.diff_snpmers,
            x.shared_snpmers,
            args.c.try_into().unwrap(),
            args.snpmer_threshold,
            args.snpmer_error_rate
        )).collect::<Vec<_>>();

        let best_hit = ss_hits.into_iter().max_by_key(|x| x.shared_minimizers + x.shared_snpmers);
        if let Some(hit) = best_hit{
            let mut map = mapping_boundaries_map.lock().unwrap();
            let vec = map.entry(hit.i2).or_insert(vec![]);
            vec.push((hit.start2, hit.end2, rid));
        }
    });


    for (contig_id, boundaries_and_rid) in mapping_boundaries_map.into_inner().unwrap().into_iter() {
        let unitig_length = unitig_graph.nodes.get(&contig_id).unwrap().length();
        let mut map_vec = vec![];
        for mapping in boundaries_and_rid.iter() {
            let map_len = mapping.1 - mapping.0;
            map_vec.push(map_len);
        }

        let intervals = boundaries_and_rid
            .iter()
            .map(|x| Interval {
                start: x.0 as u32,
                stop: x.1 as u32,
                val: true,
            })
            .collect::<Vec<Interval<u32, bool>>>();
        let lapper = Lapper::new(intervals);
        let mapping_indices = boundaries_and_rid.iter().map(|x| x.2).collect::<Vec<usize>>();

        let mean_depth = map_vec.iter().sum::<usize>() as f64 / unitig_length as f64;
        let map_info = MappingInfo {
            mapping_boundaries: lapper,
            mean_depth: mean_depth,
            median_depth: mean_depth,
            present: true,
            length: unitig_length,
            mapped_indices: mapping_indices,
        };
        let mut_node = unitig_graph.nodes.get_mut(&contig_id).unwrap();
        mut_node.set_mapping_info(map_info);
        mut_node.approx_depth = Some(mean_depth);
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct HitInfo {
    pub index: usize,
    pub contig_id: usize,
    pub pos: usize,
}

pub struct Anchors {
    anchors: Vec<Anchor>,
    max_mult: usize,
}

fn _get_splitmers(snpmers: &[(usize, u64)], k: u64) -> Vec<(usize, u64)> {
    let mask = !(3 << (k - 1));
    snpmers
        .iter()
        .map(|x| (x.0, x.1 as u64 & mask))
        .collect::<Vec<(usize, u64)>>()
}

pub fn cov_mapping_breakpoints<T>(mapped: &T) -> Vec<Breakpoints>
where
    T: NodeMapping,
{
    if mapped.mapping_boundaries().intervals.is_empty() {
        return vec![];
    }
    let mut breakpoints = vec![];
    let depths = mapped.mapping_boundaries().depth().collect::<Vec<_>>();
    if depths.len() < 3 {
        return vec![];
    }
    // <ooooo|-------
    if depths[0].start > 100 {
        breakpoints.push(Breakpoints {
            pos1: 0,
            pos2 : depths[0].start as usize,
            cov: 0,
        });
    }
    //<xxxxx|--------
    //Cut bad left endpoints. 
    let depth_start_right = if mapped.reference_length() > 100 + depths[0].stop as usize{
        mapped
            .mapping_boundaries()
            .count(depths[0].stop + 99, depths[0].stop + 100)
    } else {
        0
    };
    if depths[0].stop > 100 && (depths[1].val > 3 || depth_start_right > 3) && depths[0].val == 1 {
        breakpoints.push(Breakpoints {
            pos1: depths[0].start as usize,
            pos2: depths[0].stop as usize,
            cov: depths[0].val as usize,
        });
    }
    // -----|xxxx|----
    for i in 1..depths.len() - 1 {
        let interval = &depths[i];
        let last_cov = depths[i - 1].val;
        let next_cov = depths[i + 1].val;
        let start = interval.start;
        let stop = interval.stop as usize;
        let cov = interval.val as usize;
        let cond1 = last_cov > 3 && next_cov > 3 && cov == 1;
        let cond2;
        if start > 200 && stop + 200 < mapped.reference_length() {
            let left_count = mapped.mapping_boundaries().count(start - 200, start - 198);
            let right_count = mapped.mapping_boundaries().count(stop as u32 + 198, stop as u32 + 200);
            cond2 = left_count > 3 && right_count > 3 && cov == 1;
        } else {
            cond2 = false;
        }
        let cond3 = (last_cov > (cov as u32 * 5) || next_cov > (cov as u32 * 5)) && cov < 3;
        if start > 200
            && stop + 200 < mapped.reference_length()
            && (cond1 || cond2 || cond3)
        {
            breakpoints.push(Breakpoints {
                pos1: start as usize,
                pos2: stop as usize,
                cov: cov,
            });
        }
    }

    // --------|xxxxx>
    if depths[depths.len() - 1].start > 100 {
        let depth_stop_left = mapped.mapping_boundaries().count(
            depths[depths.len() - 1].start - 100,
            depths[depths.len() - 1].start - 99,
        );
        if (depth_stop_left > 3 || depths[depths.len() - 2].val > 3)
            && depths[depths.len() - 1].val == 1
        {
            breakpoints.push(Breakpoints {
                pos1: depths[depths.len() - 1].start as usize,
                pos2: depths[depths.len() - 1].stop as usize,
                cov: 0,
            });
        }
    }

    // -----|ooooo>
    if depths[depths.len() - 1].stop as usize + 100 < mapped.reference_length() {
        breakpoints.push(Breakpoints {
            pos1: depths[depths.len() - 1].stop as usize,
            pos2: mapped.reference_length(),
            cov: depths[depths.len() - 1].val as usize,
        });
    }

    return breakpoints;
}
