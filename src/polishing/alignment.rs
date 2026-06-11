use crate::cli::Cli;
use crate::constants::GAPS;
use crate::constants::SUB_MATRIX;
use crate::types::*;
use bio_seq::prelude::*;
use block_aligner::{cigar::*, scan_block::*, scores::*};
use fxhash::FxHashMap;
use crate::seeding::minimizer_seeds_positions;

fn get_length_from_cigar(cigar: &Vec<OpLen>) -> (usize, usize) {
    let mut add_length_ref = 0;
    let mut add_length_q = 0;
    for op_len in cigar {
        match op_len.op {
            Operation::M | Operation::X | Operation::Eq => {
                add_length_ref += op_len.len;
                add_length_q += op_len.len;
            }
            Operation::I => {
                add_length_q += op_len.len;
            }
            Operation::D => {
                add_length_ref += op_len.len;
            }
            _ => {}
        }
    }
    return (add_length_ref, add_length_q);
}

// Just extend the ends to get the full alignment range
pub fn extend_ends_chain(
    q_seq_redir: &Seq<Dna>,
    r_seq: &Seq<Dna>,
    overlap: &TwinOverlap,
    args: &Cli,
) -> (usize, usize, usize, usize){
    let k = args.kmer_size as u32;
    let chain = overlap.minimizer_chain.as_ref().unwrap();
    let q_seq;
    let _q_seq_rev;
    let qlen = q_seq_redir.len() as u32;
    let end_ind = chain.len() - 1;
    if overlap.chain_reverse {
        _q_seq_rev = Some(q_seq_redir.to_revcomp());
        q_seq = _q_seq_rev.as_ref().unwrap();
    } else {
        q_seq = q_seq_redir;
        _q_seq_rev = None;
    }

    //Main chain-extend between ends
    let prev_q_pos = if overlap.chain_reverse {
        qlen - chain[end_ind].pos1 - k
    } else {
        chain[0].pos1
    };

    let prev_r_pos = if overlap.chain_reverse {
        chain[end_ind].pos2
    } else {
        chain[0].pos2
    };


    // ----- LEFT EXTEND ---- 
    let left_ref = if prev_q_pos > prev_r_pos {
        0
    } else {
        prev_r_pos - prev_q_pos
    };

    let left_query = if prev_q_pos > prev_r_pos {
        prev_q_pos - prev_r_pos
    } else {
        0
    };

    let q_left_slice = &q_seq[left_query as usize..prev_q_pos as usize];
    let r_left_slice = &r_seq[left_ref as usize..prev_r_pos as usize];
    let q_left_slice_u8 = dna_slice_to_u8(q_left_slice);
    let r_left_slice_u8 = dna_slice_to_u8(r_left_slice);
    //Reverse and do forward local
    let q_left_slice_u8_rev = q_left_slice_u8
        .iter()
        .rev()
        .map(|&x| x)
        .collect::<Vec<u8>>();
    let r_left_slice_u8_rev = r_left_slice_u8
        .iter()
        .rev()
        .map(|&x| x)
        .collect::<Vec<u8>>();
    let left_cigar_flipped =
        align_seq_to_ref_slice(&r_left_slice_u8_rev, &q_left_slice_u8_rev, &GAPS, Some(100));
    let (add_length_ref_l, add_length_q_l) = get_length_from_cigar(&left_cigar_flipped);
    let left_start_q = prev_q_pos as usize - add_length_q_l;
    let left_start_r = prev_r_pos as usize - add_length_ref_l;

    // --- RIGHT EXTEND ----
    // qlen = 3, end = 1, k = 2, gap = 0
    let end_q_pos = if overlap.chain_reverse {
        qlen - chain[0].pos1 - k
    } else {
        chain[end_ind].pos1
    };

    let end_r_pos = if overlap.chain_reverse {
        chain[0].pos2 
    } else {
        chain[end_ind].pos2
    };

    let right_gap = (qlen - end_q_pos - k).min(r_seq.len() as u32 - end_r_pos - k) as usize;
    let pqpos = (end_q_pos + k) as usize;
    let prpos = (end_r_pos + k) as usize;

    let q_right_slice = &q_seq[pqpos..pqpos + right_gap];
    let r_right_slice = &r_seq[prpos..prpos + right_gap];
    let q_right_slice_u8 = dna_slice_to_u8(q_right_slice);
    let r_right_slice_u8 = dna_slice_to_u8(r_right_slice);
    let right_cigar = align_seq_to_ref_slice(&r_right_slice_u8, &q_right_slice_u8, &GAPS, Some(100));
    let (add_length_ref_r, add_length_q_r) = get_length_from_cigar(&right_cigar);

    let q_end = pqpos + add_length_q_r;
    let q_start = left_start_q;
    let r_end = prpos + add_length_ref_r;
    let r_start = left_start_r;
    return (q_start, q_end, r_start, r_end);
}

//Seed-chain-extend
pub fn get_full_alignment(
    q_seq_redir: &Seq<Dna>,
    r_seq: &Seq<Dna>,
    overlap: &TwinOverlap,
    args: &Cli,
) -> Option<AlignmentResult> {
    //TODO
    log::trace!("Getting alignment between r{} and r{}", overlap.i1, overlap.i2);
    let k = args.kmer_size as u32;
    let chain = overlap.minimizer_chain.as_ref().unwrap();
    let q_seq;
    let _q_seq_rev;
    let qlen = q_seq_redir.len() as u32;
    let end_ind = chain.len() - 1;
    if overlap.chain_reverse {
        _q_seq_rev = Some(q_seq_redir.to_revcomp());
        q_seq = _q_seq_rev.as_ref().unwrap();
    } else {
        q_seq = q_seq_redir;
        _q_seq_rev = None;
    }

    //Main chain-extend between ends
    let mut prev_q_pos = if overlap.chain_reverse {
        qlen - chain[end_ind].pos1 - k
    } else {
        chain[0].pos1
    };

    let mut prev_r_pos = if overlap.chain_reverse {
        chain[end_ind].pos2
    } else {
        chain[0].pos2
    };

    let mut cigar_vec;

    // ----- LEFT EXTEND ---- 
    let left_ref = if prev_q_pos > prev_r_pos {
        0
    } else {
        prev_r_pos - prev_q_pos
    };

    let q_left_slice = &q_seq[0..prev_q_pos as usize];
    let r_left_slice = &r_seq[left_ref as usize..prev_r_pos as usize];
    let q_left_slice_u8 = dna_slice_to_u8(q_left_slice);
    let r_left_slice_u8 = dna_slice_to_u8(r_left_slice);
    //Reverse and do forward local
    let q_left_slice_u8_rev = q_left_slice_u8
        .iter()
        .rev()
        .map(|&x| x)
        .collect::<Vec<u8>>();
    let r_left_slice_u8_rev = r_left_slice_u8
        .iter()
        .rev()
        .map(|&x| x)
        .collect::<Vec<u8>>();
    let left_cigar =
        align_seq_to_ref_slice(&r_left_slice_u8_rev, &q_left_slice_u8_rev, &GAPS, Some(10));
    //let left_cigar = vec![]; //TODO
    let (add_length_ref_l, add_length_q_l) = get_length_from_cigar(&left_cigar);
    let left_start_q = prev_q_pos as usize - add_length_q_l;
    let left_start_r = prev_r_pos as usize - add_length_ref_l;
    let left_cigar = left_cigar.into_iter().rev().collect::<Vec<OpLen>>();
    cigar_vec = left_cigar;

    //TODO
    //log::trace!("Left cigar: {}", fmt(&cigar_vec));

    // --- CHAIN + EXTEND ---- 
    cigar_vec.push(OpLen {
        op: Operation::M,
        len: k as usize,
    });

    for i in 1..chain.len() {
        
        let q_pos = if overlap.chain_reverse {
            qlen - chain[end_ind - i].pos1 - k
        } else {
            chain[i].pos1
        };
        assert!(q_pos >= prev_q_pos);
        let r_pos = if overlap.chain_reverse {
            chain[end_ind - i].pos2
        } else {
            chain[i].pos2
        };

        //Empty slice q
        let mut consecutive_q = false;
        let mut consecutive_r = false;
        if q_pos - prev_q_pos <= k {
            consecutive_q = true;
        }
        if r_pos - prev_r_pos <= k {
            consecutive_r = true;
        }
        if consecutive_q || consecutive_r {
            continue;
        }

        let r1 = (prev_q_pos + k) as usize..q_pos as usize;
        let r2 = (prev_r_pos + k) as usize..r_pos as usize;

        let q_slice_u8 = dna_slice_to_u8(&q_seq[r1]);
        let r_slice_u8 = dna_slice_to_u8(&r_seq[r2]);
        let cigar = align_seq_to_ref_slice(&r_slice_u8, &q_slice_u8, &GAPS, None);
        extend_cigar(&mut cigar_vec, cigar);
        extend_cigar(&mut cigar_vec, vec![OpLen {
            op: Operation::M,
            len: k as usize,
        }]);
        prev_q_pos = q_pos;
        prev_r_pos = r_pos;
    }

    //TODO
    //log::trace!("Chain cigar: {}", fmt(&cigar_vec));
    
    // --- RIGHT EXTEND ----
    // qlen = 3, end = 1, k = 2, gap = 0
    let right_gap = (qlen - prev_q_pos - k).min(r_seq.len() as u32 - prev_r_pos - k) as usize;
    let pqpos = (prev_q_pos + k) as usize;
    let prpos = (prev_r_pos + k) as usize;

    let q_right_slice = &q_seq[pqpos..pqpos + right_gap];
    let r_right_slice = &r_seq[prpos..prpos + right_gap];
    let q_right_slice_u8 = dna_slice_to_u8(q_right_slice);
    let r_right_slice_u8 = dna_slice_to_u8(r_right_slice);
    let right_cigar = align_seq_to_ref_slice(&r_right_slice_u8, &q_right_slice_u8, &GAPS, Some(10));
    //let right_cigar = vec![]; //TODO

    //  log::trace!(
    //      "Right cigar: {}",
    //      fmt(&right_cigar)
    //  );

    let (add_length_ref_r, add_length_q_r) = get_length_from_cigar(&right_cigar);
    extend_cigar(&mut cigar_vec, right_cigar);

    // --FOR Development--
    // log::trace!(
    //     "{}, {}, KMERS {}, {}, START {}, {}, start: {}, END {}, {}, STR: {}",
    //     fmt(&cigar_vec[0..10]),
    //     overlap.chain_reverse,
    //     &q_seq[q_start..q_start + 50].to_string(),
    //     &r_seq[r_start..r_start + 50].to_string(),
    //     &q_seq[left_start_q..left_start_q + 50].to_string(),
    //     &r_seq[left_start_r..left_start_r + 50].to_string(),
    //     left_start_q,
    //     q_seq[pqpos..pqpos + right_gap].to_string(),
    //     r_seq[prpos..prpos + right_gap].to_string(),
    //     q_seq[0..10].to_string(),
    // );

    let (cigar_length_r, cigar_length_q) = get_length_from_cigar(&cigar_vec);
    let q_end = pqpos + add_length_q_r;
    let q_start = left_start_q;
    let r_end = prpos + add_length_ref_r;
    let r_start = left_start_r;
    assert!(cigar_length_r == (r_end - r_start));
    assert!(cigar_length_q == (q_end - q_start));

        //Seed-chain-extend implementation

    return Some(AlignmentResult {
        cigar: OpLenVec::new(cigar_vec),
        q_start,
        r_start,
        q_end,
        r_end,
    });
}


#[inline]
fn extend_cigar(cigar: &mut Vec<OpLen>, new_cigar: Vec<OpLen>) {
    for op_len in new_cigar {
        if cigar.last().unwrap().op == op_len.op {
            cigar.last_mut().unwrap().len += op_len.len;
        } 
        else {
            cigar.push(op_len);
        }
    }
}

pub fn fmt(cigar: &[OpLen]) -> String {
    let mut s = String::new();
    for &op_len in cigar.iter() {
        let c = match op_len.op {
            Operation::M => 'M',
            Operation::Eq => 'M',
            Operation::X => 'X',
            Operation::I => 'I',
            Operation::D => 'D',
            _ => continue,
        };
        s.push_str(&format!("{}{}", op_len.len, c));
    }
    return s;
}

pub fn align_seq_to_ref_slice(
    reference_sliced: &[u8],
    query_sliced: &[u8],
    gaps: &Gaps,
    xdrop: Option<i32>,
) -> Vec<OpLen> {
    let mut cigar; 
    let bs = 16;
    if xdrop.is_some(){
        let mut a = Block::<true, true>::new(query_sliced.len(), reference_sliced.len(), bs);
        let score_mat = SUB_MATRIX;
        let reference_pad = PaddedBytes::from_bytes::<NucMatrix>(&reference_sliced, bs);
        let query_pad = PaddedBytes::from_bytes::<NucMatrix>(&query_sliced, bs);
        a.align(
            &query_pad,
            &reference_pad,
            &score_mat,
            *gaps,
            //MIN_BLOCK_SIZE..=MAX_BLOCK_SIZE,
            4 as usize..=8 as usize,
            xdrop.unwrap_or(i32::MAX),
        );
        let res = a.res();
        cigar = Cigar::new(res.query_idx, res.reference_idx);
        a.trace().cigar_eq(
            &query_pad,
            &reference_pad,
            res.query_idx,
            res.reference_idx,
            &mut cigar,
        );

    }
    else{
        let max_bs = block_size_from_seq(query_sliced.len());
        let min_bs = 4;
        let mut a = Block::<true, false>::new(query_sliced.len(), reference_sliced.len(), max_bs);
        let score_mat = SUB_MATRIX;
        let reference_pad = PaddedBytes::from_bytes::<NucMatrix>(&reference_sliced, max_bs);
        let query_pad = PaddedBytes::from_bytes::<NucMatrix>(&query_sliced, max_bs);
        a.align(
            &query_pad,
            &reference_pad,
            &score_mat,
            *gaps,
            min_bs..=max_bs,
            xdrop.unwrap_or(i32::MAX),
        );
        let res = a.res();
        cigar = Cigar::new(res.query_idx, res.reference_idx);
        a.trace().cigar_eq(
        &query_pad,
        &reference_pad,
        res.query_idx,
        res.reference_idx,
        &mut cigar,
    );

    }
    let mut cig = cigar.to_vec();
    //xdrop for blockaligner is weird because its fuzzy; remove this.
    if xdrop.is_some(){
        while cig.len() > 1 && cig.last().unwrap().op != Operation::M {
            cig.pop();
        }
    }
    cig.iter_mut().for_each(|op_len| {
        if op_len.op == Operation::Eq || op_len.op == Operation::X {
            op_len.op = Operation::M;
        }
    });

    return cig;
}

#[inline]
pub fn block_size_from_seq(seq_length: usize) -> usize {
    let mut block_size = 16;


    if seq_length > 100{
        block_size = 32;
    }
    else if seq_length > 500{
        block_size = 128;
    }
    
    return block_size;
}



pub fn align_seq_to_ref_slice_local(
    query_sliced: &[u8],
    reference_sliced: &[u8],
    gaps: &Gaps,
) -> (usize, usize, Vec<OpLen>) {
    let mut cigar; 
    let max_block = 256;
    let mut a = Block::<true, true, true>::new(query_sliced.len(), reference_sliced.len(), max_block);
    let score_mat = SUB_MATRIX;
    let reference_pad = PaddedBytes::from_bytes::<NucMatrix>(&reference_sliced, max_block);
    let query_pad = PaddedBytes::from_bytes::<NucMatrix>(&query_sliced, max_block);
    a.align(
        &query_pad,
        &reference_pad,
        &score_mat,
        *gaps,
        max_block..=max_block,
        i32::MAX,
    );
    let res = a.res();
    cigar = Cigar::new(res.query_idx, res.reference_idx);
    a.trace().cigar_eq(
        &query_pad,
        &reference_pad,
        res.query_idx,
        res.reference_idx,
        &mut cigar,
    );

    // log::trace!("cig:{}, queryend:{}, refend:{}, score:{}", fmt(&cigar.to_vec()), 
    //     res.query_idx, 
    //     res.reference_idx,
    //     res.score,
    // );

    //Unlikely, very low score indicates that we should not trust these alignments, just concatenate instead. 
    if res.score < 10 && query_sliced.len().min(reference_sliced.len()) > 100 {
        return (query_sliced.len(), 1, vec![]);
    }
    
    return (res.query_idx, res.reference_idx, cigar.to_vec());
}

const SEED_K: usize = 21;
const SEED_W: usize = 25;
const MIN_ANCHORS: usize = 10;
const MAX_HITS_PER_KMER: usize = 4;

/// Minimizer index built once over a reference sequence; reused for many query alignments.
pub struct SeedChainAligner<'a> {
    reference: &'a [u8],
    index: FxHashMap<u64, Vec<u32>>,
    use_seeds: bool,
}

impl<'a> SeedChainAligner<'a> {
    pub fn new(reference: &'a [u8]) -> Self {
        let mut ref_kmers: Vec<u64> = Vec::new();
        let mut ref_positions: Vec<u64> = Vec::new();
        minimizer_seeds_positions(reference, &mut ref_kmers, &mut ref_positions, SEED_W, SEED_K);

        let use_seeds = ref_kmers.len() >= MIN_ANCHORS;
        let mut index: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
        for (hash, pos) in ref_kmers.iter().zip(ref_positions.iter()) {
            index.entry(*hash).or_default().push(*pos as u32);
        }
        SeedChainAligner { reference, index, use_seeds }
    }

    /// Align `query` against the stored reference.
    /// Returns `(query_end, ref_end, cigar)` — same contract as `align_seq_to_ref_slice_local`.
    /// Falls back to full-DP local alignment when too few seeds are found.
    pub fn align(&self, query: &[u8], gaps: &Gaps) -> (usize, usize, Vec<OpLen>) {
        if !self.use_seeds {
            return align_seq_to_ref_slice_local(query, self.reference, gaps);
        }

        // --- Anchor lookup ---
        let mut q_kmers: Vec<u64> = Vec::new();
        let mut q_positions: Vec<u64> = Vec::new();
        minimizer_seeds_positions(query, &mut q_kmers, &mut q_positions, SEED_W, SEED_K);

        let mut anchors: Vec<(u32, u32)> = Vec::new();
        for (hash, q_pos) in q_kmers.iter().zip(q_positions.iter()) {
            if let Some(r_positions) = self.index.get(hash) {
                if r_positions.len() <= MAX_HITS_PER_KMER {
                    for &r_pos in r_positions {
                        anchors.push((*q_pos as u32, r_pos));
                    }
                }
            }
        }

        if anchors.len() < MIN_ANCHORS {
            return align_seq_to_ref_slice_local(query, self.reference, gaps);
        }

        // --- Chain ---
        let chain = chain_anchors_diagonal(&anchors, 50);
        if chain.len() < MIN_ANCHORS {
            return align_seq_to_ref_slice_local(query, self.reference, gaps);
        }

        // --- Left extend (reversed local, same as get_full_alignment) ---
        let (q0, r0) = chain[0];
        let r_left_start = if q0 > r0 { 0u32 } else { r0 - q0 };
        let q_left_rev: Vec<u8> = query[..q0 as usize].iter().rev().copied().collect();
        let r_left_rev: Vec<u8> =
            self.reference[r_left_start as usize..r0 as usize].iter().rev().copied().collect();
        let left_cigar_fwd = if q_left_rev.is_empty() || r_left_rev.is_empty() {
            vec![]
        } else {
            align_seq_to_ref_slice(&r_left_rev, &q_left_rev, gaps, Some(10))
        };
        let (left_r_len, left_q_len) = get_length_from_cigar(&left_cigar_fwd);
        let left_start_q = q0 as usize - left_q_len;
        let left_start_r = r0 as usize - left_r_len;
        let left_cigar: Vec<OpLen> = left_cigar_fwd.into_iter().rev().collect();

        // --- Chain extend ---
        let mut cigar_vec: Vec<OpLen> = left_cigar;
        let first_anchor_op = OpLen { op: Operation::M, len: SEED_K };
        if cigar_vec.is_empty() {
            cigar_vec.push(first_anchor_op);
        } else {
            extend_cigar(&mut cigar_vec, vec![first_anchor_op]);
        }

        let mut prev_q = q0;
        let mut prev_r = r0;
        for &(q_pos, r_pos) in &chain[1..] {
            // skip anchors that are too close (overlap with previous k-mer span)
            if q_pos.saturating_sub(prev_q) <= SEED_K as u32
                || r_pos.saturating_sub(prev_r) <= SEED_K as u32
            {
                continue;
            }
            let q_gap = &query[(prev_q + SEED_K as u32) as usize..q_pos as usize];
            let r_gap =
                &self.reference[(prev_r + SEED_K as u32) as usize..r_pos as usize];
            if !q_gap.is_empty() && !r_gap.is_empty() {
                let gap_cigar = align_seq_to_ref_slice(r_gap, q_gap, gaps, None);
                extend_cigar(&mut cigar_vec, gap_cigar);
            } else if !q_gap.is_empty() {
                extend_cigar(&mut cigar_vec, vec![OpLen { op: Operation::I, len: q_gap.len() }]);
            } else if !r_gap.is_empty() {
                extend_cigar(&mut cigar_vec, vec![OpLen { op: Operation::D, len: r_gap.len() }]);
            }
            extend_cigar(&mut cigar_vec, vec![OpLen { op: Operation::M, len: SEED_K }]);
            prev_q = q_pos;
            prev_r = r_pos;
        }

        // --- Right extend ---
        let after_q = (prev_q + SEED_K as u32) as usize;
        let after_r = (prev_r + SEED_K as u32) as usize;
        let right_gap = query.len().saturating_sub(after_q)
            .min(self.reference.len().saturating_sub(after_r));
        if right_gap > 0 {
            let right_cigar = align_seq_to_ref_slice(
                &self.reference[after_r..after_r + right_gap],
                &query[after_q..after_q + right_gap],
                gaps,
                Some(10),
            );
            extend_cigar(&mut cigar_vec, right_cigar);
        }

        let (total_r, total_q) = get_length_from_cigar(&cigar_vec);
        let ref_end = left_start_r + total_r;
        let query_end = left_start_q + total_q;
        (query_end, ref_end, cigar_vec)
    }
}

/// Diagonal-based anchor chaining.
/// Sorts anchors, finds the median diagonal, keeps anchors within `bw` bandwidth,
/// then enforces monotone reference positions.
fn chain_anchors_diagonal(anchors: &[(u32, u32)], bw: i64) -> Vec<(u32, u32)> {
    let mut sorted = anchors.to_vec();
    sorted.sort_unstable();
    sorted.dedup();

    let mut diags: Vec<i64> = sorted.iter()
        .map(|&(q, r)| r as i64 - q as i64)
        .collect();
    diags.sort_unstable();
    let median_diag = diags[diags.len() / 2];

    let mut chained: Vec<(u32, u32)> = Vec::new();
    let mut last_r: u32 = 0;
    let mut last_score = i64::MAX;
    let mut last_q: u32 = 0;
    for &(q_pos, r_pos) in &sorted {
        let diag = r_pos as i64 - q_pos as i64;
        let score = (diag - median_diag).abs();
        if score <= bw && r_pos >= last_r {
            // Multiple anchors for the same k-mer
            if q_pos == last_q{
                if score > last_score {
                    continue;
                }
            } else {
                last_score = score;
            }
            last_r = r_pos + SEED_K as u32;
            chained.push((q_pos, r_pos));
        }
    }
    chained
}