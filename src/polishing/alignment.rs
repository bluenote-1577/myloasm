use crate::cli::Cli;
use crate::constants::GAPS;
use crate::constants::MAX_BLOCK_SIZE;
use crate::constants::MIN_BLOCK_SIZE;
use crate::constants::SUB_MATRIX;
use crate::types::*;
use bio_seq::prelude::*;
use bio_seq::seq::ReverseComplement;
use block_aligner::{cigar::*, scan_block::*, scores::*};

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

//Seed-chain-extend
pub fn get_full_alignment(
    q_seq_redir: &Seq<Dna>,
    r_seq: &Seq<Dna>,
    overlap: &TwinOverlap,
    args: &Cli,
) -> Option<AlignmentResult> {
    let k = args.kmer_size as u32;
    let chain = overlap.minimizer_chain.as_ref().unwrap();
    let q_seq;
    let q_seq_rev;
    let q_start;
    let r_start;
    let qlen = q_seq_redir.len() as u32;
    let end_ind = chain.len() - 1;
    if overlap.chain_reverse {
        q_seq_rev = Some(q_seq_redir.revcomp());
        q_seq = q_seq_rev.as_ref().unwrap();
        q_start = (qlen - chain[end_ind].pos1 - k) as usize;
        r_start = chain[end_ind].pos2 as usize;
    } else {
        q_seq = q_seq_redir;
        q_seq_rev = None;
        q_start = chain[0].pos1 as usize;
        r_start = chain[0].pos2 as usize;
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
        align_seq_to_ref_slice(&r_left_slice_u8_rev, &q_left_slice_u8_rev, &GAPS, Some(30));
    let (add_length_ref_l, add_length_q_l) = get_length_from_cigar(&left_cigar);
    let left_start_q = prev_q_pos as usize - add_length_q_l;
    let left_start_r = prev_r_pos as usize - add_length_ref_l;
    let left_cigar = left_cigar.into_iter().rev().collect::<Vec<OpLen>>();
    cigar_vec = left_cigar;

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
            let diff = (q_pos - prev_q_pos) as i32 - (r_pos - prev_r_pos) as i32;

            //Non diagonal overlapping k-mer match; skip
            if diff == 0 {
                if cigar_vec.last().unwrap().op == Operation::M {
                    cigar_vec.last_mut().unwrap().len += (q_pos - prev_q_pos) as usize;
                } else {
                    cigar_vec.push(OpLen {
                        op: Operation::M,
                        len: (q_pos - prev_q_pos) as usize,
                    });
                }
                prev_q_pos = q_pos;
                prev_r_pos = r_pos;
            }
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
    
    // --- RIGHT EXTEND ----
    // qlen = 3, end = 1, k = 2, gap = 0
    let right_gap = (qlen - prev_q_pos - k).min(r_seq.len() as u32 - prev_r_pos - k) as usize;
    let pqpos = (prev_q_pos + k) as usize;
    let prpos = (prev_r_pos + k) as usize;

    let q_right_slice = &q_seq[pqpos..pqpos + right_gap];
    let r_right_slice = &r_seq[prpos..prpos + right_gap];
    let q_right_slice_u8 = dna_slice_to_u8(q_right_slice);
    let r_right_slice_u8 = dna_slice_to_u8(r_right_slice);
    let right_cigar = align_seq_to_ref_slice(&r_right_slice_u8, &q_right_slice_u8, &GAPS, Some(30));
    let (add_length_ref_r, add_length_q_r) = get_length_from_cigar(&right_cigar);
    extend_cigar(&mut cigar_vec, right_cigar);

    log::trace!(
        "{}, {}, KMERS {}, {}, START {}, {}, start: {}, END {}, {}, STR: {}",
        fmt(&cigar_vec[0..10]),
        overlap.chain_reverse,
        &q_seq[q_start..q_start + 50].to_string(),
        &r_seq[r_start..r_start + 50].to_string(),
        &q_seq[left_start_q..left_start_q + 50].to_string(),
        &r_seq[left_start_r..left_start_r + 50].to_string(),
        left_start_q,
        q_seq[pqpos..pqpos + right_gap].to_string(),
        r_seq[prpos..prpos + right_gap].to_string(),
        q_seq[0..10].to_string(),
    );

    let (cigar_length_r, cigar_length_q) = get_length_from_cigar(&cigar_vec);
    let q_end = pqpos + add_length_q_r;
    let q_start = left_start_q;
    let r_end = prpos + add_length_ref_r;
    let r_start = left_start_r;
    assert!(cigar_length_r == (r_end - r_start));
    assert!(cigar_length_q == (q_end - q_start));

    //Seed-chain-extend implementation

    return Some(AlignmentResult {
        cigar: cigar_vec,
        q_start,
        r_start,
        q_end,
        r_end,
    });
}

fn extend_cigar(cigar: &mut Vec<OpLen>, new_cigar: Vec<OpLen>) {
    for op_len in new_cigar {
        if cigar.last().unwrap().op == op_len.op {
            cigar.last_mut().unwrap().len += op_len.len;
        } else {
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
            Operation::X => 'M',
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
    if xdrop.is_some(){
        let mut a = Block::<true, true>::new(query_sliced.len(), reference_sliced.len(), MAX_BLOCK_SIZE);
        let score_mat = SUB_MATRIX;
        let reference_pad = PaddedBytes::from_bytes::<NucMatrix>(&reference_sliced, MAX_BLOCK_SIZE);
        let query_pad = PaddedBytes::from_bytes::<NucMatrix>(&query_sliced, MAX_BLOCK_SIZE);
        a.align(
            &query_pad,
            &reference_pad,
            &score_mat,
            *gaps,
            MIN_BLOCK_SIZE..=MAX_BLOCK_SIZE,
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
        let mut a = Block::<true, false>::new(query_sliced.len(), reference_sliced.len(), MAX_BLOCK_SIZE);
        let score_mat = SUB_MATRIX;
        let reference_pad = PaddedBytes::from_bytes::<NucMatrix>(&reference_sliced, MAX_BLOCK_SIZE);
        let query_pad = PaddedBytes::from_bytes::<NucMatrix>(&query_sliced, MAX_BLOCK_SIZE);
        a.align(
            &query_pad,
            &reference_pad,
            &score_mat,
            *gaps,
            MIN_BLOCK_SIZE..=MAX_BLOCK_SIZE,
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
    return cigar.to_vec();
}

