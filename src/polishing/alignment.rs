use crate::constants::SUB_MATRIX;
use crate::constants::MIN_BLOCK_SIZE;
use crate::constants::MAX_BLOCK_SIZE;
use crate::polishing::consensus::*;
use block_aligner::{cigar::*, scan_block::*, scores::*};
use std::time::Instant;

pub fn align_seq_to_ref_slice(
    reference_sliced: &[u8],
    query_sliced: &[u8],
    gaps: &Gaps,
) -> Cigar{
    let start = Instant::now();
    let mut a = Block::<true, true>::new(query_sliced.len(), reference_sliced.len(), MAX_BLOCK_SIZE);
    let score_mat = SUB_MATRIX;
    let reference_pad =
        PaddedBytes::from_bytes::<NucMatrix>(&reference_sliced, MAX_BLOCK_SIZE);
    let query_pad = PaddedBytes::from_bytes::<NucMatrix>(&query_sliced, MAX_BLOCK_SIZE);
    log::info!("Padded bytes took {:?}", start.elapsed());
    let start = Instant::now();
    a.align(
        &query_pad,
        &reference_pad,
        &score_mat,
        *gaps,
        MIN_BLOCK_SIZE..=MAX_BLOCK_SIZE,
        200,
    );
    log::info!("Alignment took {:?}", start.elapsed());
    let start = Instant::now();
    let res = a.res();
    let mut cigar = Cigar::new(res.query_idx, res.reference_idx);
    a.trace().cigar_eq(
        &query_pad,
        &reference_pad,
        res.query_idx,
        res.reference_idx,
        &mut cigar,
    );
    log::info!("Cigar took {:?}", start.elapsed());
    return cigar;
}

pub fn align_seq_to_ref_end_to_end(
    reference: &PaddedBytes,
    query: &PaddedBytes,
    q: &HomopolymerCompressedSeq,
    gaps: &Gaps,
    consensus_builder: &mut ConsensusBuilder,
) {
    let mut a = Block::<true, true>::new(query.len(), reference.len(), MAX_BLOCK_SIZE);
    let score_mat = SUB_MATRIX;
    a.align(
        &query,
        &reference,
        &score_mat,
        *gaps,
        MIN_BLOCK_SIZE..=MAX_BLOCK_SIZE,
        200,
    );
    let res = a.res();
    let mut cigar = Cigar::new(res.query_idx, res.reference_idx);
    a.trace().cigar_eq(
        &query,
        &reference,
        res.query_idx,
        res.reference_idx,
        &mut cigar,
    );
    consensus_builder.process_alignment(cigar, q, 0, 0);
}
