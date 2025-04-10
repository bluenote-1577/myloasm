use crate::constants::*;
use crate::mapping;
use crate::polishing::alignment;
use crate::types::SmallTwinOl;
use crate::types::*;
use crate::cli::Cli;
use crate::seeding;
use bio_seq::prelude::*;
use fxhash::FxHashSet;
use block_aligner::cigar::*;
use fxhash::FxHashMap;
use rayon::prelude::*;
use rust_lapper::Interval;
use rust_spoa::poa_consensus;
use std::sync::Mutex;

#[derive(Debug, Default, Clone)]
pub struct BaseConsensusSimple {
    pub match_count: u32,
    pub prev_ins_weight: u32,
    pub prev_nonins_weight: u32,
}

pub struct PoaConsensusBuilder {
    seq: Vec<Mutex<Vec<Vec<u8>>>>,
    qual: Vec<Mutex<Vec<Vec<u8>>>>,
    breakpoints: Vec<usize>,
    genome_length: usize,
    window_overlap_len: usize,
    contig_name: String,
    bp_len: usize,
}

impl PoaConsensusBuilder {
    pub fn spoa_blocks(self) -> Vec<Vec<u8>> {
        let consensus_max_length = 800;
        let alignment_type = 1; // 0 local, 1 global, 2, semiglobal
        let match_score = 5;
        let mismatch_score = -2;
        let gap_open = -2;
        let gap_extend = -1;

        let consensuses = Mutex::new(vec![]);
        //parallel iter over seq and qual
        self.seq
            .into_par_iter()
            .enumerate()
            .zip(self.qual.into_par_iter())
            .for_each(|((i, seq), qual)| {
                let mut seqs = seq.lock().unwrap();
                let mut quals = qual.lock().unwrap();

                let max_len = seqs.iter().map(|x| x.len()).max().unwrap_or(0).min(self.window_overlap_len + self.bp_len);

                if seqs.len() <= 3 {
                    seqs.retain(|y| y.len() > max_len * 4 / 5);
                    quals.retain(|y| y.len() > max_len * 4 / 5);
                }

                //log::trace!("Processing block {} of {}: {:?}",i, self.contig_name, seqs.iter().map(|x| x.len()).collect::<Vec<_>>());

                let mut qual_map = FxHashMap::default();
                for (i, qual) in quals.iter().enumerate() {
                    let mean_qual = qual.iter().map(|x| *x as i32).sum::<i32>() / qual.len() as i32;
                    qual_map.insert(i, -(qual.len() as i32) * mean_qual);
                }

                //Sort seqs and quals by the qual_map
                let mut sorted_indices = qual_map.into_iter().collect::<Vec<_>>();
                sorted_indices.sort_by_key(|x| x.1);
                let sorted_indices = sorted_indices.into_iter().map(|x| x.0).collect::<Vec<_>>();

                let mut seqs = sorted_indices
                    .iter()
                    .map(|&i| std::mem::take(&mut seqs[i]))
                    .collect::<Vec<_>>();
                let mut quals = sorted_indices
                    .iter()
                    .map(|&i| std::mem::take(&mut quals[i]))
                    .collect::<Vec<_>>();
                
                seqs.truncate(MAX_OL_POLISHING);
                quals.truncate(MAX_OL_POLISHING);

                let cons;
                if seqs.len() == 0{
                    cons = vec![];
                }
                else if seqs.len() == 1{
                    // - 1 because last byte is null terminator
                    cons = seqs[0][0..seqs[0].len()-1].to_vec();
                }
                else{
                    cons = poa_consensus(
                        &seqs,
                        &quals,
                        consensus_max_length,
                        alignment_type,
                        match_score,
                        mismatch_score,
                        gap_open,
                        gap_extend,
                    );
                }

                //log::trace!("Consensus for block {} complete", i);

                consensuses.lock().unwrap().push((i, cons));
            });

        log::trace!("Consensus building complete for {}", self.contig_name);
        let mut consensuses = consensuses.into_inner().unwrap();
        consensuses.sort_by_key(|x| x.0);
        let consensuses = consensuses.into_iter().map(|x| x.1).collect::<Vec<_>>();
        let consensuses =
            PoaConsensusBuilder::modify_join_consensus(consensuses, self.window_overlap_len, &self.contig_name);
        return consensuses;
    }

    pub fn modify_join_consensus(mut cons: Vec<Vec<u8>>, window_len: usize, contig_name: &str) -> Vec<Vec<u8>> {
        if cons.len() == 0 {
            return vec![];
        }

        let window_cut_ratio = 2;
        let breakpoints = Mutex::new(vec![]);

        (0..cons.len() - 1).into_par_iter().for_each(|i| {
            if cons[i].len() == 0 || cons[i + 1].len() == 0 {
                let false_bp1 = 0;
                let false_bp2 = 0;
                breakpoints.lock().unwrap().push((i, false_bp1, false_bp2, 0));
                return;
            }
            let ol_len = cons[i].len().min(cons[i + 1].len().min(window_len));
            let cut = ol_len.min(window_len / window_cut_ratio);

            //          < ol_len >
            //          ALIGN CUT
            // --------|-----|-->
            //     <--|-----|---------
            //     CUT ALIGN
            let overhang_i = &cons[i][cons[i].len() - ol_len..][cut/2..ol_len - cut/2];
            let overhang_j = &cons[i + 1][0..ol_len][cut/2..ol_len-cut/2];

            //dbg!(&overhang_i, &overhang_j);

            //i is query
            let (i_end, j_end, _) = alignment::align_seq_to_ref_slice_local(overhang_i, overhang_j, &GAPS_LAX_INDEL);

            //breakpoints.push((query_pos, ref_pos));
            breakpoints.lock().unwrap().push((i, i_end, j_end, cut));
        });
        let mut breakpoints = breakpoints.into_inner().unwrap();
        breakpoints.sort_by_key(|x| x.0);
        let breakpoints = breakpoints.into_iter().map(|x| (x.1, x.2, x.3)).collect::<Vec<_>>();

        let mut new_consensus = vec![];
        let mut new_cons_i = std::mem::take(&mut cons[0]);
        for (i, bp) in breakpoints.into_iter().enumerate() {

            let mut new_cons_j = std::mem::take(&mut cons[i + 1]);
            let ol_len = new_cons_i.len().min(new_cons_j.len().min(window_len));
            let hang;

            //I believe this happens when the overlap alignments are discordant...
            if bp.0 > ol_len{
                log::debug!("Potential error in consensus joining at block {} for contig {}", i, contig_name);
                hang = ol_len;
            }
            else{
                hang = ol_len - bp.0;
            }

            let cut = bp.2;
            let break_pos_i = new_cons_i.len() - hang + cut/2;
            let break_pos_j = bp.1 + cut/2;

            new_cons_i.truncate(break_pos_i);
            new_cons_j = new_cons_j.split_off(break_pos_j);
            new_consensus.push(new_cons_i);
            new_cons_i = new_cons_j;
        }
        new_consensus.push(new_cons_i);

        return new_consensus;
    }

    pub fn new(genome_length: usize, contig_name: String) -> Self {
        PoaConsensusBuilder {
            seq: Vec::new(),
            qual: Vec::new(),
            breakpoints: Vec::new(),
            contig_name,
            genome_length,
            bp_len: 0,
            window_overlap_len: 0,
        }
    }

    fn new_test(genome_length: usize) -> Self{
        PoaConsensusBuilder {
            seq: Vec::new(),
            qual: Vec::new(),
            breakpoints: Vec::new(),
            contig_name: "test".to_string(),
            genome_length,
            bp_len: 0,
            window_overlap_len: 0,
        }
    }

    pub fn generate_breakpoints(&mut self, bp_length: usize, window_overlap_len: usize) {
        self.bp_len = bp_length;
        self.window_overlap_len = window_overlap_len;
        let mut breakpoints = Vec::new();
        let mut pos = bp_length;
        while pos < self.genome_length {
            breakpoints.push(pos);
            pos += bp_length;
        }
        breakpoints.push(self.genome_length);
        self.breakpoints = breakpoints;
        self.seq = (0..self.breakpoints.len())
            .map(|_| Mutex::new(vec![]))
            .collect();
        self.qual = (0..self.breakpoints.len())
            .map(|_| Mutex::new(vec![]))
            .collect();
    }

    fn populate_block(
        &self,
        seq: &[u8],
        qual: &[u8],
        current_block_string: &mut Vec<u8>,
        current_block_qual: &mut Vec<u8>,
        bp_i: usize,
        current_position_q: usize,
        breakpoint_q: usize,
    ) {
        current_block_string.push(0);
        current_block_qual.push(0);
        let mut seq_lock = self.seq[bp_i].lock().unwrap();
        let mut qual_lock = self.qual[bp_i].lock().unwrap();
        seq_lock.push(std::mem::take(current_block_string));
        qual_lock.push(std::mem::take(current_block_qual));
        current_block_string.extend(&seq[breakpoint_q..current_position_q]);
        current_block_qual.extend(&qual[breakpoint_q..current_position_q]);
    }

    pub fn add_seq(
        &self,
        seq: Vec<u8>,
        qual: Vec<u8>,
        cigar: &OpLenVec,
        align_start_r: usize,
        align_start_q: usize,
    ) {
        //Process CIGAR, store the aligned strings fore very breakpoint block
        let start_bp_i = self.breakpoints.iter().position(|&x| x > align_start_r);
        if start_bp_i.is_none() {
            return;
        }
        let mut bp_i = start_bp_i.unwrap();

        let mut current_block_string = vec![];
        let mut current_block_qual = vec![];
        let mut current_position_q = align_start_q;
        let mut current_position_r = align_start_r;
        let mut breakpoint = self.breakpoints[bp_i];
        let mut breakpoint_q = align_start_q;
        let mut window_breakpoint = self.breakpoints[bp_i] + self.window_overlap_len;
        let mut past_breakpoint = false;
        for (op, len) in cigar.iter() {

            match op {
                Operation::M | Operation::Eq | Operation::X => {
                    for _ in 0..len {
                        current_block_string.push(seq[current_position_q]);
                        current_block_qual.push(qual[current_position_q]);
                        current_position_q += 1;
                        current_position_r += 1;

                        if current_position_r >= breakpoint {
                            breakpoint_q = current_position_q;
                            bp_i += 1;
                            if let Some(bp) = self.breakpoints.get(bp_i) {
                                breakpoint = *bp;
                            } else {
                                return;
                            }
                            past_breakpoint = true;
                        }
                        if current_position_r >= window_breakpoint {
                            self.populate_block(
                                &seq,
                                &qual,
                                &mut current_block_string,
                                &mut current_block_qual,
                                bp_i-1,
                                current_position_q,
                                breakpoint_q,
                            );
                            window_breakpoint = breakpoint + self.window_overlap_len;
                            past_breakpoint = false;
                        }
                    }
                }
                Operation::I => {
                    for _ in 0..len {
                        current_block_string.push(seq[current_position_q]);
                        current_block_qual.push(qual[current_position_q]);
                        current_position_q += 1;
                    }
                }
                Operation::D => {
                    for _ in 0..len {
                        current_position_r += 1;
                        if current_position_r >= breakpoint {
                            breakpoint_q = current_position_q;
                            bp_i += 1;
                            if let Some(bp) = self.breakpoints.get(bp_i) {
                                breakpoint = *bp;
                            } else {
                                return;
                            }
                            past_breakpoint = true;
                        }
                        if current_position_r >= window_breakpoint {
                            self.populate_block(
                                &seq,
                                &qual,
                                &mut current_block_string,
                                &mut current_block_qual,
                                bp_i - 1,
                                current_position_q,
                                breakpoint_q,
                            );
                            window_breakpoint = breakpoint + self.window_overlap_len;
                            past_breakpoint = false;
                        }
                    }
                }
                _ => {}
            }
        }

        let last_bp_i = if past_breakpoint { bp_i - 1 } else { bp_i };
        if current_block_qual.len() > 0 {
            self.populate_block(
                &seq,
                &qual,
                &mut current_block_string,
                &mut current_block_qual,
                last_bp_i,
                current_position_q,
                current_position_q
            );
        }
    }

    pub fn process_mapping_boundaries(
        &mut self,
        mapping_boundaries: &[&Interval<u32, SmallTwinOl>],
        twin_reads: &[TwinRead],
    ) {
        mapping_boundaries.par_iter().for_each(|interval| {
            let ol = &interval.val;
            let ar = &ol.alignment_result.as_ref().unwrap();
            //TODO assume qualities exist
            let query_seq = &twin_reads[ol.query_id as usize].dna_seq;
            let query_quals = &twin_reads[ol.query_id as usize].qual_seq.as_ref().unwrap();
            let query_seq_u8: Vec<u8>;
            let query_quals_u8: Vec<u8>;
            if ol.reverse {
                query_seq_u8 = query_seq
                    .to_revcomp()
                    .iter()
                    .map(|x| x.to_char().to_ascii_uppercase() as u8)
                    .collect();
                query_quals_u8 = query_quals
                    .to_revcomp()
                    .iter()
                    .map(|x| (x as u8) * 3 + 33)
                    .collect();
            } else {
                query_seq_u8 = query_seq
                    .iter()
                    .map(|x| x.to_char().to_ascii_uppercase() as u8)
                    .collect();
                query_quals_u8 = query_quals.iter().map(|x| (x as u8) * 3 + 33).collect();
            }
            self.add_seq(
                query_seq_u8,
                query_quals_u8,
                &ar.cigar,
                ar.r_start,
                ar.q_start,
            );
        })
    }
}

pub fn join_circular_ends(seq: &mut Vec<u8>, overlap_len: usize, hang1: usize, hang2: usize, contig_name: &str, args: &Cli) {

    let seq_len = seq.len();
    if seq_len < 1000{
        return
    }

    //let overlap_length = edge.overlap.overlap_len_bases + edge.overlap.hang1 + edge.overlap.hang2 + 1000;
    let overlap_length = overlap_len + hang1 + hang2 + 500;

    if overlap_length > seq_len{
        log::warn!("Overlap length is greater than sequence length; something went wrong during polishing.");
        return
    }

    let overhang1 = seq[seq_len-overlap_length..seq_len].to_vec();
    let overhang2 = seq[0..overlap_length].to_vec();

    let tr1 = seeding::get_twin_read_syncmer(overhang1, None, args.kmer_size, args.c, &FxHashSet::default(), String::new()).unwrap();
    let tr2 = seeding::get_twin_read_syncmer(overhang2, None, args.kmer_size, args.c, &FxHashSet::default(), String::new()).unwrap();

    let mut lax_args = args.clone();
    lax_args.min_ol = args.min_ol/2;

    let tr_options = CompareTwinReadOptions{
        compare_snpmers: false,
        retain_chain: false,
        force_one_to_one_alignments: true,
        ..Default::default()
    };

    let overlaps = mapping::compare_twin_reads(&tr1, &tr2, None, None, 0, 1, &tr_options, args);

    if overlaps.is_empty(){
        log::debug!("Circular contig with hash id {} was not able to be end-polished correctly", &contig_name);
        return
    }

    let best_overlap = overlaps.iter().max_by_key(|x| x.end1 - x.start1).unwrap();

    //dbg!(&best_overlap);
    let end_trim = overlap_length - best_overlap.start1;
    let start_trim = best_overlap.start2;
    let new_seq = seq[start_trim..seq_len-end_trim].to_vec();
    *seq = new_seq;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_breakpoints() {
        let mut builder = PoaConsensusBuilder::new_test(100);
        builder.generate_breakpoints(10, 0);
        assert_eq!(
            builder.breakpoints,
            vec![10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
        );
    }

    #[test]
    fn test_add_seq_1() {
        let mut builder = PoaConsensusBuilder::new_test(10);
        builder.generate_breakpoints(10, 0);
        let seq = b"ACGTACGTACGT".to_vec();
        let qual = vec![30; 12];
        let cigar = vec![
            OpLen {
                op: Operation::M,
                len: 4,
            },
            OpLen {
                op: Operation::I,
                len: 4,
            },
            OpLen {
                op: Operation::M,
                len: 4,
            },
        ];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);
        assert_eq!(builder.seq.len(), 1);
        assert_eq!(builder.qual.len(), 1);
        assert_eq!(
            builder.seq[0].lock().unwrap()[0],
            b"ACGTACGTACGT\0".to_vec()
        );
    }

    #[test]
    fn test_add_seq_2() {
        let mut builder = PoaConsensusBuilder::new_test(100);
        builder.generate_breakpoints(10, 0);
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTA".to_vec();
        let qual = vec![30; 29];
        let cigar = vec![OpLen {
            op: Operation::M,
            len: 29,
        }];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);
        assert_eq!(builder.seq.len(), 10);
        assert_eq!(builder.qual.len(), 10);
        assert_eq!(builder.seq[0].lock().unwrap()[0], b"ACGTACGTAC\0".to_vec());
        assert_eq!(builder.seq[1].lock().unwrap()[0], b"GTACGTACGT\0".to_vec());
        assert_eq!(builder.seq[2].lock().unwrap()[0], b"ACGTACGTA\0".to_vec());
    }

    #[test]
    fn test_add_seq_2_window() {
        let mut builder = PoaConsensusBuilder::new_test(100);
        builder.generate_breakpoints(10, 2);
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTA".to_vec();
        let qual = vec![30; 29];
        let cigar = vec![OpLen {
            op: Operation::M,
            len: 29,
        }];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);
        assert_eq!(builder.seq.len(), 10);
        assert_eq!(builder.qual.len(), 10);
        assert_eq!(
            builder.seq[0].lock().unwrap()[0],
            b"ACGTACGTACGT\0".to_vec()
        );
        assert_eq!(
            builder.seq[1].lock().unwrap()[0],
            b"GTACGTACGTAC\0".to_vec()
        );
        assert_eq!(builder.seq[2].lock().unwrap()[0], b"ACGTACGTA\0".to_vec());
    }

    #[test]
    fn test_add_seq_3() {
        let mut builder = PoaConsensusBuilder::new_test(100);
        builder.generate_breakpoints(10, 0);
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec();
        let qual = vec![30; 30];
        let cigar = vec![OpLen {
            op: Operation::M,
            len: 30,
        }];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);
        assert_eq!(builder.seq.len(), 10);
        assert_eq!(builder.qual.len(), 10);
        assert_eq!(builder.seq[0].lock().unwrap()[0], b"ACGTACGTAC\0".to_vec());
        assert_eq!(builder.seq[1].lock().unwrap()[0], b"GTACGTACGT\0".to_vec());
        assert_eq!(builder.seq[2].lock().unwrap()[0], b"ACGTACGTAC\0".to_vec());

        //Only 3 of the 10 blocks are filled
        assert_eq!(builder.seq[3].lock().unwrap().len(), 0);
        assert_eq!(builder.seq[4].lock().unwrap().len(), 0);
    }

    #[test]
    fn test_add_seq_del() {
        let mut builder = PoaConsensusBuilder::new_test(100);
        builder.generate_breakpoints(10, 0);
        let seq = b"ACCGTACGTACGTACGTACGTACGTAC".to_vec();
        let qual = vec![30; 30];
        let cigar = vec![
            OpLen {
                op: Operation::M,
                len: 2,
            },
            OpLen {
                op: Operation::D,
                len: 3,
            },
            OpLen {
                op: Operation::M,
                len: 25,
            },
        ];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);
        assert_eq!(builder.seq.len(), 10);
        assert_eq!(builder.qual.len(), 10);
        assert_eq!(builder.seq[0].lock().unwrap()[0], b"ACCGTAC\0".to_vec());
        assert_eq!(builder.seq[1].lock().unwrap()[0], b"GTACGTACGT\0".to_vec());
        assert_eq!(builder.seq[2].lock().unwrap()[0], b"ACGTACGTAC\0".to_vec());

        //Only 3 of the 10 blocks are filled
        assert_eq!(builder.seq[3].lock().unwrap().len(), 0);
        assert_eq!(builder.seq[4].lock().unwrap().len(), 0);
    }

    #[test]
    fn test_add_seq_ins() {
        let mut builder = PoaConsensusBuilder::new_test(100);
        builder.generate_breakpoints(10, 0);
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTATTC".to_vec();
        let qual = vec![30; 32];
        let cigar = vec![
            OpLen {
                op: Operation::M,
                len: 29,
            },
            OpLen {
                op: Operation::I,
                len: 2,
            },
            OpLen {
                op: Operation::M,
                len: 1,
            },
        ];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);
        assert_eq!(builder.seq.len(), 10);
        assert_eq!(builder.qual.len(), 10);
        assert_eq!(builder.seq[0].lock().unwrap()[0], b"ACGTACGTAC\0".to_vec());
        assert_eq!(builder.seq[1].lock().unwrap()[0], b"GTACGTACGT\0".to_vec());
        assert_eq!(
            builder.seq[2].lock().unwrap()[0],
            b"ACGTACGTATTC\0".to_vec()
        );

        //Only 3 of the 10 blocks are filled
        assert_eq!(builder.seq[3].lock().unwrap().len(), 0);
        assert_eq!(builder.seq[4].lock().unwrap().len(), 0);
    }

    #[test]
    fn test_poa_cons() {
        for window_len in [0, 0] {
            let mut builder = PoaConsensusBuilder::new_test(100);
            builder.generate_breakpoints(10, window_len);
            let seq = b"ACGTACGTACGTACGTACGTACGTACGTATTC".to_vec();
            let qual = vec![30; 32];
            let cigar = vec![
                OpLen {
                    op: Operation::M,
                    len: 29,
                },
                OpLen {
                    op: Operation::I,
                    len: 2,
                },
                OpLen {
                    op: Operation::M,
                    len: 1,
                },
            ];
            builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);

            let seq = b"ACGTACGTACGTACGTACGTACGTACGTATTC".to_vec();
            let qual = vec![30; 32];
            let cigar = vec![
                OpLen {
                    op: Operation::M,
                    len: 29,
                },
                OpLen {
                    op: Operation::I,
                    len: 2,
                },
                OpLen {
                    op: Operation::M,
                    len: 1,
                },
            ];
            builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);

            let seq = b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec();
            let qual = vec![30; 32];
            let cigar = vec![
                OpLen {
                    op: Operation::M,
                    len: 29,
                },
                OpLen {
                    op: Operation::M,
                    len: 1,
                },
            ];
            builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);

            dbg!(&"HERE");

            let consensuses = builder.spoa_blocks();
            assert_eq!(consensuses.len(), 10);
            assert_eq!(consensuses[0], b"ACGTACGTAC".to_vec());
            assert_eq!(consensuses[1], b"GTACGTACGT".to_vec());
            assert_eq!(consensuses[2], b"ACGTACGTATTC".to_vec());
            assert_eq!(consensuses[3], b"".to_vec());
        }
    }

    #[test]
    fn test_poa_triplets() {
        let mut builder = PoaConsensusBuilder::new_test(100);
        builder.generate_breakpoints(10, 0);

        let seq = b"ATCG".to_vec();
        let qual = vec![50; 3];
        let cigar = vec![OpLen {
            op: Operation::M,
            len: 3,
        }];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);

        let seq = b"ATG".to_vec();
        let qual = vec![50; 3];
        let cigar = vec![OpLen {
            op: Operation::M,
            len: 3,
        }];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);

        let seq = b"ATG".to_vec();
        let qual = vec![50; 3];
        let cigar = vec![OpLen {
            op: Operation::M,
            len: 3,
        }];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);
        let cons = builder.spoa_blocks();
        dbg!(&cons[0]);
    }

    #[test]
    fn test_poa_with_without_window() {
        let window_lengths = [0, 3];
        for window_length in window_lengths {
            let mut builder = PoaConsensusBuilder::new_test(100);
            builder.generate_breakpoints(5, window_length);

            let seq = b"CCCCCTTTTTGGGGGAAAAA".to_vec();
            let qual = vec![50; 20];
            let cigar = vec![OpLen {
                op: Operation::M,
                len: 20,
            }];
            builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);

            let seq = b"CCCCCATTTTTGGGGGAAAAA".to_vec();
            let qual = vec![50; 21];
            let cigar = vec![
                OpLen {
                    op: Operation::M,
                    len: 3,
                },
                OpLen {
                    op: Operation::I,
                    len: 1,
                },
                OpLen {
                    op: Operation::M,
                    len: 7,
                },
                OpLen {
                    op: Operation::M,
                    len: 10,
                },
            ];
            builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);

            let seq = b"CCCCCTTTTTGGGGGAAAAA".to_vec();
            let qual = vec![50; 20];
            let cigar = vec![OpLen {
                op: Operation::M,
                len: 20,
            }];
            builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);
            //dbg!(&builder.seq[1].lock().unwrap());
            let cons = builder.spoa_blocks();
            for cons in cons.iter() {
                if cons.len() != 0{
                    println!("{:?}", String::from_utf8_lossy(&cons));
                }
            }

            if window_length == 0 {
                //insertion at border
                assert!(cons[0].len() == 6);
            } else {
                //no insertion
                assert!(cons.iter().map(|x| x.len()).sum::<usize>() == 20);
            }
        }
    }

    #[test]
    fn test_poa_cons_quality() {
        let mut builder = PoaConsensusBuilder::new_test(100);

        builder.generate_breakpoints(10, 1);

        let seq = b"ACGTTTACGTACGTACGTACGTACGTACGTAC".to_vec();
        let qual = vec![60; 32];
        let cigar = vec![
            OpLen {
                op: Operation::M,
                len: 3,
            },
            OpLen {
                op: Operation::I,
                len: 2,
            },
            OpLen {
                op: Operation::M,
                len: 27,
            },
        ];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);

        let seq = b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec();
        let qual = vec![50; 32];
        let cigar = vec![
            OpLen {
                op: Operation::M,
                len: 29,
            },
            OpLen {
                op: Operation::M,
                len: 1,
            },
        ];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);

        let seq = b"ACGTACGTACGTACGTACGTACGTACGTAC".to_vec();
        let qual = vec![50; 32];
        let cigar = vec![
            OpLen {
                op: Operation::M,
                len: 29,
            },
            OpLen {
                op: Operation::M,
                len: 1,
            },
        ];
        builder.add_seq(seq, qual, &OpLenVec::new(cigar), 0, 0);

        for x in builder.seq[2].lock().unwrap().iter() {
            println!("{:?}", x.to_ascii_uppercase());
        }
        let consensuses = builder.spoa_blocks();
        assert_eq!(consensuses.len(), 10);
        // assert_eq!(consensuses[0], b"ACGTACGTAC".to_vec());
        // assert_eq!(consensuses[1], b"GTACGTACGT".to_vec());
        // assert_eq!(consensuses[2], b"ACGTACGTAC".to_vec());
        // assert_eq!(consensuses[3], b"".to_vec());
        assert_eq!(consensuses.iter().map(|x| x.len()).sum::<usize>(), 30);
    }

    #[test]
    fn test_dna_consensus() {
        let mut seqs = vec![];
        let mut quals = vec![];

        // generated each string by adding small tweaks to the expected consensus "AATGCCCGTT"
        for seq in [
            "ATTGCCCGTT\0",
            "AATGCCGTT\0",
            "AATGCCCGAT\0",
            "AACGCCCGTC\0",
            "AGTGCTCGTT\0",
            "AATGCTCGTT\0",
        ]
        .iter()
        {
            seqs.push((*seq).bytes().map(|x| x as u8).collect::<Vec<u8>>());
        }

        //generate quality scores
        for qual in vec![
            "1111111111\0",
            "111111111\0",
            "1111111111\0",
            "1111111111\0",
            "1111111111\0",
            "1111111111\0",
        ]
        .iter()
        {
            quals.push((*qual).bytes().map(|x| x as u8).collect::<Vec<u8>>());
        }

        let consensus = poa_consensus(&seqs, &quals, 20, 1, 5, -4, -3, -1);

        let expected = "AATGCCCGTT".to_string().into_bytes();
        assert_eq!(consensus, expected);
    }

    #[test]
    fn test_poa_cons_real_quality() {
        let mut builder = PoaConsensusBuilder::new_test(10000);

        builder.generate_breakpoints(700, 50);

        let seqs: Vec<Vec<u8>> = vec![b"CTGATAACTCTAAATCTTGCTGTCAGGCAT".to_vec(),
        b"CTGATAACTCAAACTGCTGTCAGGCATTGCCAGAACAGCAAGATAATGAATGCCAAATTCATCTATCTTATTGAATACGATGTCACTCAATATAGACTGGTCCAGGAAAGAAGAAGGAACCTGCTTATCGAAATCCTTCTTATGAAAATCACGAGAGAACAGACAGTTGGCACCATGATAAACATGCAAGTTCACAATATTATCATAATAAACATTGTCTACCTCAACACCATCATCATTATAAGAACTCTTGTAGACCTTATAACTTGTTGGGTTAACCTGTACATAGAGATGATATTTCCCATCGCCCCTTACAACAACCGTATCACGTTTAATCAGCGTATTCTGATTCAGCGCCACGGAAGCATGGTTGTGGATGAACTGCTGCAAATAAGACTTATCTTCCGTCTTCACCAGCTTAACCACATCCCCATTCTGCACTTTAAACTGGAAAAGTACAATGCTCAGCCTGTTTTACGATCGGATATTTTGCAATATTAGC".to_vec(),
        b"CTGATAACTCAAACTGCTGTCAGGCATTGCCAGAACAGCAAGATAATGAATGCCAAATTCATCTATCTTATTGAATACGATGTCACTCAATATAGACTGGTCCAGGAAAGAAGAAGGAACCTGCTTATCGAAATCCTTCTTATGAAAATCACGAGAAGAAACAAGCAGTTGGCACCATGATAAACATGCAAGTTCACAATATTATCATAATAAACATTGTCTACCTCAACACCATCATCATTATAAGAACTCTTGTAGACCTTATAACTTGTTGGGTTAACCTGTACATAGAGATGATATTTCCCATCGCCCCTTACAACAACCGTATCACGTTTAATCAGCGTATTCTGATTCAGCGCCACGGAAGCATGGTTGTGGATGAACTGCTGCAAATAAGACTTATCTTCCGTCTTCACCAGCTTAACCACATCCCCATTCTGCACTTTAAACTGGAAAATATGCTCAGCCTGTTTTACGATCGGATATTTTGCAATATTAGC".to_vec(),
        b"CTGATAACTCAAACTGCTGTCAGGCATTGCCAGAACAGCAAGATAATGAATGCCAAATTCATCTCTGTTGCTTGAATACGATGTCACTCAATATAGACTGGTCCAGAAAGAAGAAGGAACCTGCTTATCGAAATCCTTCTTATGAAAATCACAAGAAAGAGAACAGACAGTTGGCACCATGATAAACATGCAAGTTCACAATATTATCATAATAAACATTGTCTACCTCAACACCATCATCATTATAAGAACTCTTGTAGACCTTATAACTTGTTGGGTTAACCTGTACATAGAGATGATATTTCCCATCGCCCCTTACAACAACCGTATCACGTTTAATCAGCGTATTCTGATTCAGCGCCACGGAAGCATGGTTGTGGATGAACTGCTGCAAATAAGACTTATCTTCCGTCTTCACCAGCTTAACCACATCCCCATTCTGCACTTTAAACTGAAAAATATTCCTCAGCCTGTTTTACGATCGGATATTTTGCAATATTAGC".to_vec(),
        b"CTGATAACTCAAACTGCTGTCAGGCATTGCCAGAACAGCAAGATAATGAATGCCAAATTCATCTATCTTATTGAATACGATGTCACTCAATATAGACTGGTCCAGGAAAGAAGAAGGAACCTGCTTATCGAAATCCTTCTTATGAAAATCACGAGAGAACAGACAGTTGGCACCATGATAAACATGCAAGTTCACAATATTATCATAATAAACATTGTCTACCTCAACACCATCATCATTATAAGAACTCTTGTAGACCTTATAACTTGTTGGGTTAACCTGTACATAGAGATGATATTTCCCATCGCCCCTTACAACAACCGTATCACGTTTAATCAGCGTATTCTGATTCAGCGCCACGGAAGCATGGTTGTGGATGAACTGCTGCAAATAAGACTTATCTTCCGTCTTCACCAGCTTAACCACATCCCCATTCTGCACTTTAAACTGGAAAATATGCTCAGCCTGTTTTACGATCGGATATTTTGCAATATTAGC".to_vec(),
        b"CATTGCCAGAACAGCAAGATAATGAATGCCAAATTCATCTATCTTATTGAATACGATGTCACTCAATATAGACTGGTCCAGGAAAGAAGAAGGAACCTGCTTATCGAAATCCTTCTTATGAAAATCACGAGAGAACAGACAGTTGGCACCATGATAAACATGCAAGTTCACAATATTATCATAATAAACATTGTCTACCTCAACACCATCATCATTATCAACTCTTGTAGACCTTATAACTTGTTGGGTTAACCTGTACATAGAGATGATATTTCCCATCGCCCCTTACAACAACCGTATCACGTTTAATCAGCGTATTCTGATTTCAGCGCCACAAAGTCCTGGTTGTGGATGAACTGCTGCAAATAAGACTTATCTTCCGTCTTCACCAGCTTAACCACATCCCCCATTCTGCACTTACAACTGGAAAATATGCTCAGCCTGTTTTACGATCGGATATTTTGCAATATTAGC".to_vec(),];
        let mut quals = seqs.iter().map(|x| vec![50; x.len()]).collect::<Vec<_>>();
        quals[0] = vec![100; 30];

        let cigars = seqs
            .iter()
            .map(|x| {
                vec![OpLen {
                    op: Operation::M,
                    len: x.len(),
                }]
            })
            .collect::<Vec<_>>();

        for i in 0..6 {
            builder.add_seq(seqs[i].clone(), quals[i].clone(), &OpLenVec::new(cigars[i].clone()), 0, 0);
        }

        let consensuses = builder.spoa_blocks();
        println!("Consensuses: {:?}", consensuses[0]);
        assert!(consensuses[0].len() > 480 && consensuses[0].len() < 520);
    }
    
    #[test]
    fn test_modify_new(){
        let consensuses = vec![
            b"GCAACGTATGT".to_vec(), // last C is error
            b"ACGTATGTGTGTGT".to_vec() // first T is errors
        ];

        let new_cons = PoaConsensusBuilder::modify_join_consensus(consensuses, 8, "test");

        let mut final_consensus = vec![];
        for cons in new_cons{
            println!("{:?}", String::from_utf8_lossy(&cons));
            final_consensus.extend(cons);
        }

        assert_eq!(final_consensus, b"GCAACGTATGTGTGTGT".to_vec());

    }

    #[test]
    fn test_modify_new_del(){
        let consensuses = vec![
            b"AAAAATTTA".to_vec(), // last A is error
            b"CTTTGGGG".to_vec() // first T is errors
        ];

        let new_cons = PoaConsensusBuilder::modify_join_consensus(consensuses, 5, "test");

        let mut final_consensus = vec![];
        for cons in new_cons{
            println!("{:?}", String::from_utf8_lossy(&cons));
            final_consensus.extend(cons);
        }

        assert_eq!(final_consensus, b"AAAAATTTGGGG".to_vec());

    }

    #[test]
    fn circular_join_basic_test(){
        let mut args = Cli::default();
        args.c = 5;
        args.kmer_size = 10;
        {
            let random_string = b"GCATGCGTTCAACGTAGGCCGTACTAGCTGCGTAATCGACGGAATGGCAGTATCGCGATAACGCTTGAAACGCTACGAGCCATAGCGGTATCGTAGCAACGCTAATCGGCATAGCTATCGATGCAGTCGCTATAGCTAGCTAGCGATCGGCCGATAGCGATCGATCGGCTAGCGGCATCGATAGCGGCCGATCGCGATCAGCATGGCCGATGCGATCGCGTATCAGCGCGATCGAGCCGATCGATCGCGTCCGATGCATGCAACGATCGGCATATCACGCGCGATCGACTAGCGATCGATCGCGTACGCATCGATCGAGCGATCGACTGATCGCTAGCTGCATGCATACGCTAGCTGCAGCTAGCATCGATCGCTATGCTAGCTAGCATCGAGCTGATCGTAGCATCGATCGATCGATCGATCGATCGAGCTATCGATCGATACGCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGCGATCGACTGCGATCGCTAGCTAGCTAGCTATGCTAGCTAGCTGCTAGTCGACGATCGATCGATCGATCGATCTAGCTAGCATCGCTAGCTGATCGTAGCTAGCTAGCATCGATCGA".to_vec();
            let mut seq = random_string.clone();
            seq.extend(vec![b'G'; 250]);
            seq.extend(vec![b'A'; 250]);
            seq.extend(vec![b'T'; 250]);
            seq.extend(vec![b'C'; 250]);
            seq.extend(random_string.clone());

            dbg!(seq.len());

            join_circular_ends(&mut seq, random_string.len(), 10, 10, "test", &args);

            assert_eq!(seq.len(), 1000 + random_string.len());
        }
    }
}
