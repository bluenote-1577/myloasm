use crate::cli::Cli;
use crate::constants::*;
use crate::mapping;
use crate::polishing::alignment;
use crate::seeding;
use crate::seeding::estimate_sequence_identity;
use crate::types::SmallTwinOl;
use crate::types::*;
use abpoa::{AlignMode, Aligner, OutputMode, Parameters, Scoring, SequenceBatch};
use bio_seq::prelude::*;
use block_aligner::cigar::*;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use rayon::prelude::*;
use rust_lapper::Interval;
use rust_spoa::poa_consensus;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;
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
    genome: Vec<u8>,
    args: Cli,
}

impl PoaConsensusBuilder {
    pub fn spoa_blocks(self) -> Vec<Vec<u8>> {
        let consensus_max_length = 1500;
        let alignment_type = 1; // 0 local, 1 global, 2, semiglobal
        let match_score = 5;
        let mismatch_score = -2;
        let gap_open = -2;
        let gap_extend = -1;

        let consensuses = Mutex::new(vec![]);

        // TODO DEBUG
        // Open bufwriter for writing for this contig
        // let file = File::create("consensus_latest.txt").unwrap();
        // let writer = Mutex::new(BufWriter::new(file));

        // let file = File::create("complete_consensus.txt").unwrap();
        // let complete_writer = Mutex::new(BufWriter::new(file));

        //parallel iter over seq and qual
        self.seq
            .into_par_iter()
            .enumerate()
            .zip(self.qual.into_par_iter())
            .for_each(|((i, seq), qual)| {
                let mut seqs = seq.lock().unwrap();
                let mut quals = qual.lock().unwrap();
                let debug = i % 100 == 0 && false;
                let start_time_everything = std::time::Instant::now();

                assert!(seqs.len() == quals.len());
                let preprocess_time = std::time::Instant::now();

                let max_len = seqs.iter().map(|x| x.len()).max().unwrap_or(0).min(self.window_overlap_len + self.bp_len);

                if seqs.len() <= 3 {
                    seqs.retain(|y| y.len() > max_len * 4 / 5);
                    quals.retain(|y| y.len() > max_len * 4 / 5);
                }

                //log::trace!("Processing block {} of {}: {:?}",i, self.contig_name, seqs.iter().map(|x| x.len()).collect::<Vec<_>>());

                let mut qual_map = FxHashMap::default();
                let mut length_dist = vec![];
                for (i, qual) in quals.iter().enumerate() {
                    let mean_qual = estimate_sequence_identity(Some(&qual[0..qual.len()-1])).unwrap_or(0.0);

                    if mean_qual < self.args.min_qual_polishing {
                        continue;
                    }

                    // IMPORTANT!!!!
                    // If the entire string is qual 0, it STALLS the polishing step (spoa does not finish)
                    if qual[0..qual.len()-1].iter().all(|&q| q == 33) {
                        continue;
                    }

                    //let mean_qual = qual.iter().map(|x| (*x - 33)as f64).sum::<f64>() / qual.len() as f64;
                    length_dist.push(qual.len());
                    //qual_map.insert(i, -(qual.len() as i32) * mean_qual);
                    qual_map.insert(i, (mean_qual) as i32);
                }
                length_dist.sort();
                let median_length_dist;
                let twenty_five_length;
                let seventy_five_length;
                if qual_map.len() == 0{
                    median_length_dist = self.window_overlap_len + self.bp_len;
                    twenty_five_length = self.window_overlap_len + self.bp_len;
                    seventy_five_length = self.window_overlap_len + self.bp_len;
                } 
                else{
                    median_length_dist = length_dist[length_dist.len() / 2];
                    twenty_five_length = length_dist[length_dist.len() / 4];
                    seventy_five_length = length_dist[length_dist.len() * 3 / 4];
                }
                log::trace!("Median length distribution: {} ({}-{}) with {} seqs at block {}", median_length_dist, twenty_five_length, seventy_five_length, qual_map.len(), i);

                for (i, mean_qual) in qual_map.iter_mut() {
                    //De-prioritize reads that are very different in length from the median
                    let length_penalty = (seventy_five_length as i32 - quals[*i].len() as i32).abs();
                    *mean_qual -= length_penalty;
                }

                //Sort seqs and quals by the qual_map
                let mut sorted_indices = qual_map.into_iter().collect::<Vec<_>>();
                sorted_indices.sort_by_key(|x| -x.1);
                let mut sorted_indices = sorted_indices.into_iter().map(|x| x.0).collect::<Vec<_>>();
                if seqs.len() > 20{
                    let trim = sorted_indices.len()  / 20;
                    if sorted_indices.len() > 10{
                        sorted_indices = sorted_indices[..sorted_indices.len() - trim].to_vec();
                    }
                    else{
                        log::debug!("Many low quality sequences ({}) for block {} of contig {}, proceeding with available sequences", seqs.len() - sorted_indices.len(), i, self.contig_name);
                    }
                }

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

                //Write seqs and quals to consensus_latest.txt for debugging
                // let writer = &mut writer.lock().unwrap();
                // writeln!(writer, ">Block {} of contig {} with {} seqs", i, self.contig_name, seqs.len()).unwrap();
                // for seq in seqs.iter() {
                //     writeln!(writer, "{}", std::str::from_utf8(seq).unwrap()).unwrap();
                // }
                // for qual in quals.iter() {
                //     writeln!(writer, "{}", std::str::from_utf8(qual).unwrap()).unwrap();
                // }

                let mut max_run = 0;
                if seqs.len() > 0{
                    let mut current_base = seqs[0][0];
                    let mut current_run = 0;
                    for base in seqs[0].iter() {
                        if *base == current_base {
                            current_run += 1;
                        }
                        else{
                            current_base = *base;
                            current_run = 1;
                        }
                        if current_run > max_run {
                            max_run = current_run;
                        }
                    }
                }

                if debug {
                    log::debug!("Processing time for block {} of contig {} with {} seqs and max homopolymer run {}: {:?}", 
                    i, self.contig_name, seqs.len(), max_run, preprocess_time.elapsed());
                }

                let use_hpc = if max_run >= 10 { true } else { false };
                let use_abpoa = seventy_five_length - twenty_five_length < 100;

                let mut cons;
                if seqs.len() == 0{
                    let start = i * self.bp_len;
                    let end = (i+1) * self.bp_len + self.window_overlap_len;
                    let end = end.min(self.genome.len());
                    if start > end{
                        cons = vec![];
                    }
                    else{
                        cons = self.genome[start..end].to_vec();
                    }
                }
                else if seqs.len() == 1{
                    // - 1 because last byte is null terminator
                    cons = seqs[0][0..seqs[0].len()-1].to_vec();
                }
                else if use_hpc || self.args.hpc {
                    let hpc = hpc_compress_all(&seqs, &quals);

                    let start_time = std::time::Instant::now();
                    cons = if use_abpoa{
                        abpoa_consensus_impl(&hpc.seqs, &hpc.quals, match_score, mismatch_score, gap_open, gap_extend)
                    } else {
                        poa_consensus(&hpc.seqs, &hpc.quals, consensus_max_length, alignment_type, match_score, mismatch_score, gap_open, gap_extend)
                    };
                    if cons.is_empty() && !hpc.seqs.is_empty() {
                        cons = hpc.seqs[0][..hpc.seqs[0].len() - 1].to_vec();
                    }
                    if debug {
                        log::debug!("POA consensus time for block {} of contig {} with {} HPC seqs: {:?}",
                        i, self.contig_name, hpc.seqs.len(), start_time.elapsed());
                    }

                    // Trim on the HPC consensus, then expand
                    let start_time = std::time::Instant::now();
                    if quals.len() > 15 && self.args.new_polish_trimming {
                        cons = trim_consensus_by_coverage(&cons, &hpc.seqs);
                    }
                    if debug {
                        log::debug!("Trimming time for block {} of contig {} with HPC consensus length {}: {:?}", 
                        i, self.contig_name, cons.len(), start_time.elapsed());
                    }

                    let start_time = std::time::Instant::now();
                    cons = expand_hpc_consensus(&cons, &hpc);
                    if debug {
                        log::debug!("HPC expansion time for block {} of contig {} with HPC consensus length {}: {:?}", 
                        i, self.contig_name, cons.len(), start_time.elapsed());
                    }
                }
                else{
                    let start_time = std::time::Instant::now();
                    cons = if use_abpoa{
                        abpoa_consensus_impl(&seqs, &quals, match_score, mismatch_score, gap_open, gap_extend)
                    } else {
                        poa_consensus(&seqs, &quals, consensus_max_length, alignment_type, match_score, mismatch_score, gap_open, gap_extend)
                    };
                    if cons.is_empty() && !seqs.is_empty() {
                        cons = seqs[0][..seqs[0].len() - 1].to_vec();
                    }
                    if debug{
                        log::debug!("POA consensus time for block {} of contig {} with {} seqs: {:?}",
                        i, self.contig_name, seqs.len(), start_time.elapsed());
                    }

                    // Trim consensus ends based on coverage
                    if quals.len() > 15 && self.args.new_polish_trimming{
                        let start_time = std::time::Instant::now();
                        cons = trim_consensus_by_coverage(&cons, &seqs);
                        if debug{
                            log::debug!("Trimming time for block {} of contig {} with consensus length {}: {:?}", 
                            i, self.contig_name, cons.len(), start_time.elapsed());
                        }
                    }
                }

                if cons.len() as i32 - seventy_five_length as i32 > 300 && quals.len() >= 5 && (cons.len() as i32 - (self.bp_len + self.window_overlap_len) as i32) > 150{
                    log::debug!("Warning: Consensus for block at approximately {} of {} is much longer than expected ({} vs {}). This may indicate low coverage or poor consensus.", 
                    i * (self.bp_len), self.contig_name, cons.len(), seventy_five_length);
                }

                //log::trace!("Consensus for block {} complete", i);

                if debug{
                    log::debug!("Total processing time for block {} of contig {} with {} seqs: {:?}", 
                    i, self.contig_name, seqs.len(), start_time_everything.elapsed());
                }
                consensuses.lock().unwrap().push((i, cons));

                //writeln!(complete_writer.lock().unwrap(), "Block {} completed for {} with {} seqs", i, self.contig_name, seqs.len()).unwrap();
            });

        log::trace!("Consensus building complete for {}", self.contig_name);
        let mut consensuses = consensuses.into_inner().unwrap();
        consensuses.sort_by_key(|x| x.0);
        let consensuses = consensuses.into_iter().map(|x| x.1).collect::<Vec<_>>();
        // if log trace level
        if log::log_enabled!(log::Level::Trace) {
            //mkdir cons/
            let cons_path = Path::new(&self.args.output_dir).join("cons");
            std::fs::create_dir_all(&cons_path).unwrap();

            let filename = cons_path.join(format!("cons_test_{}.fa", self.contig_name));
            let mut fasta_writer = BufWriter::new(File::create(&filename).unwrap());
            for (i, cons) in consensuses.iter().enumerate() {
                writeln!(fasta_writer, ">contig_{}_block_{}", self.contig_name, i).unwrap();
                writeln!(fasta_writer, "{}", std::str::from_utf8(cons).unwrap()).unwrap();
            }
        }

        let consensuses = PoaConsensusBuilder::modify_join_consensus(
            consensuses,
            self.window_overlap_len,
            self.bp_len,
            &self.contig_name,
        );

        if log::log_enabled!(log::Level::Trace) {
            //mkdir cons/
            let cons_path = Path::new(&self.args.output_dir).join("cons");
            std::fs::create_dir_all(&cons_path).unwrap();

            let filename = cons_path.join(format!("cons_test_{}_polished.fa", self.contig_name));
            let mut fasta_writer = BufWriter::new(File::create(&filename).unwrap());
            for (i, cons) in consensuses.iter().enumerate() {
                writeln!(fasta_writer, ">contig_{}_block_{}", self.contig_name, i).unwrap();
                writeln!(fasta_writer, "{}", std::str::from_utf8(cons).unwrap()).unwrap();
            }
        }

        return consensuses;
    }

    pub fn modify_join_consensus(
        mut cons: Vec<Vec<u8>>,
        window_len: usize,
        bp_length: usize,
        contig_name: &str,
    ) -> Vec<Vec<u8>> {
        if cons.len() == 0 {
            return vec![];
        }

        let window_cut_ratio = 1000000;
        let breakpoints = Mutex::new(vec![]);

        (0..cons.len() - 1).into_par_iter().for_each(|i| {
            if cons[i].len() == 0 || cons[i + 1].len() == 0 {
                let false_bp1 = 0;
                let false_bp2 = 0;
                breakpoints
                    .lock()
                    .unwrap()
                    .push((i, false_bp1, false_bp2, 0));
                return;
            }
            let ol_len = cons[i].len().min(cons[i + 1].len().min(window_len));
            let cut = ol_len.min(window_len / window_cut_ratio);

            //          < ol_len >
            //          ALIGN CUT
            // --------|-----|-->
            //     <--|-----|---------
            //     CUT ALIGN
            let overhang_i = &cons[i][cons[i].len() - ol_len..][cut / 2..ol_len - cut / 2];
            let overhang_j = &cons[i + 1][0..ol_len][cut / 2..ol_len - cut / 2];

            //dbg!(&overhang_i, &overhang_j);

            //i is query
            //log::trace!("Aligning overhangs for blocks {} and {}", i, i + 1);
            let (i_end, j_end, _) =
                alignment::align_seq_to_ref_slice_local(overhang_i, overhang_j, &GAPS_LAX_INDEL);
            //log::debug!("Alignment for block {} complete, i_end {} j_end {}", i);
            log::trace!(
                "Alignment for blocks {}, {} and contig {} complete, i_end {} j_end {}",
                i,
                i + 1,
                contig_name,
                i_end,
                j_end
            );

            //breakpoints.push((query_pos, ref_pos));
            breakpoints.lock().unwrap().push((i, i_end, j_end, cut));
        });
        let mut breakpoints = breakpoints.into_inner().unwrap();
        breakpoints.sort_by_key(|x| x.0);
        let breakpoints = breakpoints
            .into_iter()
            .map(|x| (x.1, x.2, x.3))
            .collect::<Vec<_>>();

        let mut new_consensus = vec![];
        let mut new_cons_i = std::mem::take(&mut cons[0]);
        let num_bps = breakpoints.len();
        for (i, bp) in breakpoints.into_iter().enumerate() {
            let mut new_cons_j = std::mem::take(&mut cons[i + 1]);
            let ol_len = new_cons_i.len().min(new_cons_j.len().min(window_len));
            let hang;

            //I believe this happens when the overlap alignments are discordant or near the ends of contigs
            if bp.0 > ol_len {
                if i != num_bps - 1 {
                    log::trace!("Potential error in consensus joining at block {} for contig {}. Iblock len {}, Jblock len {}", (i+1) * bp_length - window_len, contig_name, new_cons_i.len(), new_cons_j.len());
                }
                hang = 0;
            } else {
                hang = ol_len - bp.0;
            }

            let cut = bp.2;
            let break_pos_i = new_cons_i.len() - hang + cut / 2;
            let break_pos_j = bp.1 + cut / 2;

            new_cons_i.truncate(break_pos_i);
            new_cons_j = new_cons_j.split_off(break_pos_j);
            new_consensus.push(new_cons_i);
            new_cons_i = new_cons_j;
        }
        new_consensus.push(new_cons_i);

        return new_consensus;
    }

    pub fn new(
        genome_length: usize,
        contig_name: String,
        genome_string_u8: Vec<u8>,
        args: &Cli,
    ) -> Self {
        PoaConsensusBuilder {
            seq: Vec::new(),
            qual: Vec::new(),
            breakpoints: Vec::new(),
            contig_name,
            genome_length,
            bp_len: 0,
            window_overlap_len: 0,
            genome: genome_string_u8,
            args: args.clone(),
        }
    }

    #[cfg(test)]
    fn new_test(genome_length: usize) -> Self {
        PoaConsensusBuilder {
            seq: Vec::new(),
            qual: Vec::new(),
            breakpoints: Vec::new(),
            contig_name: "test".to_string(),
            genome_length,
            bp_len: 0,
            window_overlap_len: 0,
            genome: vec![],
            args: Cli::default(),
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

    pub fn num_blocks(&self) -> usize {
        self.breakpoints.len()
    }

    /// Returns (num_blocks_with_seqs, total_seqs, max_seqs_in_block, median_seqs_in_block, suspicious_blocks)
    /// suspicious_blocks = blocks where max seq length > 1000
    pub fn get_block_stats(&self) -> (usize, usize, usize, usize, usize, usize) {
        let mut counts: Vec<usize> = Vec::new();
        let mut suspicious_blocks = 0usize;
        let mut max_length = 0usize;

        for block in self.seq.iter() {
            let seqs = block.lock().unwrap();
            counts.push(seqs.len());
            let max_len = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
            if max_len > 1000 {
                suspicious_blocks += 1;
            }
            if max_len > max_length {
                max_length = max_len;
            }
        }

        let num_blocks_with_seqs = counts.iter().filter(|&&c| c > 0).count();
        let total_seqs: usize = counts.iter().sum();
        let max_seqs = *counts.iter().max().unwrap_or(&0);
        counts.sort();
        let median_seqs = if counts.is_empty() {
            0
        } else {
            counts[counts.len() / 2]
        };
        (
            num_blocks_with_seqs,
            total_seqs,
            max_seqs,
            median_seqs,
            suspicious_blocks,
            max_length,
        )
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
        {
            let mut seq_lock = self.seq[bp_i].lock().unwrap();
            let mut qual_lock = self.qual[bp_i].lock().unwrap();
            seq_lock.push(std::mem::take(current_block_string));
            qual_lock.push(std::mem::take(current_block_qual));
        }
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
                                bp_i - 1,
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
                current_position_q,
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

            let query_seq = &twin_reads[ol.query_id as usize].dna_seq;
            let query_seq_u8: Vec<u8>;

            if ol.reverse {
                query_seq_u8 = query_seq
                    .to_revcomp()
                    .iter()
                    .map(|x| x.to_char().to_ascii_uppercase() as u8)
                    .collect();
            } else {
                query_seq_u8 = query_seq
                    .iter()
                    .map(|x| x.to_char().to_ascii_uppercase() as u8)
                    .collect();
            }

            let query_quals_opt = &twin_reads[ol.query_id as usize].qual_seq;
            //query_quals = &twin_reads[ol.query_id as usize].qual_seq.as_ref().unwrap();
            let mut query_quals_u8: Vec<u8>;
            if let Some(query_quals) = query_quals_opt {
                let query_quals_u8_binned: Vec<u8>;
                if ol.reverse {
                    query_quals_u8_binned = query_quals
                        .to_revcomp()
                        .iter()
                        .map(|x| (x as u8) * 3 + 33)
                        .collect();
                } else {
                    query_quals_u8_binned =
                        query_quals.iter().map(|x| (x as u8) * 3 + 33).collect();
                }

                let bin_size = QUALITY_SEQ_BIN;
                query_quals_u8 = query_quals_u8_binned
                    .iter()
                    .flat_map(|x| vec![*x; bin_size])
                    .collect::<Vec<u8>>();

                if query_quals_u8.len() > query_seq.len() {
                    query_quals_u8.truncate(query_seq.len());
                } else if query_quals_u8.len() < query_seq.len() {
                    let last_qual = query_quals_u8[query_quals_u8.len() - 1];
                    query_quals_u8.extend(vec![last_qual; query_seq.len() - query_quals_u8.len()]);
                }
            } else {
                query_quals_u8 = vec![70; query_seq_u8.len()];
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

pub fn join_circular_ends(
    seq: &mut Vec<u8>,
    overlap_len: usize,
    hang1: usize,
    hang2: usize,
    contig_name: &str,
    args: &Cli,
) {
    let seq_len = seq.len();
    if seq_len < 1000 {
        return;
    }

    //let overlap_length = edge.overlap.overlap_len_bases + edge.overlap.hang1 + edge.overlap.hang2 + 1000;
    let mut overlap_length = overlap_len + hang1 + hang2 + (seq.len() / 10).min(500);

    if overlap_length > seq_len {
        //log::debug!("Overlap length is greater than sequence length for contig {}; something went wrong during polishing -- unsuccessfully circularized.", contig_name);
        overlap_length = seq_len * 3 / 4;
    }

    let overhang1 = seq[seq_len - overlap_length..seq_len].to_vec();
    let overhang2 = seq[0..overlap_length].to_vec();

    let tr1 = seeding::get_twin_read_syncmer(
        overhang1,
        None,
        args.kmer_size,
        args.c,
        &FxHashSet::default(),
        String::new(),
    )
    .unwrap();
    let tr2 = seeding::get_twin_read_syncmer(
        overhang2,
        None,
        args.kmer_size,
        args.c,
        &FxHashSet::default(),
        String::new(),
    )
    .unwrap();

    let mut lax_args = args.clone();
    lax_args.min_ol = args.min_ol / 2;

    let tr_options = CompareTwinReadOptions {
        compare_snpmers: false,
        retain_chain: false,
        //force_one_to_one_alignments: true,
        force_query_nonoverlap: true,
        ..Default::default()
    };

    let overlaps =
        mapping::compare_twin_reads(&tr1, &tr2, None, None, None, 0, 1, &tr_options, args);

    if overlaps.is_empty() {
        log::debug!(
            "Circular contig with hash id {} was not able to be end-polished correctly",
            &contig_name
        );
        return;
    }

    let best_overlap = overlaps.iter().max_by_key(|x| x.end1 - x.start1).unwrap();

    //dbg!(&best_overlap);
    let end_trim = overlap_length - best_overlap.start1;
    let start_trim = best_overlap.start2;
    log::trace!(
        "Overlap length {}, end_trim {}, start_trim {}, contig {}",
        overlap_length,
        end_trim,
        start_trim,
        contig_name
    );
    log::trace!(
        "Best overlap for circular contig {}: {:?}. New length: {}",
        contig_name,
        best_overlap,
        seq_len - end_trim - start_trim
    );
    let new_seq = seq[start_trim..seq_len - end_trim].to_vec();
    *seq = new_seq;
}

/// Trims consensus sequence ends based on alignment coverage.
/// Finds regions where fewer than half of the sequences cover the consensus,
/// and trims those regions from the ends.
pub fn trim_consensus_by_coverage(consensus: &[u8], seqs: &[Vec<u8>]) -> Vec<u8> {
    if consensus.is_empty() || seqs.is_empty() {
        return consensus.to_vec();
    }

    // Align each sequence to the consensus and collect coverage intervals
    let mut intervals: Vec<(usize, usize)> = Vec::with_capacity(seqs.len());
    let aligner = alignment::SeedChainAligner::new(consensus);

    for seq in seqs {
        // Skip sequences that are too short
        if seq.len() < 2 {
            continue;
        }

        // Remove null terminator if present
        let seq_clean = if seq.last() == Some(&0) {
            &seq[..seq.len() - 1]
        } else {
            &seq[..]
        };

        if seq_clean.is_empty() {
            continue;
        }

        // Align sequence to consensus (seed-chain-extend; falls back to full DP if needed)
        let (_query_end, ref_end, cigar) = aligner.align(seq_clean, &GAPS);

        // If alignment failed or is too poor, skip this sequence
        if cigar.is_empty() {
            continue;
        }

        // Parse CIGAR to get the reference (consensus) interval covered
        let (ref_start, ref_stop) = parse_cigar_ref_interval(&cigar, ref_end);

        if ref_start < ref_stop && ref_stop <= consensus.len() {
            intervals.push((ref_start, ref_stop));
        }
    }

    if intervals.is_empty() {
        return consensus.to_vec();
    }

    // Find trim positions using median of starts and ends
    let (trim_start, trim_end) = find_trim_positions(&intervals, consensus.len());

    log::trace!(
        "Trimming consensus from {} bases to {} bases (trimming {}-{} and {}-{})",
        consensus.len(),
        trim_end - trim_start,
        0,
        trim_start,
        trim_end,
        consensus.len()
    );

    // Return trimmed consensus
    if trim_start < trim_end && trim_end <= consensus.len() {
        consensus[trim_start..trim_end].to_vec()
    } else {
        consensus.to_vec()
    }
}

/// Parses CIGAR to find the reference interval that was aligned.
/// Returns (ref_start, ref_end) relative to the reference sequence.
fn parse_cigar_ref_interval(cigar: &[OpLen], ref_end: usize) -> (usize, usize) {
    // Work backwards from ref_end to find ref_start
    let mut ref_consumed = 0;
    let mut _query_consumed = 0;

    for op_len in cigar {
        match op_len.op {
            Operation::M | Operation::X | Operation::Eq => {
                ref_consumed += op_len.len;
                _query_consumed += op_len.len;
            }
            Operation::D => {
                ref_consumed += op_len.len;
            }
            Operation::I => {
                _query_consumed += op_len.len;
            }
            _ => {}
        }
    }

    // ref_end is the position where alignment ends
    // ref_start is ref_end - ref_consumed
    let ref_start = ref_end.saturating_sub(ref_consumed);

    (ref_start, ref_end)
}

/// Finds trim positions by taking the median of interval starts and ends.
/// Trims ends where fewer than half of sequences have coverage.
fn find_trim_positions(intervals: &[(usize, usize)], consensus_len: usize) -> (usize, usize) {
    if intervals.is_empty() {
        return (0, consensus_len);
    }

    let half_count = intervals.len() / 2;

    // Sort starts and ends separately
    let mut starts: Vec<usize> = intervals.iter().map(|(s, _)| *s).collect();
    let mut ends: Vec<usize> = intervals.iter().map(|(_, e)| *e).collect();

    starts.sort_unstable();
    ends.sort_unstable();

    // Take the position where at least half of sequences start covering
    // This is approximately the median start position
    let trim_start = if starts.len() > half_count {
        starts[half_count]
    } else {
        0
    };

    // Take the position where at least half of sequences stop covering
    // This is approximately the median end position
    let trim_end = if ends.len() > half_count {
        ends[half_count]
    } else {
        consensus_len
    };

    // Sanity check
    if trim_start >= trim_end {
        return (0, consensus_len);
    }

    (trim_start, trim_end)
}

/// abpoa-based POA consensus using the same call convention as `poa_consensus`.
///
/// `seqs` and `quals` are null-terminated (stripped internally before passing to abpoa).
/// `mismatch_score`, `gap_open`, and `gap_extend` follow the spoa sign convention (negative).
/// A new `Aligner` is created per call; this is intentional because `Aligner` is `!Send+!Sync`
/// and this function is called inside a rayon parallel closure.
fn abpoa_consensus_impl(
    seqs: &[Vec<u8>],
    quals: &[Vec<u8>],
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> Vec<u8> {
    // Strip null terminators — abpoa takes plain byte slices
    let clean: Vec<&[u8]> = seqs
        .iter()
        .map(|s| {
            if s.last() == Some(&0) {
                &s[..s.len() - 1]
            } else {
                &s[..]
            }
        })
        .collect();

    // Convert Phred+33 quality bytes to non-negative i32 weights
    let weights: Vec<Vec<i32>> = quals
        .iter()
        .zip(clean.iter())
        .map(|(q, s)| {
            let q_data = if q.last() == Some(&0) {
                &q[..q.len() - 1]
            } else {
                &q[..]
            };
            (0..s.len())
                .map(|i| {
                    if i < q_data.len() {
                        (q_data[i] as i32 - 33).max(0)
                    } else {
                        0
                    }
                })
                .collect()
        })
        .collect();
    let weight_refs: Vec<&[i32]> = weights.iter().map(|w| w.as_slice()).collect();

    // abpoa uses positive penalty magnitudes; negate spoa-convention values
    let mismatch_pen = (-mismatch_score).max(1);
    let gap_open_pen = (-gap_open).max(1);
    let gap_ext_pen = (-gap_extend).max(1);

    let mut params = match Parameters::configure() {
        Ok(p) => p,
        Err(_) => return vec![],
    };
    params
        .set_outputs(OutputMode::CONSENSUS)
        .set_align_mode(AlignMode::Global)
        .set_use_quality(true)
        .set_use_read_ids(false)
        .set_consensus(abpoa::ConsensusAlgorithm::HeaviestBundle, 1, 0.00)
        .unwrap();

    if params
        .set_scoring_scheme(Scoring::affine(
            match_score,
            mismatch_pen,
            gap_open_pen,
            gap_ext_pen,
        ))
        .is_err()
    {
        return vec![];
    }

    let mut aligner = match Aligner::with_params(params) {
        Ok(a) => a,
        Err(_) => return vec![],
    };

    let batch = match SequenceBatch::from_sequences(&clean)
        .and_then(|b| b.with_quality_weights(&weight_refs))
    {
        Ok(b) => b,
        Err(_) => return vec![],
    };

    match aligner.msa(batch) {
        Ok(result) => {
            if result.clusters.is_empty() {
                vec![]
            } else {
                result.clusters[0].consensus.as_bytes().to_vec()
            }
        }
        _ => {
            vec![]
        }
    }
}

/// All homopolymer-compressed reads for one polishing window.
/// `seqs` and `quals` are null-terminated and can be passed directly to `poa_consensus`
/// and `trim_consensus_by_coverage` without any extra allocation.
struct HpcSeqsCompressed {
    seqs: Vec<Vec<u8>>,         // null-terminated HPC sequences
    quals: Vec<Vec<u8>>,        // null-terminated median qualities (one per HPC base)
    run_lengths: Vec<Vec<u32>>, // [read][hpc_base] -> original run length
}

/// Compress all `seqs`/`quals` (both null-terminated) into an [`HpcSeqsCompressed`].
fn hpc_compress_all(seqs: &[Vec<u8>], quals: &[Vec<u8>]) -> HpcSeqsCompressed {
    let mut hpc_seqs = Vec::with_capacity(seqs.len());
    let mut hpc_quals = Vec::with_capacity(quals.len());
    let mut hpc_runs = Vec::with_capacity(seqs.len());

    for (seq, qual) in seqs.iter().zip(quals.iter()) {
        let seq_data = if seq.last() == Some(&0) {
            &seq[..seq.len() - 1]
        } else {
            seq
        };
        let qual_data = if qual.last() == Some(&0) {
            &qual[..qual.len() - 1]
        } else {
            qual
        };
        let qual_data = &qual_data[..qual_data.len().min(seq_data.len())];

        let mut hpc_seq = Vec::new();
        let mut hpc_qual = Vec::new();
        let mut runs = Vec::new();

        let mut i = 0;
        while i < seq_data.len() {
            let base = seq_data[i];
            let run_start = i;
            while i < seq_data.len() && seq_data[i] == base {
                i += 1;
            }
            let run_len = i - run_start;
            let q_end = (run_start + run_len).min(qual_data.len());
            let run_quals = &qual_data[run_start.min(qual_data.len())..q_end];
            let median_qual = if run_quals.is_empty() {
                48u8
            } else {
                let mut sorted = run_quals.to_vec();
                sorted.sort_unstable();
                sorted[sorted.len() / 2]
            };
            hpc_seq.push(base);
            hpc_qual.push(median_qual);
            runs.push(run_len as u32);
        }

        hpc_seq.push(0u8);
        hpc_qual.push(0u8);

        hpc_seqs.push(hpc_seq);
        hpc_quals.push(hpc_qual);
        hpc_runs.push(runs);
    }

    HpcSeqsCompressed {
        seqs: hpc_seqs,
        quals: hpc_quals,
        run_lengths: hpc_runs,
    }
}

/// Expand an HPC consensus back to full length by voting on run lengths from aligned reads.
///
/// For each position in `hpc_cons`, each read that aligns there votes for its original run
/// length weighted by its median quality.  The run length with the highest total weight wins.
///
/// Note: this performs one alignment per read.  In the future this can be unified with
/// `trim_consensus_by_coverage` to share alignments.
fn expand_hpc_consensus(hpc_cons: &[u8], hpc: &HpcSeqsCompressed) -> Vec<u8> {
    if hpc_cons.is_empty() {
        return vec![];
    }

    // votes[r] = (run_length, quality_weight) from reads that matched consensus position r
    let mut votes: Vec<Vec<(u32, u8)>> = vec![vec![]; hpc_cons.len()];

    let aligner = alignment::SeedChainAligner::new(hpc_cons);
    for (read_idx, seq) in hpc.seqs.iter().enumerate() {
        let seq_data = if seq.last() == Some(&0) {
            &seq[..seq.len() - 1]
        } else {
            &seq[..]
        };
        if seq_data.is_empty() {
            continue;
        }

        let (query_end, ref_end, cigar) = aligner.align(seq_data, &GAPS);
        if cigar.is_empty() {
            continue;
        }

        // Derive alignment start positions from CIGAR totals
        let mut ref_consumed = 0usize;
        let mut q_consumed = 0usize;
        for op in &cigar {
            match op.op {
                Operation::M | Operation::X | Operation::Eq => {
                    ref_consumed += op.len;
                    q_consumed += op.len;
                }
                Operation::D => {
                    ref_consumed += op.len;
                }
                Operation::I => {
                    q_consumed += op.len;
                }
                _ => {}
            }
        }
        let mut r_pos = ref_end.saturating_sub(ref_consumed);
        let mut q_pos = query_end.saturating_sub(q_consumed);

        let runs = &hpc.run_lengths[read_idx];
        let quals = &hpc.quals[read_idx]; // null at [runs.len()]; safe to index up to runs.len()-1

        for op in &cigar {
            match op.op {
                Operation::M | Operation::X | Operation::Eq => {
                    for _ in 0..op.len {
                        if q_pos < runs.len() && r_pos < votes.len() {
                            votes[r_pos].push((runs[q_pos], quals[q_pos]));
                        }
                        q_pos += 1;
                        r_pos += 1;
                    }
                }
                Operation::I => {
                    q_pos += op.len;
                }
                Operation::D => {
                    r_pos += op.len;
                }
                _ => {}
            }
        }
    }

    // For each HPC base, repeat it the weighted-mode number of times
    let mut expanded = Vec::new();
    for (r, &base) in hpc_cons.iter().enumerate() {
        let run_len = if votes[r].is_empty() {
            1u32
        } else {
            let mut scores: FxHashMap<u32, u32> = FxHashMap::default();
            for &(rl, w) in &votes[r] {
                *scores.entry(rl).or_insert(0) += w as u32;
            }
            scores
                .into_iter()
                .max_by_key(|(_, v)| *v)
                .map(|(k, _)| k)
                .unwrap_or(1)
        };
        for _ in 0..run_len {
            expanded.push(base);
        }
    }

    expanded
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
            let qual = vec![50; 32];
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
            let qual = vec![50; 32];
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
                if cons.len() != 0 {
                    println!("{:?}", String::from_utf8_lossy(&cons));
                }
            }

            if window_length == 0 {
                //insertion at border
                assert!(cons[0].len() == 5);
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

        let seqs: Vec<Vec<u8>> = vec![
        b"CTGATAACTCAAACTGCTGTCAGGCATTGCCAGAACAGCAAGATAATGAATGCCAAATTCATCTATCTTATTGAATACGATGTCACTCAATATAGACTGGTCCAGGAAAGAAGAAGGAACCTGCTTATCGAAATCCTTCTTATGAAAATCACGAGAGAACAGACAGTTGGCACCATGATAAACATGCAAGTTCACAATATTATCATAATAAACATTGTCTACCTCAACACCATCATCATTATAAGAACTCTTGTAGACCTTATAACTTGTTGGGTTAACCTGTACATAGAGATGATATTTCCCATCGCCCCTTACAACAACCGTATCACGTTTAATCAGCGTATTCTGATTCAGCGCCACGGAAGCATGGTTGTGGATGAACTGCTGCAAATAAGACTTATCTTCCGTCTTCACCAGCTTAACCACATCCCCATTCTGCACTTTAAACTGGAAAAGTACAATGCTCAGCCTGTTTTACGATCGGATATTTTGCAATATTAGC".to_vec(),
        b"CTGATAACTCAAACTGCTGTCAGGCATTGCCAGAACAGCAAGATAATGAATGCCAAATTCATCTATCTTATTGAATACGATGTCACTCAATATAGACTGGTCCAGGAAAGAAGAAGGAACCTGCTTATCGAAATCCTTCTTATGAAAATCACGAGAAGAAACAAGCAGTTGGCACCATGATAAACATGCAAGTTCACAATATTATCATAATAAACATTGTCTACCTCAACACCATCATCATTATAAGAACTCTTGTAGACCTTATAACTTGTTGGGTTAACCTGTACATAGAGATGATATTTCCCATCGCCCCTTACAACAACCGTATCACGTTTAATCAGCGTATTCTGATTCAGCGCCACGGAAGCATGGTTGTGGATGAACTGCTGCAAATAAGACTTATCTTCCGTCTTCACCAGCTTAACCACATCCCCATTCTGCACTTTAAACTGGAAAATATGCTCAGCCTGTTTTACGATCGGATATTTTGCAATATTAGC".to_vec(),
        b"CTGATAACTCAAACTGCTGTCAGGCATTGCCAGAACAGCAAGATAATGAATGCCAAATTCATCTCTGTTGCTTGAATACGATGTCACTCAATATAGACTGGTCCAGAAAGAAGAAGGAACCTGCTTATCGAAATCCTTCTTATGAAAATCACAAGAAAGAGAACAGACAGTTGGCACCATGATAAACATGCAAGTTCACAATATTATCATAATAAACATTGTCTACCTCAACACCATCATCATTATAAGAACTCTTGTAGACCTTATAACTTGTTGGGTTAACCTGTACATAGAGATGATATTTCCCATCGCCCCTTACAACAACCGTATCACGTTTAATCAGCGTATTCTGATTCAGCGCCACGGAAGCATGGTTGTGGATGAACTGCTGCAAATAAGACTTATCTTCCGTCTTCACCAGCTTAACCACATCCCCATTCTGCACTTTAAACTGAAAAATATTCCTCAGCCTGTTTTACGATCGGATATTTTGCAATATTAGC".to_vec(),
        b"CTGATAACTCAAACTGCTGTCAGGCATTGCCAGAACAGCAAGATAATGAATGCCAAATTCATCTATCTTATTGAATACGATGTCACTCAATATAGACTGGTCCAGGAAAGAAGAAGGAACCTGCTTATCGAAATCCTTCTTATGAAAATCACGAGAGAACAGACAGTTGGCACCATGATAAACATGCAAGTTCACAATATTATCATAATAAACATTGTCTACCTCAACACCATCATCATTATAAGAACTCTTGTAGACCTTATAACTTGTTGGGTTAACCTGTACATAGAGATGATATTTCCCATCGCCCCTTACAACAACCGTATCACGTTTAATCAGCGTATTCTGATTCAGCGCCACGGAAGCATGGTTGTGGATGAACTGCTGCAAATAAGACTTATCTTCCGTCTTCACCAGCTTAACCACATCCCCATTCTGCACTTTAAACTGGAAAATATGCTCAGCCTGTTTTACGATCGGATATTTTGCAATATTAGC".to_vec(),
        b"CATTGCCAGAACAGCAAGATAATGAATGCCAAATTCATCTATCTTATTGAATACGATGTCACTCAATATAGACTGGTCCAGGAAAGAAGAAGGAACCTGCTTATCGAAATCCTTCTTATGAAAATCACGAGAGAACAGACAGTTGGCACCATGATAAACATGCAAGTTCACAATATTATCATAATAAACATTGTCTACCTCAACACCATCATCATTATCAACTCTTGTAGACCTTATAACTTGTTGGGTTAACCTGTACATAGAGATGATATTTCCCATCGCCCCTTACAACAACCGTATCACGTTTAATCAGCGTATTCTGATTTCAGCGCCACAAAGTCCTGGTTGTGGATGAACTGCTGCAAATAAGACTTATCTTCCGTCTTCACCAGCTTAACCACATCCCCCATTCTGCACTTACAACTGGAAAATATGCTCAGCCTGTTTTACGATCGGATATTTTGCAATATTAGC".to_vec(),
        b"CTGATAACTCTAAATCTTGCTGTCAGGCAT".to_vec(),
        ];
        let mut quals = seqs.iter().map(|x| vec![50; x.len()]).collect::<Vec<_>>();

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
            builder.add_seq(
                seqs[i].clone(),
                quals[i].clone(),
                &OpLenVec::new(cigars[i].clone()),
                0,
                0,
            );
        }

        let consensuses = builder.spoa_blocks();
        println!("Consensuses: {:?}", consensuses[0]);
        assert!(consensuses[0].len() > 480 && consensuses[0].len() < 520);
    }

    #[test]
    fn test_modify_new() {
        let consensuses = vec![
            b"GCAACGTATGT".to_vec(),    // last C is error
            b"ACGTATGTGTGTGT".to_vec(), // first T is errors
        ];

        let new_cons = PoaConsensusBuilder::modify_join_consensus(consensuses, 8, 20, "test");

        let mut final_consensus = vec![];
        for cons in new_cons {
            println!("{:?}", String::from_utf8_lossy(&cons));
            final_consensus.extend(cons);
        }

        assert_eq!(final_consensus, b"GCAACGTATGTGTGTGT".to_vec());
    }

    #[test]
    fn test_modify_new_del() {
        let consensuses = vec![
            b"AAAAATTTA".to_vec(), // last A is error
            b"CTTTGGGG".to_vec(),  // first T is errors
        ];

        let new_cons = PoaConsensusBuilder::modify_join_consensus(consensuses, 5, 20, "test");

        let mut final_consensus = vec![];
        for cons in new_cons {
            println!("{:?}", String::from_utf8_lossy(&cons));
            final_consensus.extend(cons);
        }

        assert_eq!(final_consensus, b"AAAAATTTGGGG".to_vec());
    }

    #[test]
    fn circular_join_basic_test() {
        let mut args = Cli::default();
        args.c = 5;
        args.kmer_size = 17;
        {
            let random_string = b"GCATGCGTTCAACGTAGGCCGTACTAGCTGCGTAATCGACGGAATGGCAGTATCGCGATAACGCTTGAAACGCTACGAGCCATAGCGGTATCGTAGCAACGCTAATCGGCATAGCTATCGATGCAGTCGCTATAGCTAGCTAGCGATCGGCCGATAGCGATCGATCGGCTAGCGGCATCGATAGCGGCCGATCGCGATCAGCATGGCCGATGCGATCGCGTATCAGCGCGATCGAGCCGATCGATCGCGTCCGATGCATGCAACGATCGGCATATCACGCGCGATCGACTAGCGATCGATCGCGTACGCATCGATCGAGCGATCGACTGATCGCTAGCTGCATGCATACGCTAGCTGCAGCTAGCATCGATCGCTATGCTAGCTAGCATCGAGCTGATCGTAGCATCGATCGATCGATCGATCGATCGAGCTATCGATCGATACGCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGCGATCGACTGCGATCGCTAGCTAGCTAGCTATGCTAGCTAGCTGCTAGTCGACGATCGATCGATCGATCGATCTAGCTAGCATCGCTAGCTGATCGTAGCTAGCTAGCATCGATCGA".to_vec();
            let mut seq = random_string.clone();
            seq.extend(vec![b'G'; 250]);
            seq.extend(vec![b'A'; 250]);
            seq.extend(vec![b'T'; 250]);
            seq.extend(vec![b'C'; 250]);
            seq.extend(random_string.clone());

            dbg!(seq.len());

            join_circular_ends(&mut seq, random_string.len(), 11, 11, "test", &args);

            assert_eq!(seq.len(), 1000 + random_string.len());
        }
    }

    #[test]
    fn circular_join_basic_test_fuzzy() {
        let mut args = Cli::default();
        args.c = 5;
        args.kmer_size = 17;
        {
            let random_string = b"GCATGCGTTCAACGTAGGCCGTACTAGCTGCGTAATCGACGGAATGGCAGTATCGCGATAACGCTTGAAACGCTACGAGCCATAGCGGTATCGTAGCAACGCTAATCGGCATAGCTATCGATGCAGTCGCTATAGCTAGCTAGCGATCGGCCGATAGCGATCGATCGGCTAGCGGCATCGATAGCGGCCGATCGCGATCAGCATGGCCGATGCGATCGCGTATCAGCGCGATCGAGCCGATCGATCGCGTCCGATGCATGCAACGATCGGCATATCACGCGCGATCGACTAGCGATCGATCGCGTACGCATCGATCGAGCGATCGACTGATCGCTAGCTGCATGCATACGCTAGCTGCAGCTAGCATCGATCGCTATGCTAGCTAGCATCGAGCTGATCGTAGCATCGATCGATCGATCGATCGATCGAGCTATCGATCGATACGCGATCGATCGATCGCGATCGATCGATCGATCGCGATCGATCGCGATCGACTGCGATCGCTAGCTAGCTAGCTATGCTAGCTAGCTGCTAGTCGACGATCGATCGATCGATCGATCTAGCTAGCATCGCTAGCTGATCGTAGCTAGCTAGCATCGATCGA".to_vec();
            let mut seq = b"ACGTA".to_vec();
            seq.extend(&random_string.clone());
            seq.extend(vec![b'G'; 250]);
            seq.extend(vec![b'A'; 250]);
            seq.extend(vec![b'T'; 250]);
            seq.extend(vec![b'C'; 250]);
            seq.extend(random_string.clone());
            seq.extend(b"GTGTGG");

            dbg!(seq.len());

            join_circular_ends(&mut seq, random_string.len(), 11, 11, "test", &args);

            assert_eq!(seq.len(), 1000 + random_string.len());
        }
    }

    // This test fails, hence we should use HPC instead....
    //#[test]
    fn _test_poa_homopolymer_distribution() {
        // 25-base high-entropy flanks surrounding a poly-A homopolymer
        let left_flank = b"ACGTACGATCGCATGCTACGATCGA"; // 25 bases
        let right_flank = b"TGCATCGATCGATCGTAGCATGCAT"; // 25 bases

        // Triangular distribution of homopolymer lengths centered at 10.
        // Range 5-15; lengths 5 and 15 have 0 sequences so are omitted.
        //   len:   6   7   8   9  10  11  12  13  14
        //   count: 2   4   6   8  10   8   6   4   2
        let distribution: &[(usize, usize)] = &[
            (6, 2),
            (7, 4),
            (8, 6),
            (9, 8),
            (10, 10),
            (11, 8),
            (12, 6),
            (13, 4),
            (14, 2),
        ];

        let mut seqs: Vec<Vec<u8>> = vec![];
        let mut quals: Vec<Vec<u8>> = vec![];

        for &(hp_len, count) in distribution {
            let mut seq = Vec::new();
            seq.extend_from_slice(left_flank);
            seq.extend(std::iter::repeat(b'A').take(hp_len));
            seq.extend_from_slice(right_flank);
            seq.push(0u8); // null terminator expected by poa_consensus

            let mut qual = vec![50u8; seq.len() - 1]; // same length as seq excluding null
            qual.push(0u8); // null terminator required by poa_consensus

            for _ in 0..count {
                seqs.push(seq.clone());
                quals.push(qual.clone());
            }
        }

        // Use same scoring as production spoa_blocks
        let consensus = poa_consensus(&seqs, &quals, 100, 1, 5, -2, -2, -1);

        // Find the longest poly-A run in the consensus
        let mut max_run = 0usize;
        let mut cur_run = 0usize;
        for &b in &consensus {
            if b == b'A' {
                cur_run += 1;
                max_run = max_run.max(cur_run);
            } else {
                cur_run = 0;
            }
        }

        assert!(
            max_run >= 9 && max_run <= 11,
            "Expected poly-A homopolymer length ~10, got {}. Consensus: {}",
            max_run,
            std::str::from_utf8(&consensus).unwrap_or("(invalid utf8)")
        );
    }

    #[test]
    fn test_poa_cons_stall() {
        let mut builder = PoaConsensusBuilder::new_test(10000);

        builder.generate_breakpoints(700, 50);

        let seqs: Vec<Vec<u8>> = vec![b"GAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA".to_vec(),
        b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT".to_vec(),];
        let quals = seqs.iter().map(|x| vec![33; x.len()]).collect::<Vec<_>>();
        //quals[0][5] = 50;

        let cigars = seqs
            .iter()
            .map(|x| {
                vec![OpLen {
                    op: Operation::M,
                    len: x.len(),
                }]
            })
            .collect::<Vec<_>>();

        for i in 0..seqs.len() {
            builder.add_seq(
                seqs[i].clone(),
                quals[i].clone(),
                &OpLenVec::new(cigars[i].clone()),
                0,
                0,
            );
        }

        let consensuses = builder.spoa_blocks();
        println!("Consensuses: {:?}", consensuses[0]);
    }

    #[test]
    fn test_hpc_compress_all() {
        // Two reads: one with a 3-base poly-A, one with a 2-base poly-A
        let seqs = vec![b"GCAAATGC\0".to_vec(), b"GCAATGC\0".to_vec()];
        // Qualities: 50 for every real base, null at end
        let quals: Vec<Vec<u8>> = seqs
            .iter()
            .map(|s| {
                let mut q = vec![50u8; s.len() - 1];
                q.push(0u8);
                q
            })
            .collect();

        let hpc = hpc_compress_all(&seqs, &quals);

        // Both reads compress to GCATGC\0 — no more adjacent identical bases
        assert_eq!(&hpc.seqs[0], b"GCATGC\0");
        assert_eq!(&hpc.seqs[1], b"GCATGC\0");

        // Run lengths for read 0: G=1, C=1, A=3, T=1, G=1, C=1
        assert_eq!(hpc.run_lengths[0], vec![1, 1, 3, 1, 1, 1]);
        // Run lengths for read 1: G=1, C=1, A=2, T=1, G=1, C=1
        assert_eq!(hpc.run_lengths[1], vec![1, 1, 2, 1, 1, 1]);

        // Median quality of the poly-A run in read 0 (all 50) should be 50
        assert_eq!(hpc.quals[0][2], 50);
    }

    #[test]
    fn test_hpc_expand_recovers_majority_run_length() {
        // Reads: 3 with poly-A length 5, 2 with poly-A length 3.
        // After compress → POA → expand, the poly-A should be length 5.
        let make_seq = |poly_len: usize| {
            let mut s = b"GCATGC".to_vec();
            s.extend(vec![b'A'; poly_len]);
            s.extend_from_slice(b"TGCATG");
            s.push(0u8);
            s
        };

        let seqs: Vec<Vec<u8>> = [1, 2, 3, 4, 5, 6, 7, 3, 4, 2, 3]
            .iter()
            .map(|&n| make_seq(n))
            .collect();
        let quals: Vec<Vec<u8>> = seqs
            .iter()
            .map(|s| {
                let mut q = vec![50u8; s.len() - 1];
                q.push(0u8);
                q
            })
            .collect();

        let consensus_without_hpc = poa_consensus(&seqs, &quals, 50, 1, 5, -2, -2, -1);
        // There will be 5 As in the consensus...
        assert_eq!(consensus_without_hpc, b"GCATGCAAAAAAATGCATG".to_vec());

        let hpc = hpc_compress_all(&seqs, &quals);
        // All HPC seqs are identical: GCATGCATGCATG\0
        let hpc_cons = poa_consensus(&hpc.seqs, &hpc.quals, 50, 1, 5, -2, -2, -1);
        let expanded = expand_hpc_consensus(&hpc_cons, &hpc);

        // Longest poly-A run in the expanded sequence should be 3
        let mut max_run = 0usize;
        let mut cur = 0usize;
        for &b in &expanded {
            if b == b'A' {
                cur += 1;
                max_run = max_run.max(cur);
            } else {
                cur = 0;
            }
        }
        assert_eq!(
            max_run,
            3,
            "Expected poly-A of length 3, got {}. Expanded: {}",
            max_run,
            std::str::from_utf8(&expanded).unwrap_or("?")
        );
    }

    // Helper: build null-terminated seqs/quals from plain string slices
    fn make_null_terminated(seqs_str: &[&str], qual_byte: u8) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
        let seqs: Vec<Vec<u8>> = seqs_str
            .iter()
            .map(|s| {
                let mut v = s.as_bytes().to_vec();
                v.push(0u8);
                v
            })
            .collect();
        let quals: Vec<Vec<u8>> = seqs
            .iter()
            .map(|s| {
                let mut q = vec![qual_byte; s.len() - 1];
                q.push(0u8);
                q
            })
            .collect();
        (seqs, quals)
    }

    #[test]
    fn test_abpoa_basic() {
        // Same 6-sequence set used in test_dna_consensus; expected consensus is "AATGCCCGTT"
        let (seqs, quals) = make_null_terminated(
            &[
                "ATTGCCCGTT",
                "AATGCCGTT",
                "AATGCCCGAT",
                "AACGCCCGTC",
                "AGTGCTCGTT",
                "AATGCTCGTT",
            ],
            50,
        );
        let cons = abpoa_consensus_impl(&seqs, &quals, 5, -2, -2, -1);
        assert_eq!(
            cons,
            b"AATGCCCGTT".to_vec(),
            "abpoa basic consensus mismatch: got {}",
            std::str::from_utf8(&cons).unwrap_or("?")
        );
    }

    #[test]
    fn test_abpoa_homopolymer_distribution() {
        // Triangular distribution of poly-A lengths centered at 10 (same as test_poa_homopolymer_distribution)
        let left_flank = "ACGTACGATCGCATGCTACGATCGA";
        let right_flank = "TGCATCGATCGATCGTAGCATGCAT";
        let distribution: &[(usize, usize)] = &[
            (6, 2),
            (7, 4),
            (8, 6),
            (9, 8),
            (10, 10),
            (11, 8),
            (12, 6),
            (13, 4),
            (14, 2),
        ];

        let mut seqs: Vec<Vec<u8>> = vec![];
        let mut quals: Vec<Vec<u8>> = vec![];
        for &(hp_len, count) in distribution {
            let seq_str = format!("{}{}{}", left_flank, "A".repeat(hp_len), right_flank);
            for _ in 0..count {
                let mut s = seq_str.as_bytes().to_vec();
                s.push(0u8);
                let mut q = vec![50u8; s.len() - 1];
                q.push(0u8);
                seqs.push(s);
                quals.push(q);
            }
        }

        let cons = abpoa_consensus_impl(&seqs, &quals, 5, -2, -2, -1);

        let mut max_run = 0usize;
        let mut cur = 0usize;
        for &b in &cons {
            if b == b'A' {
                cur += 1;
                max_run = max_run.max(cur);
            } else {
                cur = 0;
            }
        }
        assert!(
            max_run >= 9 && max_run <= 11,
            "Expected poly-A homopolymer length ~10, got {}. Consensus: {}",
            max_run,
            std::str::from_utf8(&cons).unwrap_or("?")
        );
    }
}
