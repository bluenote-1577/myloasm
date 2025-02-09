use rust_lapper::{Interval, Lapper};
use smallvec::SmallVec;
use fxhash::FxHashMap;
use block_aligner::cigar::*;
use rust_spoa::poa_consensus;

#[derive(Debug, Clone, Default)]
pub struct HomopolymerCompressedSeq{
    pub seq: Vec<u8>,
    pub qualities: Vec<u8>,
    pub homopolymer_lengths: Option<Vec<u16>>,
}

impl HomopolymerCompressedSeq{
    pub fn new_slice(seq: &[u8], qualities: &[u8], hpc:bool) -> Self {
        return HomopolymerCompressedSeq::new(seq, qualities.to_vec(), hpc);
    }
    pub fn new(seq: &[u8], mut qualities: Vec<u8>, hpc: bool) -> Self {
        let sequence;
        qualities.iter_mut().for_each(|x| *x -= 33);
        let hompolymer_lengths_opt;
        if hpc{
            let mut compressed_seq = Vec::new();
            let mut homopolymer_lengths = Vec::new();
            let mut current_base = seq[0];
            let mut current_length:u16 = 1;
            for i in 1..seq.len(){
                if seq[i] == current_base{
                    current_length = current_length.saturating_add(1);
                } else {
                    compressed_seq.push(current_base);
                    homopolymer_lengths.push(current_length);
                    current_base = seq[i];
                    current_length = 1;
                }
            }
            compressed_seq.push(current_base);
            homopolymer_lengths.push(current_length);
            sequence = compressed_seq;
            hompolymer_lengths_opt = Some(homopolymer_lengths);
        }
        else{
            sequence = seq.to_vec();
            hompolymer_lengths_opt = None;
        }

        //TODO this doesn't work for HPC right now. 
        assert!(sequence.len() == qualities.len());

        Self { seq: sequence, qualities: qualities, homopolymer_lengths: hompolymer_lengths_opt}
    }
}


#[derive(Debug, Default, Clone)]
pub struct InsertionCounts {
    // Maps inserted sequence to count of times it was seen
    sequences: FxHashMap<SmallVec<[u8;2]>, u32>,
    // AT -> {index 0 : {1: 10, 2 : 1, 3 : 1},index 1: {2:8, 1: 2}} -> ATT
    sequence_hpc_lengths: Option<FxHashMap<SmallVec<[u8;2]>, SmallVec<[FxHashMap<u16, u16>;2]>>>,
}

impl InsertionCounts {
    fn new(hpc:bool) -> Self {
        let sequence_hpc_lengths = if hpc { Some(FxHashMap::default()) } else { None };
        Self {
            sequences: FxHashMap::default(),
            sequence_hpc_lengths,
        }
    }

    fn add_insertion(&mut self, sequence: &[u8], qualities: &[u8], hpc_lengths: Option<&[u16]>) {
        let c = self.sequences.entry(SmallVec::from_slice(sequence)).or_insert(0);
        let mean_but_actually_geom_mean = qualities.iter().map(|x| *x as u32).sum::<u32>() / qualities.len() as u32;
        *c = c.saturating_add(mean_but_actually_geom_mean);

        if let Some(hpc_lengths) = hpc_lengths{
            let hpc_lengths_ins = self.sequence_hpc_lengths.as_mut().unwrap().entry(SmallVec::from_slice(sequence)).or_insert(SmallVec::new());
            if hpc_lengths_ins.len() == 0{
                for _ in 0..sequence.len(){
                    hpc_lengths_ins.push(FxHashMap::default());
                }
            }
            for i in 0..sequence.len(){
                let count_i = hpc_lengths[i];
                hpc_lengths_ins[i].entry(count_i).and_modify(|e| *e = e.saturating_add(1)).or_insert(1);
            }
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct BaseConsensus {
    pub a_count: u32,
    pub c_count: u32,
    pub g_count: u32,
    pub t_count: u32,
    pub homopolymer_lengths: Option<FxHashMap<u16, u16>>,
    pub deletion_count: u32,
    pub insertions: InsertionCounts,
    pub prev_ins_weight: u32,
    pub prev_nonins_weight: u32,
}

impl BaseConsensus {
    fn new(hpc: bool) -> Self {
        let homopolymer_lengths = if hpc { Some(FxHashMap::default()) } else { None };
        Self {
            a_count: 0,
            c_count: 0,
            g_count: 0,
            t_count: 0,
            deletion_count: 0,
            homopolymer_lengths,
            prev_ins_weight: 0,
            prev_nonins_weight: 0,
            insertions: InsertionCounts::new(hpc),
        }
    }

    fn increment_base(&mut self, base: u8, quality: u8, hpc_length: Option<u16>) {
        let quality = quality as u32;
        match base.to_ascii_uppercase() {
            b'A' => self.a_count = self.a_count.saturating_add(quality),
            b'C' => self.c_count = self.c_count.saturating_add(quality),
            b'G' => self.g_count = self.g_count.saturating_add(quality),
            b'T' => self.t_count = self.t_count.saturating_add(quality),
            _ => (), // Ignore other characters
        }
        if let Some(hpc_length) = hpc_length{
            let opt_hpc = self.homopolymer_lengths.as_mut();
            let count = opt_hpc.unwrap().entry(hpc_length).or_insert(0);
            *count = count.saturating_add(1);
        }
    }

    fn increment_deletion(&mut self, quality: u8) {
        self.deletion_count = self.deletion_count.saturating_add(quality as u32);
    }

    fn add_insertion(&mut self, sequence: &[u8], qualities: &[u8], hpc_lengths: Option<&[u16]>) {
        self.insertions.add_insertion(sequence, qualities, hpc_lengths);
    }
}

#[derive(Debug)]
pub struct ConsensusBuilder {
    pub consensus: Vec<BaseConsensus>,
    pub overlaps: Vec<(u32,u32, ReadCigarIndex)>,
    pub hpc: bool,
}

impl ConsensusBuilder {
    pub fn new(reference_length: usize, hpc: bool) -> Self {
        Self {
            //+1 because end insertion
            consensus: vec![BaseConsensus::new(hpc); reference_length + 1],
            overlaps: vec![],
            hpc
        }
    }

    pub fn add_overlap(&mut self, read_index: usize, reference_index: usize, reverse: bool, cigar_index: Vec<OpLen>, alignment_start: usize) {
        let read_cigar_index = ReadCigarIndex::new(read_index, reference_index, reverse, cigar_index, alignment_start);
        self.overlaps.push((reference_index as u32, reference_index as u32 + 1, read_cigar_index));
    }

    fn _weighted_av_hpc(hpc_lengths: &FxHashMap<u16, u16>) -> u16 {
        return hpc_lengths.iter().max_by_key(|(_, v)| *v).map(|(k, _)| *k).unwrap_or(0);
    }

    fn weighted_av_hpc(hpc_lengths: &FxHashMap<u16, u16>) -> u16 {
        let mut medians = vec![];
        let mut count = 0;
        for (k, v) in hpc_lengths{
            count += v;
            for _ in 0..*v{
                medians.push(*k);
            }
        }
        if count == 0{
            return 0;
        }
        return medians[(count / 2) as usize];
    }

    //Simple alg for now: check if ins_weight + del > non_ins_weight / 10. 
    // Take ~100 bp window around (or until confident site) -> SPOA.
    // Continue scanning after this 100bp window.
    pub fn spoa_consensus_test(self) -> Vec<u8>{
        let lapper = Lapper::new(self.overlaps.into_iter().map(|(start, stop, val)| Interval{start, stop, val}).collect());
        let mut i = 0;
        return vec![];
    }

    pub fn produce_consensus(self) -> Vec<u8> {

        let mut consensus_str = Vec::new();
        for consensus in self.consensus{
            log::trace!("POS: {}, INS: {:?}, COUNTS: {} {} {} {}, DEL: {}, INS_WEIGHT {} NO_INS_WEIGHT {}", consensus_str.len(), &consensus.insertions, consensus.a_count, consensus.c_count, consensus.g_count, consensus.t_count, consensus.deletion_count, consensus.prev_ins_weight, consensus.prev_nonins_weight);
            let mut max_count = 0;
            let mut max_base = b'N';
            if consensus.a_count > max_count {
                max_count = consensus.a_count;
                max_base = b'A';
            }
            if consensus.c_count > max_count {
                max_count = consensus.c_count;
                max_base = b'C';
            }
            if consensus.g_count > max_count {
                max_count = consensus.g_count;
                max_base = b'G';
            }
            if consensus.t_count > max_count {
                max_count = consensus.t_count;
                max_base = b'T';
            }
            if consensus.deletion_count > max_count {
                max_count = consensus.deletion_count;
                max_base = b'-';
            }
            let max_hpc_length = if !self.hpc { 1 } else {ConsensusBuilder::weighted_av_hpc(&consensus.homopolymer_lengths.unwrap())};
            
            let mut inserted_seq = None;
            if consensus.prev_ins_weight > consensus.prev_nonins_weight{
                let mut max_ins_count = 0;
                for (insert_seq, count) in consensus.insertions.sequences{
                    if count > max_ins_count {
                        inserted_seq = Some(insert_seq);
                        max_ins_count = count;
                    }
                }
            }
            if let Some(inserted_seq) = inserted_seq{
                let mut decompressed_ins_seq = SmallVec::new();
                if self.hpc{
                    for i in 0..inserted_seq.len(){
                        let hpc_lengths = &consensus.insertions.sequence_hpc_lengths.as_ref().unwrap()[&inserted_seq][i];
                        let max_hpc_length = ConsensusBuilder::weighted_av_hpc(hpc_lengths);
                        for _ in 0..max_hpc_length{
                            decompressed_ins_seq.push(inserted_seq[i]);
                        }
                    }
                }
                else{
                    decompressed_ins_seq = inserted_seq;
                }
                consensus_str.extend(decompressed_ins_seq);
            }
            if max_base != b'-' {
                if max_count == 0{
                    consensus_str.push(b'N');
                    continue
                }
                //println!("POS: {}, BASE: {}, HPC: {}", consensus_str.len(), max_base as char, max_hpc_length);
                for _ in 0..max_hpc_length{
                    consensus_str.push(max_base);
                }
            }

        }

        return consensus_str;
    }

    pub fn process_alignment_nohpc(&mut self, cigar: &Vec<OpLen>, query_sequence: &[u8], query_qualities: &[u8], start_ref: usize, start_query: usize) {
        let mut ref_pos = start_ref;
        let mut query_pos = start_query;
        let mut iter = cigar.iter();
        let mut prev_state = Operation::Sentinel;

        while let Some(oplen) = iter.next() {
            let prev_pos_query = if query_pos > 0 { query_pos - 1 } else { 0 };
            match oplen.op {
                Operation::Eq | Operation::X | Operation::M => {
                    if prev_state == Operation::I{
                        let quality_base = query_qualities[prev_pos_query];
                        self.consensus[ref_pos].prev_ins_weight = self.consensus[ref_pos].prev_ins_weight.saturating_add(quality_base as u32);
                    }
                    else{
                        let quality_base = query_qualities[prev_pos_query];
                        self.consensus[ref_pos].prev_nonins_weight = self.consensus[ref_pos].prev_nonins_weight.saturating_add(quality_base as u32);
                    }
                    //1 Match or mismatch
                    for ind in 0..oplen.len{
                        if ind > 0{
                            let quality_before = query_qualities[query_pos - 1];
                            self.consensus[ref_pos].prev_nonins_weight = self.consensus[ref_pos].prev_nonins_weight.saturating_add(quality_before as u32);
                        }

                        let quality_base = query_qualities[query_pos];
                        let query_base = query_sequence[query_pos];
                        let hpc_length = None;
                        self.consensus[ref_pos].increment_base(query_base, quality_base, hpc_length);
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                Operation::D => {
                    let quality_base_before = query_qualities[prev_pos_query];
                    if prev_state == Operation::I{
                        self.consensus[ref_pos].prev_ins_weight = self.consensus[ref_pos].prev_ins_weight.saturating_add(quality_base_before as u32);
                    }
                    else{
                        self.consensus[ref_pos].prev_nonins_weight = self.consensus[ref_pos].prev_nonins_weight.saturating_add(quality_base_before as u32);
                    }
                    // Deletion relative to reference
                    for ind in 0..oplen.len {
                        let quality_base_before = query_qualities[prev_pos_query];
                        if ind > 0{
                            self.consensus[ref_pos].prev_nonins_weight = self.consensus[ref_pos].prev_nonins_weight.saturating_add(quality_base_before as u32);
                        }
                        let query_pos = query_pos.min(query_qualities.len() - 1);
                        let quality_base = query_qualities[query_pos];
                        self.consensus[ref_pos].increment_deletion(quality_base);
                        ref_pos += 1;
                    }
                }
                Operation::I => {
                    // Insertion relative to reference
                    let insertion = &query_sequence[query_pos..query_pos + oplen.len];
                    let qualities = &query_qualities[query_pos..query_pos + oplen.len];
                    let hpc_lengths = None;
                    self.consensus[ref_pos].add_insertion(insertion, qualities, hpc_lengths);
                    query_pos += oplen.len;
                }
                _ => {
                    // Skip other CIGAR operations
                }
            }
            prev_state = oplen.op;
        }
    }

    pub fn process_alignment(&mut self, cigar: Vec<OpLen>, query_hpc_seq: &HomopolymerCompressedSeq, start_ref: usize, start_query: usize) {
        let mut ref_pos = start_ref;
        let mut query_pos = start_query;
        let mut iter = cigar.iter();
        let query_sequence = &query_hpc_seq.seq;
        let query_qualities = &query_hpc_seq.qualities;
        let query_hpc_lengths = &query_hpc_seq.homopolymer_lengths;
        let mut prev_state = Operation::Sentinel;

        while let Some(oplen) = iter.next() {
            let prev_pos_query = if query_pos > 0 { query_pos - 1 } else { 0 };
            match oplen.op {
                Operation::Eq | Operation::X | Operation::M => {
                    if prev_state == Operation::I{
                        let quality_base = query_qualities[prev_pos_query];
                        self.consensus[ref_pos].prev_ins_weight = self.consensus[ref_pos].prev_ins_weight.saturating_add(quality_base as u32);
                    }
                    else{
                        let quality_base = query_qualities[prev_pos_query];
                        self.consensus[ref_pos].prev_nonins_weight = self.consensus[ref_pos].prev_nonins_weight.saturating_add(quality_base as u32);
                    }
                    //1 Match or mismatch
                    for ind in 0..oplen.len{
                        if ind > 0{
                            let quality_before = query_qualities[query_pos - 1];
                            self.consensus[ref_pos].prev_nonins_weight = self.consensus[ref_pos].prev_nonins_weight.saturating_add(quality_before as u32);
                        }

                        let quality_base = query_qualities[query_pos];
                        let query_base = query_sequence[query_pos];
                        let hpc_length = if self.hpc { Some(query_hpc_lengths.as_ref().unwrap()[query_pos]) } else { None };
                        self.consensus[ref_pos].increment_base(query_base, quality_base, hpc_length);
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                Operation::D => {
                    let quality_base_before = query_qualities[prev_pos_query];
                    if prev_state == Operation::I{
                        self.consensus[ref_pos].prev_ins_weight = self.consensus[ref_pos].prev_ins_weight.saturating_add(quality_base_before as u32);
                    }
                    else{
                        self.consensus[ref_pos].prev_nonins_weight = self.consensus[ref_pos].prev_nonins_weight.saturating_add(quality_base_before as u32);
                    }
                    // Deletion relative to reference
                    for ind in 0..oplen.len {
                        let quality_base_before = query_qualities[prev_pos_query];
                        if ind > 0{
                            self.consensus[ref_pos].prev_nonins_weight = self.consensus[ref_pos].prev_nonins_weight.saturating_add(quality_base_before as u32);
                        }
                        let query_pos = query_pos.min(query_qualities.len() - 1);
                        let quality_base = query_qualities[query_pos];
                        self.consensus[ref_pos].increment_deletion(quality_base);
                        ref_pos += 1;
                    }
                }
                Operation::I => {
                    // Insertion relative to reference
                    let insertion = &query_sequence[query_pos..query_pos + oplen.len];
                    let qualities = &query_hpc_seq.qualities[query_pos..query_pos + oplen.len];
                    let hpc_lengths = if self.hpc { Some(&query_hpc_lengths.as_ref().unwrap()[query_pos..query_pos + oplen.len]) } else { None };
                    self.consensus[ref_pos].add_insertion(insertion, qualities, hpc_lengths);
                    query_pos += oplen.len;
                }
                _ => {
                    // Skip other CIGAR operations
                }
            }
            prev_state = oplen.op;
        }
    }
}

#[derive(Debug, Clone)]
struct ReadCigarIndex{
    pub read_index: usize,
    pub reference_index: usize,
    pub reverse: bool,
    pub cigar_index: Vec<OpLen>,
    pub alignment_start: usize,
}

impl PartialEq for ReadCigarIndex{
    fn eq(&self, other: &Self) -> bool {
        self.read_index == other.read_index && self.reference_index == other.reference_index && self.alignment_start == other.alignment_start
    }
}

impl Eq for ReadCigarIndex{


}


impl ReadCigarIndex{
    fn new(read_index: usize, reference_index: usize, reverse: bool, cigar_index: Vec<OpLen>, alignment_start: usize) -> Self {
        Self { read_index, reference_index, reverse, cigar_index, alignment_start}
    }

    // a c g t - - t t t a QUERY
    // a c g t a a t t - a REF
    // 0 1 2 3 4 5 6 7 8 9
    //the query substring from [2,5] is gt
    // the query substring from [2,6] is gtt
    // the query substring from [2,7] is gttt
    // the query substring from [5,9] is ttta
    fn get_query_substring_at_ref(&self, ref_start: usize, ref_end: usize) -> (usize, usize){
        let mut query_start = 0;
        let mut query_end = 0;
        let mut ref_pos = self.alignment_start;
        let mut query_pos = 0;
        let mut cigar_iter = self.cigar_index.iter();
        while let Some(oplen) = cigar_iter.next(){
            match oplen.op{
                Operation::Eq | Operation::X => {
                    for _ in 0..oplen.len{
                        if ref_pos == ref_start{
                            query_start = query_pos;
                        }
                        if ref_pos == ref_end{
                            query_end = query_pos;
                        }
                        ref_pos += 1;
                        query_pos += 1;
                    }
                }
                Operation::D => {
                    for _ in 0..oplen.len{
                        if ref_pos == ref_start{
                            query_start = query_pos;
                        }
                        if ref_pos == ref_end{
                            query_end = query_pos;
                        }
                        ref_pos += 1;
                    }
                }
                Operation::I => {
                    for _ in 0..oplen.len{
                        if ref_pos == ref_start{
                            query_start = query_pos;
                        }
                        if ref_pos == ref_end{
                            query_end = query_pos;
                        }
                        query_pos += 1;
                    }
                }
                _ => {}
            }
        }

        return (query_start, query_end);
    }
}