//Various byte-tables and hashing methods are taken from miniprot by Heng Li. Attached below is their license:
//The MIT License

// **** miniprot LICENSE ***
//Copyright (c) 2022-     Dana-Farber Cancer Institute
//
//Permission is hereby granted, free of charge, to any person obtaining
//a copy of this software and associated documentation files (the
//"Software"), to deal in the Software without restriction, including
//without limitation the rights to use, copy, modify, merge, publish,
//distribute, sublicense, and/or sell copies of the Software, and to
//permit persons to whom the Software is furnished to do so, subject to
//the following conditions:
//
//The above copyright notice and this permission notice shall be
//included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
//BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
//ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
//******************************

use block_aligner::cigar::Operation;
use smallvec::SmallVec;
use fxhash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::collections::HashSet;
use fxhash::FxHashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::path::PathBuf;
use bio_seq::prelude::*;
use rust_lapper::Lapper;
use block_aligner::cigar::OpLen;

use crate::constants::ID_THRESHOLD_ITERS;

pub type Kmer64 = u64;
pub type Kmer32 = u32;
pub type KmerHash64 = u64;
pub type KmerHash32 = u32;
pub type Snpmer64 = u64;
pub type Splitmer64 = u64;
pub type MultiCov = [f64; ID_THRESHOLD_ITERS];

pub const BYTE_TO_SEQ: [u8; 256] = [
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

#[inline]
pub fn mm_hash_64(key: u64) -> usize {
    let mut key = key;
    key = (!key).wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    return key as usize;
}

#[inline]
pub fn mm_hash(bytes: &[u8]) -> usize {
    let mut key = usize::from_ne_bytes(bytes.try_into().unwrap()) as usize;
    key = (!key).wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    return key;
}

pub struct MMHasher {
    hash: usize,
}

impl Hasher for MMHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        self.hash = mm_hash(bytes);
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.hash as u64
    }
}

impl Default for MMHasher {
    #[inline]
    fn default() -> MMHasher {
        MMHasher { hash: 0 }
    }
}

//Implement minimap2 hashing, will test later.
pub type MMBuildHasher = BuildHasherDefault<MMHasher>;
pub type MMHashMap<K, V> = HashMap<K, V, MMBuildHasher>;
pub type MMHashSet<K> = HashSet<K, MMBuildHasher>;

// Take a bit-encoded k-mer (k <= 32) and decode it as a string of ACGT

pub fn decode_kmer(kmer: Kmer64, k: u8) -> String {
    let mut seq = String::new();
    for i in 0..k {
        let c = (kmer >> (i * 2)) & 0b11;
        seq.push(match c {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => unreachable!(),
        });
    }
    //reverse string
    seq.chars().rev().collect()
}

#[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct PreFragment {
    pub kmers_with_refpos: Vec<(Kmer64, u64)>,
    pub upper_base: usize,
    pub lower_base: usize,
    pub id: String,
}

#[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct VarmerFragment {
    pub upper: usize,
    pub lower: usize,
    pub varmers: FxHashSet<usize>,
    pub upper_base: usize,
    pub lower_base: usize,
    pub id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct Varmer {
    pub kmer: Kmer64,
    pub count: u32,
    pub pos: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, Default)]
pub struct BasePileup {
    pub ref_pos: u64,
    pub ref_base: u8,
    pub base_freqs: [u32; 4],
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, Default)]
pub struct TigdexOverlap {
    pub tig1: usize,
    pub tig2: usize,
    pub tig1_start: usize,
    pub tig1_end: usize,
    pub tig2_start: usize,
    pub tig2_end: usize,
    pub shared_tig: usize,
    pub variable_roots: usize,
    pub variable_tigs: usize,
    pub chain_reverse: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct ChainInfo {
    pub chain: Vec<Anchor>,
    pub reverse: bool,
    pub score: i32,
}

pub type EdgeIndex = usize;
pub type NodeIndex = usize;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct TigRead {
    pub tig_seq: Vec<u32>,
    pub id: String,
}

//TODO: we restructure minimizers and snpmmer positions to another vector
//and change to u32 to save space
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct TwinRead {
    //pub minimizers: Vec<(usize, u64)>,
    pub minimizer_positions: Vec<u32>,
    pub minimizer_kmers: Vec<u64>,
    //pub snpmers: Vec<(usize, u64)>,
    pub snpmer_kmers: Vec<u64>,
    pub snpmer_positions: Vec<u32>,
    pub id: String,
    pub base_id: String,
    pub k: u8,
    pub base_length: usize,
    pub dna_seq: Seq<Dna>,
    pub qual_seq: Option<Seq<QualCompact3>>,
    pub est_id: Option<f64>,
    pub min_depth_multi: Option<MultiCov>,
    pub median_depth: Option<f64>,
    pub split_chimera: bool,
    pub split_start: u32,
    pub outer: bool,
    pub snpmer_id_threshold: Option<f64>,
}


#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum QualCompact3{
    Q33 = 0b0000,
    Q36 = 0b0001,
    Q39 = 0b0010,
    Q42 = 0b0011,
    Q45 = 0b0100,
    Q48 = 0b0101,
    Q51 = 0b0110,
    Q54 = 0b0111,
    Q57 = 0b1000,
    Q60 = 0b1001,
    Q63 = 0b1010,
    Q66 = 0b1011,
    Q69 = 0b1100,
    Q72 = 0b1101,
    Q75 = 0b1110,
    Q78 = 0b1111,
}

impl Codec for QualCompact3{
    const BITS: u8 = 4;

    /// Take the two least significant bits of a `u8` and map them to the
    /// corresponding nucleotides.
    fn unsafe_from_bits(b: u8) -> Self {
        unsafe { std::mem::transmute(b & 0b1111) }
    }

    /// We can efficient verify that a byte is a valid `Dna` value if it's
    /// between 0 and 3.
    fn try_from_bits(b: u8) -> Option<Self> {

        // Round to nearest 3 and map to enum variant
        let rounded = match b {
            0..=34  => 0,  // Q33
            35..=37 => 1,  // Q36
            38..=40 => 2,  // Q39
            41..=43 => 3,  // Q42
            44..=46 => 4,  // Q45
            47..=49 => 5,  // Q48
            50..=52 => 6,  // Q51
            53..=55 => 7,  // Q54
            56..=58 => 8,  // Q57
            59..=61 => 9,  // Q60
            62..=64 => 10, // Q63
            65..=67 => 11, // Q66
            68..=70 => 12, // Q69
            71..=73 => 13, // Q72
            74..=76 => 14, // Q75
            _ => 15,       // Q78 or higher
        };

        // Use match instead of transmute for safety
        let m = match rounded {
            0 => Some(Self::Q33),
            1 => Some(Self::Q36),
            2 => Some(Self::Q39),
            3 => Some(Self::Q42),
            4 => Some(Self::Q45),
            5 => Some(Self::Q48),
            6 => Some(Self::Q51),
            7 => Some(Self::Q54),
            8 => Some(Self::Q57),
            9 => Some(Self::Q60),
            10 => Some(Self::Q63),
            11 => Some(Self::Q66),
            12 => Some(Self::Q69),
            13 => Some(Self::Q72),
            14 => Some(Self::Q75),
            15 => Some(Self::Q78),
            _ => None,  // This case should never happen given our match above
        };

        return m
    }

    /// The ASCII values of 'A', 'C', 'G', and 'T' can be translated into
    /// the numbers 0, 1, 2, and 3 using bitwise operations: `((b << 1) + b) >> 3`.
    fn unsafe_from_ascii(b: u8) -> Self {
        Self::unsafe_from_bits(b)
    }

    fn try_from_ascii(c: u8) -> Option<Self> {
        Self::try_from_bits(c)
    }

    //Not -33
    fn to_char(self) -> char {
        match self {
            QualCompact3::Q33 => '!',
            QualCompact3::Q36 => '"',
            QualCompact3::Q39 => '#',
            QualCompact3::Q42 => '$',
            QualCompact3::Q45 => '%',
            QualCompact3::Q48 => '&',
            QualCompact3::Q51 => '\'',
            QualCompact3::Q54 => '(',
            QualCompact3::Q57 => ')',
            QualCompact3::Q60 => '*',
            QualCompact3::Q63 => '+',
            QualCompact3::Q66 => ',',
            QualCompact3::Q69 => '-',
            QualCompact3::Q72 => '.',
            QualCompact3::Q75 => '/',
            QualCompact3::Q78 => '0',
        }
    }

    fn to_bits(self) -> u8 {
        self as u8
    }

    fn items() -> impl Iterator<Item = Self> {
        vec![
            QualCompact3::Q33,
            QualCompact3::Q36,
            QualCompact3::Q39,
            QualCompact3::Q42,
            QualCompact3::Q45,
            QualCompact3::Q48,
            QualCompact3::Q51,
            QualCompact3::Q54,
            QualCompact3::Q57,
            QualCompact3::Q60,
            QualCompact3::Q63,
            QualCompact3::Q66,
            QualCompact3::Q69,
            QualCompact3::Q72,
            QualCompact3::Q75,
            QualCompact3::Q78,
        ].into_iter()
    }
}

impl Complement for QualCompact3 {
    fn comp(&self) -> Self {
        return *self 
    }
}

impl TwinRead{
    pub fn clear(&mut self){
        self.minimizer_kmers.clear();
        self.snpmer_kmers.clear();
        self.minimizer_positions.clear();
        self.snpmer_positions.clear();
    }

    pub fn minimizers(&self) -> impl Iterator<Item = (u32, u64)> + '_ {
        assert!(self.minimizer_positions.len() == self.minimizer_kmers.len());
        self.minimizer_positions.iter().zip(self.minimizer_kmers.iter()).map(|(x, y)| (*x, *y))
    }

    pub fn minimizers_vec(&self) -> Vec<(u32, u64)> {
        assert!(self.minimizer_positions.len() == self.minimizer_kmers.len());
        self.minimizer_positions.iter().zip(self.minimizer_kmers.iter()).map(|(x, y)| (*x, *y)).collect()
    }

    pub fn snpmers(&self) -> impl Iterator<Item = (u32, u64)> + '_ {
        assert!(self.snpmer_positions.len() == self.snpmer_kmers.len());
        self.snpmer_positions.iter().zip(self.snpmer_kmers.iter()).map(|(x, y)| (*x, *y))
    }

    pub fn retain_mini_positions(&mut self, positions: FxHashSet<usize>) {
        retain_vec_positions(&mut self.minimizer_kmers, &positions);
        retain_vec_positions(&mut self.minimizer_positions, &positions);
        self.minimizer_kmers.shrink_to_fit();
        self.minimizer_positions.shrink_to_fit();
    }

    pub fn retain_snpmer_positions(&mut self, positions: FxHashSet<usize>) {
        retain_vec_positions(&mut self.snpmer_kmers, &positions);
        retain_vec_positions(&mut self.snpmer_positions, &positions);
        self.snpmer_kmers.shrink_to_fit();
        self.snpmer_positions.shrink_to_fit();
    }

    
    pub fn shift_and_retain(&mut self, other_read: &TwinRead, last_break: usize, bp_start: usize, k: usize){

        // new_read.minimizers = twin_read.minimizers.iter().filter(|x| x.0 >= last_break && x.0 + k - 1 < bp_start).copied().map(|x| (x.0 - last_break, x.1)).collect();
        // new_read.snpmers = twin_read.snpmers.iter().filter(|x| x.0 >= last_break && x.0 + k - 1 < bp_start).copied().map(|x| (x.0 - last_break, x.1)).collect();
        // new_read.minimizers.shrink_to_fit();
        // new_read.snpmers.shrink_to_fit();

        let mut mini_positions_filtered = Vec::new();
        let mut mini_kmers_filtered = Vec::new();
        let mut snp_positions_filtered = Vec::new();
        let mut snp_kmers_filtered = Vec::new();

        for (pos, kmer) in other_read.minimizers(){
            if pos >= last_break as u32 && pos + k as u32 - 1 < bp_start as u32{
                mini_positions_filtered.push(pos - last_break as u32);
                mini_kmers_filtered.push(kmer);
            }
        }

        for (pos, kmer) in other_read.snpmers(){
            if pos >= last_break as u32 && pos + k as u32 - 1 < bp_start as u32{
                snp_positions_filtered.push(pos - last_break as u32);
                snp_kmers_filtered.push(kmer);
            }
        }

        //shirnk to fit
        mini_positions_filtered.shrink_to_fit();
        mini_kmers_filtered.shrink_to_fit();
        snp_positions_filtered.shrink_to_fit();
        snp_kmers_filtered.shrink_to_fit();

        self.minimizer_positions = mini_positions_filtered;
        self.minimizer_kmers = mini_kmers_filtered;
        self.snpmer_positions = snp_positions_filtered;
        self.snpmer_kmers = snp_kmers_filtered;
    }
}

pub fn retain_vec_positions<T>(vec: &mut Vec<T>, positions: &FxHashSet<usize>){
    let mut i = 0;
    vec.retain(|_| {
        let keep = positions.contains(&i);
        i += 1;
        keep
    });
}


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct KmerGlobalInfo {
    pub snpmer_info: Vec<SnpmerInfo>,
    pub solid_kmers: HashSet<Kmer64>,
    pub high_freq_thresh: f64,
    pub read_files: Vec<PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct TwinReadContainer {
    pub twin_reads: Vec<TwinRead>,
    pub outer_indices: Vec<usize>,
    //Not implemented TODO
    pub tig_reads: Vec<TwinRead>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Eq, Ord, PartialOrd, Hash)]
pub struct SnpmerInfo {
    pub split_kmer: u64,
    pub mid_bases: SmallVec<[u8;2]>,
    pub counts: SmallVec<[u32;2]>,
    pub k: u8,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub struct TwinOverlap{
    pub i1: usize,
    pub i2: usize,
    pub start1: usize,
    pub end1: usize,
    pub start2: usize,
    pub end2: usize,
    pub shared_minimizers: usize,
    pub shared_snpmers: usize,
    pub snpmers_in_both: (usize, usize),
    pub diff_snpmers: usize,
    pub chain_reverse: bool,
    pub intersect: (usize, usize),
    pub chain_score: i32,
    pub minimizer_chain: Option<Vec<Anchor>>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub struct SnpmerHit {
    pub pos1: u32,
    pub pos2: u32,
    pub bases: (u8, u8)
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Hash, Eq)]
pub struct CountsAndBases{
    pub counts: SmallVec<[[u32;2];2]>,
    pub bases: SmallVec<[u8; 4]>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Eq, Hash)]
pub struct Anchor{
    pub i: u32,
    pub j: u32,
    pub pos1: u32,
    pub pos2: u32,
}

// Enum for marking the state of a node during processing
#[derive(PartialEq, Eq, Clone, Debug, Hash, Copy)]
pub enum Direction {
    Incoming,
    Outgoing
}
impl Direction{
    pub fn reverse(&self) -> Direction{
        match self{
            Direction::Incoming => Direction::Outgoing,
            Direction::Outgoing => Direction::Incoming
        }
    }
}

#[inline]
pub fn bits_to_ascii(bit_rep: u8) -> u8{
    match bit_rep{
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => unreachable!()
    }
}


#[derive(Debug, Clone)]
pub struct MappingInfo {
    pub median_depth: f64,
    pub minimum_depth: f64,
    pub mapping_boundaries: Lapper<u32, SmallTwinOl>,
    //pub kmer_counts: Vec<u32>,
    pub present: bool,
    pub length: usize,
}

pub trait NodeMapping {
    fn median_mapping_depth(&self) -> f64;
    fn min_mapping_depth(&self) -> f64;
    fn mapping_boundaries(&self) -> &Lapper<u32, SmallTwinOl>;
    fn set_mapping_info(&mut self, mapping_info: MappingInfo);
    fn mapping_info_present(&self) -> bool;
    fn reference_length(&self) -> usize;
    fn mapped_indices(&self) -> Vec<usize>;
}

#[derive(Debug, Clone,  Default)]
pub struct SmallTwinOl{
    pub query_id: u32,
    pub snpmer_identity: f32,
    pub reverse: bool,
    pub maximal_overlap: bool,
    pub alignment_result: Option<AlignmentResult>
}

impl Eq for SmallTwinOl{}

impl PartialEq for SmallTwinOl{
    fn eq(&self, other: &Self) -> bool{
        self.query_id == other.query_id && self.snpmer_identity == other.snpmer_identity && self.reverse == other.reverse && self.maximal_overlap == other.maximal_overlap
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Breakpoints {
    pub pos1: usize,
    pub pos2: usize,
    pub cov: usize,
}
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GetSequenceInfoConfig{
    pub blunted: bool,
    pub dna_seq_info: bool,
    pub best_overlap_chunk: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BubblePopResult{
    pub original_direction: Direction,
    pub end_direction: Direction,
    pub source_hash_id: NodeIndex,
    pub sink_hash_id: NodeIndex,
    pub remove_nodes: Vec<NodeIndex>,
    pub remove_edges: FxHashSet<EdgeIndex>,
}

impl BubblePopResult{
    pub fn new(original_direction: Direction, end_direction: Direction, source_hash_id: NodeIndex, sink_hash_id: NodeIndex, remove_nodes: Vec<NodeIndex>, remove_edges: FxHashSet<EdgeIndex>) -> Self{
        BubblePopResult{
            original_direction,
            end_direction,
            source_hash_id,
            sink_hash_id,
            remove_nodes,
            remove_edges
        }
    }
}


#[derive(Debug, Clone, PartialEq, Default)]
pub struct BeamSearchSoln{
    pub path: Vec<EdgeIndex>,
    pub coverages: Vec<(MultiCov, usize)>, 
    pub score: f64,
    pub path_nodes: Vec<NodeIndex>,
    pub depth: usize,
    pub current_length: usize
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct OverlapAdjMap {
    pub adj_map: FxHashMap<NodeIndex, Vec<NodeIndex>>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct OpLenVec{
    pub op_vec: Vec<Operation>,
    pub len_vec: Vec<u32>,
}

impl OpLenVec{
    pub fn new(cigar_vec: Vec<OpLen>) -> Self{
        let mut op_vec = Vec::new();
        let mut len_vec = Vec::new();
        let mut last_op = Operation::Sentinel;
        for op_len in cigar_vec{
            if op_len.op == last_op{
                *len_vec.last_mut().unwrap() += op_len.len as u32;
            }
            else if (op_len.op == Operation::X || op_len.op == Operation::Eq) && last_op == Operation::M{
                *len_vec.last_mut().unwrap() += op_len.len as u32;
            }
            else{
                op_vec.push(op_len.op);
                len_vec.push(op_len.len as u32);
                last_op = op_len.op;
            }
            if last_op == Operation::X || last_op == Operation::Eq{
                last_op = Operation::M;
            }
        }
        assert!(op_vec.len() == len_vec.len());
        OpLenVec{
            op_vec,
            len_vec,
        }
    }

    pub fn len(&self) -> usize{
        return self.op_vec.len();
    }

    pub fn iter(&self) -> impl Iterator<Item = (Operation, u32)> + '_ {
        self.op_vec.iter()
            .zip(self.len_vec.iter())
            .map(|(op, len)| (*op, *len))
    }
}


#[derive(Debug, Clone, Default)]
pub struct AlignmentResult{
    pub cigar: OpLenVec,
    pub q_start: usize,
    pub q_end: usize,
    pub r_start: usize,
    pub r_end: usize,
}

pub fn dna_seq_to_u8(slice: &Seq<Dna>) -> Vec<u8>{
    slice.iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect()
}

pub fn dna_slice_to_u8(slice: &SeqSlice<Dna>) -> Vec<u8>{
    slice.iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect()
}


pub fn quality_slice_to_u8(slice: &SeqSlice<QualCompact3>) -> Vec<u8>{
    slice.iter().map(|x| x as u8 * 3).collect()
}

pub fn quality_seq_to_u8(slice: &Seq<QualCompact3>) -> Vec<u8>{
    slice.iter().map(|x| x as u8 * 3).collect()
}

pub fn revcomp_u8(seq: &Vec<u8>) -> Vec<u8>{
    seq.iter().rev().map(|x| match x{
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N'
    }).collect()
}