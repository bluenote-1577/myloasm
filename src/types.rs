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

use smallvec::SmallVec;
use fxhash::FxHashSet;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::{BuildHasherDefault, Hasher};
use bio_seq::prelude::*;
use rust_lapper::Lapper;

pub type Kmer64 = u64;
pub type Kmer32 = u32;
pub type KmerHash64 = u64;
pub type KmerHash32 = u32;
pub type Snpmer64 = u64;
pub type Splitmer64 = u64;

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
    pub minimizers: Vec<(usize, u64)>,
    pub snpmers: Vec<(usize, u64)>,
    pub id: String,
    pub k: u8,
    pub base_length: usize,
    pub dna_seq: Seq<Dna>,
    pub est_id: Option<f64>,
}


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Hash, Eq)]
pub struct SplitKmerInfo {
    pub full_kmer: Kmer64,
    pub mid_base: u8,
    pub canonical: bool,
    pub k: u8,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct KmerGlobalInfo {
    pub snpmer_info: Vec<SnpmerInfo>,
    pub solid_kmers: FxHashSet<Kmer64>,
    pub high_freq_thresh: f64,
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
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Hash, Eq)]
pub struct CountsAndBases{
    pub counts: SmallVec<[[u32;2];4]>,
    pub bases: SmallVec<[u8; 4]>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
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
    pub mean_depth: f64,
    pub mapping_boundaries: Lapper<u32, bool>,
    pub present: bool,
    pub length: usize,
    pub mapped_indices: Vec<usize>
}

pub trait NodeMapping {
    fn median_depth(&self) -> f64;
    fn mean_depth(&self) -> f64;
    fn mapping_boundaries(&self) -> &Lapper<u32, bool>;
    fn set_mapping_info(&mut self, mapping_info: MappingInfo);
    fn mapping_info_present(&self) -> bool;
    fn reference_length(&self) -> usize;
    fn mapped_indices(&self) -> &Vec<usize>;
}

#[derive(Debug, Clone, PartialEq)]
pub struct Breakpoints {
    pub pos1: usize,
    pub pos2: usize,
    pub cov: usize,
}