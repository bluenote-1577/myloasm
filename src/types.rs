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
use std::cmp::Ordering;
use nibble_vec::*;
use std::collections::BTreeMap;
use minimum_redundancy::{Coding, Code, DecodingResult, BitsPerFragment};
use std::sync::OnceLock;


use crate::constants::ID_THRESHOLD_ITERS;
use crate::constants::MAX_GAP_CHAINING;

pub type NodeMap<K,V> = BTreeMap<K,V>;
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


#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Kmer48{
    data: [u8; 6],
}

impl Kmer48{
    // Be careful about endian
    #[inline]
    pub fn from_u64(n: u64) -> Self {
        let bytes = n.to_le_bytes();

        debug_assert!(bytes[6..8].iter().all(|&b| b == 0));

        Self {
            data: [bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5]],
        }
    }

    // Assume Kmer48 is stored as little endian. So AAA = ...0000000_11_11_11
    pub fn to_u64(self) -> u64 {
        let mut bytes = [0; 8];
        bytes[0..6].copy_from_slice(&self.data);
        u64::from_le_bytes(bytes)
    }
}

impl From<u64> for Kmer48 {
    fn from(value: u64) -> Self {
        Kmer48::from_u64(value)
    }
}

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

pub fn decode_kmer64(kmer: Kmer64, k: u8) -> String {
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

pub fn decode_kmer48(kmer: Kmer48, k: u8) -> String {
    let kmer = kmer.to_u64();
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
    pub large_indel: bool
}

pub type EdgeIndex = usize;
pub type NodeIndex = usize;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct TigRead {
    pub tig_seq: Vec<u32>,
    pub id: String,
}


pub type Percentage = f64;
pub type Fraction = f32;

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct TwinRead {
    // Delta-encoded positions (varint or Huffman-coded varint)
    pub minimizer_positions_enc: Vec<u8>,
    pub snpmer_positions_enc: Vec<u8>,
    // True if position vecs are Huffman-encoded, false if raw varint
    pub huffman_encoded: bool,
    // Full record id
    pub id: String,
    // First word of record id
    pub base_id: String,
    pub k: u8,
    pub base_length: usize,
    pub dna_seq: Seq<Dna>,
    pub qual_seq: Option<Seq<QualCompact3>>,
    pub est_id: Option<Percentage>,
    pub min_depth_multi: Option<MultiCov>,
    pub median_depth: Option<f64>,
    pub split_chimera: bool,
    pub split_start: u32,
    pub outer: bool,
    pub snpmer_id_threshold: Option<f64>,
    pub overlap_hang_length: Option<(usize, usize)>,
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

impl ComplementMut for QualCompact3{
    fn comp(&mut self) {
    }
}

impl Complement for QualCompact3 {}

#[inline]
fn reverse_bit_pairs(n: u64, k: usize) -> u64 {
    let even_mask : u64 =  0xAAAAAAAAAAAAAAAAu64;
    let odd_mask: u64 =  0x5555555555555555u64;

    let odd_bits_rev = (n & odd_mask).reverse_bits() >> (64 - 2 * k);
    let even_bits_rev = (n & even_mask).reverse_bits() >> (64 - 2 * k);

    return (odd_bits_rev >> 1) | (even_bits_rev << 1);

}

/// Encode absolute positions as variable-length delta encoding.
/// Uses continuation-bit scheme: high bit set = more bytes follow.
/// Supports full u32 range with up to 5 bytes per delta.
/// Returns the encoded byte vector.
#[inline]
pub fn encode_positions_delta(positions: &[u32]) -> Vec<u8> {
    if positions.is_empty() {
        return Vec::new();
    }

    let mut encoded = Vec::with_capacity(positions.len() * 2); // rough estimate
    let mut prev_pos = 0u32;

    for &pos in positions {
        let mut delta = pos - prev_pos;

        // Encode delta using variable-length encoding (continuation bit scheme)
        // Each byte holds 7 bits of data, high bit is continuation flag
        loop {
            let byte = (delta & 0x7F) as u8;
            delta >>= 7;

            if delta == 0 {
                // Last byte: no continuation bit
                encoded.push(byte);
                break;
            } else {
                // More bytes to follow: set continuation bit
                encoded.push(byte | 0x80);
            }
        }

        prev_pos = pos;
    }

    encoded
}

/// Decode variable-length delta-encoded positions back to absolute positions.
/// Returns a vector of absolute u32 positions.
#[inline]
pub fn decode_positions_delta(encoded: &[u8]) -> Vec<u32> {
    let mut positions = Vec::new();
    let mut current_pos = 0u32;
    let mut i = 0;

    while i < encoded.len() {
        let mut delta = 0u32;
        let mut shift = 0;

        loop {
            let byte = encoded[i];
            i += 1;

            delta |= ((byte & 0x7F) as u32) << shift;

            if byte & 0x80 == 0 {
                // No continuation bit, we're done
                break;
            }

            shift += 7;
        }

        current_pos += delta;
        positions.push(current_pos);
    }

    positions
}

/// Encode a u32 value as variable-length bytes using continuation-bit scheme.
/// High bit set = more bytes follow. Returns number of bytes written.
#[inline]
pub fn encode_varint_u32(value: u32, buf: &mut Vec<u8>) {
    let mut v = value;
    loop {
        let byte = (v & 0x7F) as u8;
        v >>= 7;
        if v == 0 {
            buf.push(byte);
            break;
        } else {
            buf.push(byte | 0x80);
        }
    }
}

/// Decode a variable-length encoded u32 from a byte slice starting at given index.
/// Returns (value, bytes_consumed).
#[inline]
pub fn decode_varint_u32(buf: &[u8], start: usize) -> (u32, usize) {
    let mut value = 0u32;
    let mut shift = 0;
    let mut i = start;

    loop {
        let byte = buf[i];
        i += 1;
        value |= ((byte & 0x7F) as u32) << shift;
        if byte & 0x80 == 0 {
            break;
        }
        shift += 7;
    }

    (value, i - start)
}

/// Encode a slice of u32 values as variable-length bytes.
#[inline]
pub fn encode_varints(values: &[u32]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(values.len() * 2);
    for &v in values {
        encode_varint_u32(v, &mut buf);
    }
    buf
}

// ---- Huffman coding for position delta bytes ----

#[derive(Clone, Copy)]
pub enum PositionKind { Minimizer, Snpmer }

static HUFFMAN_CODING_MINI: OnceLock<Coding<u8, BitsPerFragment>> = OnceLock::new();
static HUFFMAN_CODES_MINI: OnceLock<[Code; 256]> = OnceLock::new();
static HUFFMAN_CODING_SNP: OnceLock<Coding<u8, BitsPerFragment>> = OnceLock::new();
static HUFFMAN_CODES_SNP: OnceLock<[Code; 256]> = OnceLock::new();

fn huffman_coding_for(kind: PositionKind) -> Option<&'static Coding<u8, BitsPerFragment>> {
    match kind {
        PositionKind::Minimizer => HUFFMAN_CODING_MINI.get(),
        PositionKind::Snpmer => HUFFMAN_CODING_SNP.get(),
    }
}

fn huffman_codes_for(kind: PositionKind) -> Option<&'static [Code; 256]> {
    match kind {
        PositionKind::Minimizer => HUFFMAN_CODES_MINI.get(),
        PositionKind::Snpmer => HUFFMAN_CODES_SNP.get(),
    }
}

/// Build separate Huffman codings for minimizer and snpmer delta bytes,
/// sampling from the first `max_reads` reads. Stores them globally.
///
/// All 256 byte values are given a minimum frequency of 1 so that
/// unseen bytes still get valid (long) codewords instead of zero-length codes.
pub fn build_and_set_huffman_coding(reads: &[TwinRead], max_reads: usize) {
    use std::collections::HashMap;
    use minimum_redundancy::Frequencies;

    let sample = &reads[..max_reads.min(reads.len())];

    // Build minimizer Huffman coding
    let mut mini_freq = HashMap::<u8, usize>::new();
    for b in 0..=255u8 { mini_freq.insert(b, 1); }
    mini_freq.add_occurences_of(sample.iter()
        .flat_map(|r| r.minimizer_positions_enc.iter().copied()));
    let coding_mini = Coding::from_frequencies(BitsPerFragment(1), mini_freq);
    HUFFMAN_CODES_MINI.set(coding_mini.codes_for_values_array()).ok();
    HUFFMAN_CODING_MINI.set(coding_mini).ok();

    // Build snpmer Huffman coding
    let mut snp_freq = HashMap::<u8, usize>::new();
    for b in 0..=255u8 { snp_freq.insert(b, 1); }
    snp_freq.add_occurences_of(sample.iter()
        .flat_map(|r| r.snpmer_positions_enc.iter().copied()));
    let coding_snp = Coding::from_frequencies(BitsPerFragment(1), snp_freq);
    HUFFMAN_CODES_SNP.set(coding_snp.codes_for_values_array()).ok();
    HUFFMAN_CODING_SNP.set(coding_snp).ok();
}

/// Returns true if Huffman coding has been initialized.
pub fn huffman_initialized() -> bool {
    HUFFMAN_CODING_MINI.get().is_some() && HUFFMAN_CODING_SNP.get().is_some()
}

/// Encode positions: uses Huffman if initialized, otherwise varint delta.
#[inline]
pub fn encode_positions(positions: &[u32], kind: PositionKind) -> Vec<u8> {
    if let Some(codes) = huffman_codes_for(kind) {
        encode_positions_huffman(positions, codes)
    } else {
        encode_positions_delta(positions)
    }
}

/// Decode positions: dispatches based on per-read huffman flag.
#[inline]
pub fn decode_positions(encoded: &[u8], kind: PositionKind, huffman_encoded: bool) -> Vec<u32> {
    if huffman_encoded {
        let coding = huffman_coding_for(kind)
            .expect("huffman_encoded flag set but Huffman coding not initialized");
        decode_positions_huffman(encoded, coding)
    } else {
        decode_positions_delta(encoded)
    }
}

/// Count positions without full decode: dispatches based on per-read huffman flag.
#[inline]
pub fn count_positions(encoded: &[u8], huffman_encoded: bool) -> usize {
    if huffman_encoded {
        count_positions_huffman(encoded)
    } else {
        encoded.iter().filter(|&&b| b & 0x80 == 0).count()
    }
}

/// Encode positions as Huffman-coded varint bytes.
/// Format: [num_positions as varint] [packed Huffman bits of varint byte stream]
pub fn encode_positions_huffman(positions: &[u32], codes: &[Code; 256]) -> Vec<u8> {
    if positions.is_empty() {
        return Vec::new();
    }

    let mut result = Vec::new();
    encode_varint_u32(positions.len() as u32, &mut result);

    let varint_bytes = encode_positions_delta(positions);

    let mut current_byte: u8 = 0;
    let mut bits_filled: u32 = 0;

    for &vbyte in &varint_bytes {
        let code = codes[vbyte as usize];
        for bit_idx in (0..code.len).rev() {
            let bit = ((code.content >> bit_idx) & 1) as u8;
            current_byte = (current_byte << 1) | bit;
            bits_filled += 1;
            if bits_filled == 8 {
                result.push(current_byte);
                current_byte = 0;
                bits_filled = 0;
            }
        }
    }

    if bits_filled > 0 {
        current_byte <<= 8 - bits_filled;
        result.push(current_byte);
    }

    result
}

/// Decode Huffman-encoded positions back to absolute positions.
pub fn decode_positions_huffman(encoded: &[u8], coding: &Coding<u8, BitsPerFragment>) -> Vec<u32> {
    if encoded.is_empty() {
        return Vec::new();
    }

    let (num_positions, prefix_len) = decode_varint_u32(encoded, 0);
    let num_positions = num_positions as usize;

    if num_positions == 0 {
        return Vec::new();
    }

    let huffman_bytes = &encoded[prefix_len..];
    let mut bit_iter = huffman_bytes.iter()
        .flat_map(|&byte| (0..8u32).rev().map(move |i| ((byte >> i) & 1) as u32));

    let mut decoder = coding.decoder();
    let mut positions = Vec::with_capacity(num_positions);
    let mut current_pos = 0u32;
    let mut delta = 0u32;
    let mut shift = 0u32;

    while positions.len() < num_positions {
        match decoder.decode_next(&mut bit_iter) {
            DecodingResult::Value(&vbyte) => {
                delta |= ((vbyte & 0x7F) as u32) << shift;
                if vbyte & 0x80 == 0 {
                    current_pos += delta;
                    positions.push(current_pos);
                    delta = 0;
                    shift = 0;
                } else {
                    shift += 7;
                }
            }
            _ => break,
        }
    }

    positions
}

/// Count positions from Huffman-encoded data by reading the varint prefix.
#[inline]
pub fn count_positions_huffman(encoded: &[u8]) -> usize {
    if encoded.is_empty() {
        return 0;
    }
    let (count, _) = decode_varint_u32(encoded, 0);
    count as usize
}

/// Re-encode reads from varint to Huffman encoding and set the flag.
/// Must be called after build_and_set_huffman_coding.
pub fn reencode_reads_huffman(reads: &mut [TwinRead]) {
    let codes_mini = HUFFMAN_CODES_MINI.get().expect("Huffman coding not initialized");
    let codes_snp = HUFFMAN_CODES_SNP.get().expect("Huffman coding not initialized");

    for read in reads.iter_mut() {
        if read.huffman_encoded {
            continue;
        }
        if !read.minimizer_positions_enc.is_empty() {
            let positions = decode_positions_delta(&read.minimizer_positions_enc);
            let mut enc = encode_positions_huffman(&positions, codes_mini);
            enc.shrink_to_fit();
            read.minimizer_positions_enc = enc;
        }
        if !read.snpmer_positions_enc.is_empty() {
            let positions = decode_positions_delta(&read.snpmer_positions_enc);
            let mut enc = encode_positions_huffman(&positions, codes_snp);
            enc.shrink_to_fit();
            read.snpmer_positions_enc = enc;
        }
        read.huffman_encoded = true;
    }
}

/// Save Huffman coding tables to a file. Call after build_and_set_huffman_coding.
pub fn save_huffman_tables(path: &std::path::Path) -> std::io::Result<()> {
    let coding_mini = HUFFMAN_CODING_MINI.get()
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::Other, "Huffman coding not initialized"))?;
    let coding_snp = HUFFMAN_CODING_SNP.get()
        .ok_or_else(|| std::io::Error::new(std::io::ErrorKind::Other, "Huffman coding not initialized"))?;

    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    coding_mini.write(&mut file, |w, v| { w.write_all(&[*v]) })?;
    coding_snp.write(&mut file, |w, v| { w.write_all(&[*v]) })?;
    Ok(())
}

/// Load Huffman coding tables from a file and set the global state.
/// Returns true if successfully loaded, false if file doesn't exist.
pub fn load_huffman_tables(path: &std::path::Path) -> std::io::Result<bool> {
    if !path.exists() {
        return Ok(false);
    }

    let mut file = std::io::BufReader::new(std::fs::File::open(path)?);

    let coding_mini: Coding<u8, BitsPerFragment> =
        Coding::read(&mut file, |r| { let mut b = [0u8; 1]; r.read_exact(&mut b)?; Ok(b[0]) })?;
    let coding_snp: Coding<u8, BitsPerFragment> =
        Coding::read(&mut file, |r| { let mut b = [0u8; 1]; r.read_exact(&mut b)?; Ok(b[0]) })?;

    HUFFMAN_CODES_MINI.set(coding_mini.codes_for_values_array()).ok();
    HUFFMAN_CODING_MINI.set(coding_mini).ok();
    HUFFMAN_CODES_SNP.set(coding_snp.codes_for_values_array()).ok();
    HUFFMAN_CODING_SNP.set(coding_snp).ok();

    log::info!("Loaded Huffman coding tables from {:?}", path);
    Ok(true)
}

// ---- End Huffman coding section ----

/// Encode a boundary pair (start, end) as (start, gap) using variable-length encoding.
/// Appends to the provided buffer.
#[inline]
pub fn encode_boundary_pair(start: u32, end: u32, buf: &mut Vec<u8>) {
    let gap = end.saturating_sub(start);
    encode_varint_u32(start, buf);
    encode_varint_u32(gap, buf);
}

/// Iterator over decoded boundary pairs from varint-encoded bytes.
pub struct BoundaryPairIter<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> BoundaryPairIter<'a> {
    pub fn new(buf: &'a [u8]) -> Self {
        BoundaryPairIter { buf, pos: 0 }
    }
}

impl<'a> Iterator for BoundaryPairIter<'a> {
    type Item = (u32, u32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.buf.len() {
            return None;
        }

        let (start, consumed1) = decode_varint_u32(self.buf, self.pos);
        self.pos += consumed1;

        let (gap, consumed2) = decode_varint_u32(self.buf, self.pos);
        self.pos += consumed2;

        Some((start, start + gap))
    }
}

impl TwinRead{

    pub fn compact(&mut self){

        // debug capacities vs lengths
        //  log::info!("BEFORE Minimizer positions: capacity = {}, length = {}", self.minimizer_positions_enc.capacity(), self.minimizer_positions_enc.len());
        //  log::info!("BEFORE Snpmer positions: capacity = {}, length = {}", self.snpmer_positions_enc.capacity(), self.snpmer_positions_enc.len());
        //  log::info!("BEFORE DNA seq: capacity = {}, length = {}", self.dna_seq.capacity(), self.dna_seq.len());
        //  log::info!("BEFORE Qual seq: capacity = {}, length = {}", if let Some(qual_seq) = self.qual_seq.as_ref(){qual_seq.capacity()} else {0}, if let Some(qual_seq) = self.qual_seq.as_ref(){qual_seq.len()} else {0});

        let mut minimizer_positions_enc = Vec::with_capacity(self.minimizer_positions_enc.len());
        minimizer_positions_enc.extend(self.minimizer_positions_enc.iter().cloned());
        self.minimizer_positions_enc = minimizer_positions_enc;

        let mut snpmer_positions_enc = Vec::with_capacity(self.snpmer_positions_enc.len());
        snpmer_positions_enc.extend(self.snpmer_positions_enc.iter().cloned());
        self.snpmer_positions_enc = snpmer_positions_enc;

        //self.dna_seq.shrink_to_exact();
        self.dna_seq.shrink_to_fit();
        if let Some(qual_seq) = self.qual_seq.as_mut(){
            // qual_seq.shrink_to_exact();
             qual_seq.shrink_to_fit();
        }

        // debug capacities vs lengths
        //  log::info!("Minimizer positions: capacity = {}, length = {}", self.minimizer_positions_enc.capacity(), self.minimizer_positions_enc.len());
        //  log::info!("Snpmer positions: capacity = {}, length = {}", self.snpmer_positions_enc.capacity(), self.snpmer_positions_enc.len());
        //  log::info!("DNA seq: capacity = {}, length = {}", self.dna_seq.capacity(), self.dna_seq.len());
        //  log::info!("Qual seq: capacity = {}, length = {}", if let Some(qual_seq) = self.qual_seq.as_ref(){qual_seq.capacity()} else {0}, if let Some(qual_seq) = self.qual_seq.as_ref(){qual_seq.len()} else {0});
    }

    pub fn shrink_to_fit(&mut self){
        self.minimizer_positions_enc.shrink_to_fit();
        self.snpmer_positions_enc.shrink_to_fit();
        self.dna_seq.shrink_to_fit();
        if let Some(qual_seq) = self.qual_seq.as_mut(){
            qual_seq.shrink_to_fit();
        }
    }

    #[inline]
    pub fn kmer_from_position(&self, pos:u32, k: usize) -> Kmer48{
        let pos = pos as usize;
        if pos + k > self.dna_seq.len() {
            dbg!(&self.id, pos, k, self.dna_seq.len(), "Position out of bounds for k-mer extraction");
            dbg!(self.minimizer_positions_enc.len(), self.snpmer_positions_enc.len());
            dbg!(self.huffman_encoded);
            panic!("Position out of bounds for k-mer extraction");
        }

        //match various values of k from 17 - 23, odd
        //bio-seq stores CA = 0001.
        //our representation is CA = 0100.
        // Internally, we want CA = 8 HEX = 0100. So 
        let kmer = match k{
            17 => {
                reverse_bit_pairs(Kmer::<Dna, 17, u64>::unsafe_from_seqslice(&self.dna_seq[pos..pos + k]).bs, k)
            }
            19 => {
                reverse_bit_pairs(Kmer::<Dna, 19, u64>::unsafe_from_seqslice(&self.dna_seq[pos..pos + k]).bs, k)
            }
            21 => {
                reverse_bit_pairs(Kmer::<Dna, 21, u64>::unsafe_from_seqslice(&self.dna_seq[pos..pos + k]).bs, k)
            }
            23 => {
                reverse_bit_pairs(Kmer::<Dna, 23, u64>::unsafe_from_seqslice(&self.dna_seq[pos..pos + k]).bs, k)
            }
            _ => {
                panic!("Invalid kmer size")
            }
        };

        // get canonical k-mer based on sides
        let reverse_kmer = reverse_bit_pairs(kmer ^ (u64::MAX), k);
        let mid_mask = !(3 << (k - 1));
        if reverse_kmer & mid_mask < kmer & mid_mask{
            Kmer48::from_u64(reverse_kmer)
        }
        else{
            Kmer48::from_u64(kmer)
        }
    }

    pub fn clear(&mut self){
        self.minimizer_positions_enc.clear();
        self.snpmer_positions_enc.clear();
        self.minimizer_positions_enc.shrink_to_fit();
        self.snpmer_positions_enc.shrink_to_fit();
    }

    /// Decode minimizer positions (dispatches based on per-read flag)
    #[inline]
    pub fn minimizer_positions(&self) -> Vec<u32> {
        decode_positions(&self.minimizer_positions_enc, PositionKind::Minimizer, self.huffman_encoded)
    }

    /// Decode snpmer positions (dispatches based on per-read flag)
    #[inline]
    pub fn snpmer_positions(&self) -> Vec<u32> {
        decode_positions(&self.snpmer_positions_enc, PositionKind::Snpmer, self.huffman_encoded)
    }

    /// Get the count of minimizers (without full decoding)
    pub fn minimizer_count(&self) -> usize {
        count_positions(&self.minimizer_positions_enc, self.huffman_encoded)
    }

    /// Get the count of snpmers (without full decoding)
    pub fn snpmer_count(&self) -> usize {
        count_positions(&self.snpmer_positions_enc, self.huffman_encoded)
    }

    pub fn minimizer_kmers(&self) -> Vec<Kmer48> {
        let positions = self.minimizer_positions();
        positions.iter().map(|&x| self.kmer_from_position(x, self.k as usize)).collect()
    }

    pub fn snpmer_kmers(&self) -> Vec<Kmer48> {
        let positions = self.snpmer_positions();
        positions.iter().map(|&x| self.kmer_from_position(x, self.k as usize)).collect()
    }

    pub fn minimizers_vec(&self) -> Vec<(u32, Kmer48)> {
        let positions = self.minimizer_positions();
        positions.iter().map(|&x| (x, self.kmer_from_position(x, self.k as usize))).collect()
    }

    pub fn snpmers_vec(&self) -> Vec<(u32, Kmer48)> {
        let positions = self.snpmer_positions();
        positions.iter().map(|&x| (x, self.kmer_from_position(x, self.k as usize))).collect()
    }

    // Retain only the minimizers at the given INDICES, not positions.
    // Re-encodes using the SAME encoding format (huffman or varint) as the read currently has.
    pub fn retain_mini_indices(&mut self, indices: FxHashSet<usize>) {
        let positions = self.minimizer_positions();
        let mut filtered_positions = Vec::new();

        for (i, pos) in positions.iter().enumerate() {
            if indices.contains(&i) {
                filtered_positions.push(*pos);
            }
        }

        filtered_positions.shrink_to_fit();
        if self.huffman_encoded {
            let codes = HUFFMAN_CODES_MINI.get().expect("huffman_encoded but coding not initialized");
            self.minimizer_positions_enc = encode_positions_huffman(&filtered_positions, codes);
        } else {
            self.minimizer_positions_enc = encode_positions_delta(&filtered_positions);
        }
        // Do NOT change self.huffman_encoded -- both fields must stay in sync
    }

    // Retain only the snpmers at the given INDICES, not positions.
    // Re-encodes using the SAME encoding format (huffman or varint) as the read currently has.
    pub fn retain_snpmer_indices(&mut self, indices: FxHashSet<usize>) {
        let positions = self.snpmer_positions();
        let mut filtered_positions = Vec::new();

        for (i, pos) in positions.iter().enumerate() {
            if indices.contains(&i) {
                filtered_positions.push(*pos);
            }
        }

        filtered_positions.shrink_to_fit();
        if self.huffman_encoded {
            let codes = HUFFMAN_CODES_SNP.get().expect("huffman_encoded but coding not initialized");
            self.snpmer_positions_enc = encode_positions_huffman(&filtered_positions, codes);
        } else {
            self.snpmer_positions_enc = encode_positions_delta(&filtered_positions);
        }
        // Do NOT change self.huffman_encoded -- both fields must stay in sync
    }


    // Re-encodes both position fields from another read, shifted and filtered.
    // Uses current global encoding state and sets huffman_encoded accordingly.
    pub fn shift_and_retain(&mut self, other_read: &TwinRead, last_break: usize, bp_start: usize, k: usize){
        let other_mini_positions = other_read.minimizer_positions();
        let other_snp_positions = other_read.snpmer_positions();

        let mut mini_positions_filtered = Vec::new();
        let mut snp_positions_filtered = Vec::new();

        for &pos in other_mini_positions.iter(){
            if pos >= last_break as u32 && pos + k as u32 - 1 < bp_start as u32{
                mini_positions_filtered.push(pos - last_break as u32);
            }
        }

        for &pos in other_snp_positions.iter(){
            if pos >= last_break as u32 && pos + k as u32 - 1 < bp_start as u32{
                snp_positions_filtered.push(pos - last_break as u32);
            }
        }

        mini_positions_filtered.shrink_to_fit();
        snp_positions_filtered.shrink_to_fit();

        // shift_and_retain creates a new encoding from scratch for both fields
        // simultaneously, so it's safe to use auto-dispatch here
        self.minimizer_positions_enc = encode_positions(&mini_positions_filtered, PositionKind::Minimizer);
        self.snpmer_positions_enc = encode_positions(&snp_positions_filtered, PositionKind::Snpmer);
        self.huffman_encoded = huffman_initialized();
    }
}

pub fn retain_vec_indices<T>(vec: &mut Vec<T>, positions: &FxHashSet<usize>){
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
    pub solid_kmers: HashSet<Kmer48>,
    pub use_solid_kmers: bool,
    pub high_freq_thresh: f64,
    pub high_freq_kmers: HashSet<Kmer48>,
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
    pub mid_bases: [u8;2],
    pub counts: [u32;2],
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
    pub chain_score: i32,
    pub large_indel: bool,
    pub minimizer_chain: Option<Vec<Anchor>>,
}

impl TwinOverlap{
    pub fn length1(&self) -> usize{
        self.end1 - self.start1
    }

    pub fn length2(&self) -> usize{
        self.end2 - self.start2
    }
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
pub struct AnchorBuilder{
    pub pos1: u32,
    pub pos2: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Eq, Hash)]
pub struct Anchor{
    pub i: Option<u32>,
    pub j: Option<u32>,
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
    pub max_alignment_boundaries: Option<Lapper<u32, SmallTwinOl>>,
    //pub max_mapping_boundaries: Option<Lapper<u32, BareMappingOverlap>>,
    pub max_mapping_boundaries: Option<Vec<(BareInterval, BareMappingOverlap)>>,
    //pub kmer_counts: Vec<u32>,
    pub present: bool,
    pub length: usize,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BareMappingOverlap{
    pub snpmer_identity: Fraction,
}

impl Eq for BareMappingOverlap{}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct TwoCycle {
    pub read_i: usize,
    pub read_j: usize,
    pub hang_penalty: i64, 
    pub circular_length: usize,
    pub total_mini: usize,
}


#[derive(Debug, Clone,  Default)]
pub struct SmallTwinOl{
    pub query_id: u32,
    pub snpmer_identity: f32,
    pub reverse: bool,
    pub alignment_result: Option<AlignmentResult>
}

#[derive(Debug, Clone, PartialEq, Default, Eq)]
pub struct BareInterval{
    pub start: u32,
    pub stop: u32
}

impl Ord for BareInterval
{
    #[inline]
    fn cmp(&self, other: &BareInterval) -> Ordering {
        match self.start.cmp(&other.start) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => self.stop.cmp(&other.stop),
        }
    }
}

impl PartialOrd for BareInterval
{
    #[inline]
    fn partial_cmp(&self, other: &BareInterval) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for SmallTwinOl{}

impl PartialEq for SmallTwinOl{
    fn eq(&self, other: &Self) -> bool{
        self.query_id == other.query_id && self.snpmer_identity == other.snpmer_identity && self.reverse == other.reverse 
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Breakpoints {
    pub pos1: usize,
    pub pos2: usize,
    pub cov: usize,
    pub condition: i64,
}

/// Coverage statistics for a read segment
#[derive(Debug, Clone)]
pub struct CoverageStats {
    pub start: usize,
    pub end: usize,
    pub min_depth_multi: MultiCov,
    pub median_depth: f64,
    pub snpmer_id_threshold: f64,
}

/// Compact plan for splitting a read - stores only breakpoints and coverage stats
#[derive(Debug, Clone)]
pub struct SplitReadPlan {
    pub read_index: usize,
    pub breakpoints: Vec<Breakpoints>,
    pub coverage_stats: Vec<CoverageStats>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct GetSequenceInfoConfig{
    pub blunted: bool,
    pub dna_seq_info: bool,
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

pub struct BeamStartState{
    pub initial_unitig_length: usize,
    pub initial_unitig_size: usize,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct OverlapAdjMap {
    pub adj_map: FxHashMap<NodeIndex, Vec<NodeIndex>>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct OpLenVec{
    op_vec: NibbleVec<[u8; 32]>,
    // Variable-length encoded lengths
    len_vec: Vec<u8>,
    // Number of operations (needed since varints have variable size)
    count: usize,
}

impl OpLenVec{
    pub fn new(cigar_vec: Vec<OpLen>) -> Self{
        let mut op_vec = Vec::new();
        let mut lengths: Vec<u32> = Vec::new();
        let mut last_op = Operation::Sentinel;
        for op_len in cigar_vec{
            if op_len.op == last_op{
                *lengths.last_mut().unwrap() += op_len.len as u32;
            }
            else if (op_len.op == Operation::X || op_len.op == Operation::Eq) && last_op == Operation::M{
                *lengths.last_mut().unwrap() += op_len.len as u32;
            }
            else{
                op_vec.push(op_len.op as u8);
                lengths.push(op_len.len as u32);
                last_op = op_len.op;
            }
            if last_op == Operation::X || last_op == Operation::Eq{
                last_op = Operation::M;
            }
        }
        assert!(op_vec.len() == lengths.len());
        let count = op_vec.len();
        op_vec.shrink_to_fit();

        // Encode lengths as variable-length integers
        let mut len_vec = encode_varints(&lengths);
        len_vec.shrink_to_fit();

        OpLenVec{
            op_vec: NibbleVec::<[u8; 32]>::from_byte_vec(op_vec),
            len_vec,
            count,
        }
    }

    pub fn len(&self) -> usize{
        self.count
    }

    pub fn iter(&self) -> OpLenIter<'_> {
        OpLenIter {
            op_bytes: self.op_vec.as_bytes(),
            len_bytes: &self.len_vec,
            op_idx: 0,
            len_pos: 0,
            count: self.count,
        }
    }
}

/// Iterator over (Operation, u32) pairs from an OpLenVec
pub struct OpLenIter<'a> {
    op_bytes: &'a [u8],
    len_bytes: &'a [u8],
    op_idx: usize,
    len_pos: usize,
    count: usize,
}

impl<'a> Iterator for OpLenIter<'a> {
    type Item = (Operation, u32);

    fn next(&mut self) -> Option<Self::Item> {
        if self.op_idx >= self.count {
            return None;
        }

        let op = u8_to_operation(self.op_bytes[self.op_idx]);
        let (len, bytes_consumed) = decode_varint_u32(self.len_bytes, self.len_pos);

        self.op_idx += 1;
        self.len_pos += bytes_consumed;

        Some((op, len))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.count - self.op_idx;
        (remaining, Some(remaining))
    }
}

impl<'a> ExactSizeIterator for OpLenIter<'a> {}

pub fn u8_to_operation(b: u8) -> Operation{
    match b{
        0 => Operation::Sentinel,
        1 => Operation::M,
        2 => Operation::Eq,
        3 => Operation::X,
        4 => Operation::I,
        5 => Operation::D,
        _ => Operation::Sentinel,
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

#[derive(Debug, Clone, PartialEq)]
pub struct CompareTwinReadOptions{
    pub compare_snpmers: bool,
    pub retain_chain: bool,
    pub force_query_nonoverlap: bool,
    pub force_ref_nonoverlap: bool,
    pub supplementary_threshold_score: Option<f64>,
    pub supplementary_threshold_ratio: Option<f64>, 
    // When not forcing 1-to-1 alignments, allow query overlaps only if secondary threshold is below a certain amount
    pub secondary_threshold: Option<f64>,
    //Preload
    pub read1_mininimizers: Option<Vec<(u32,Kmer48)>>,
    pub read1_snpmers: Option<Vec<(u32,Kmer48)>>,
    pub max_gap: usize,
    pub double_gap: usize,
    pub maximal_only: bool,
}

impl Default for CompareTwinReadOptions{
    fn default() -> Self{
        CompareTwinReadOptions{
            compare_snpmers: true,
            retain_chain: false,
            force_query_nonoverlap: false,
            force_ref_nonoverlap: true,
            supplementary_threshold_score: Some(500.0),
            supplementary_threshold_ratio: Some(0.25),
            secondary_threshold: Some(0.50),
            read1_mininimizers: None,
            read1_snpmers: None,
            max_gap: MAX_GAP_CHAINING,
            double_gap: 10_000,
            maximal_only: false,
        }
    }
}

pub struct HeavyCutOptions<'a> 
{
    pub samples: usize,
    pub temperature: f64,
    pub steps: usize,
    pub max_forward: usize,
    pub max_reads_forward: usize,
    pub safe_length_back: usize,
    pub ol_thresh: f64,
    pub tip_threshold: usize,
    pub strain_repeat_map: Option<&'a FxHashMap<NodeIndex, FxHashSet<NodeIndex>>>,
    pub special_small: bool,
    pub max_length_search: usize,
    pub require_safety: bool,
    pub only_tips: bool,
    pub cut_tips: bool,
    pub debug: bool,
}




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn convert_from_u64() {
        let kmer = Kmer48::from_u64(1);
        assert_eq!(kmer.data, [1, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn bioseq_vs_ours(){
        let kmer_bioseq = Kmer::<Dna, 5, u64>::unsafe_from_seqslice(dna!("ACGTG"));
        let kmer_bioseq_u64 = kmer_bioseq.bs;
        println!("{:08b}", kmer_bioseq_u64);
        let kmer_ours = 0b00_01_10_11_10;
        let reversed_bioseq = reverse_bit_pairs(kmer_bioseq_u64, 5);
        println!("{:08b}", reversed_bioseq);

        assert_eq!(reverse_bit_pairs(kmer_bioseq_u64, 5), kmer_ours);
    }

    #[test]
    fn reverse_comp_kmer(){
        // ACGTG
        let kmer_ours = 0b00_01_10_11_10;
        //reverse comp is CACGT
        let goal = 0b01_00_01_10_11;
        let kmer_ours_rev_comp = reverse_bit_pairs(kmer_ours ^ (u64::MAX), 5);

        assert_eq!(kmer_ours_rev_comp, goal);
    }

    #[test]
    fn test_delta_encoding_small_deltas() {
        // Test with small deltas (< 128) - typical minimizer spacing of ~10 bases
        let positions = vec![0, 10, 20, 35, 45, 55];
        let encoded = encode_positions_delta(&positions);
        let decoded = decode_positions_delta(&encoded);

        assert_eq!(positions, decoded);
        // Each delta should be 1 byte for small values
        assert_eq!(encoded.len(), positions.len());
    }

    #[test]
    fn test_delta_encoding_large_deltas() {
        // Test with large deltas - typical SNPmer spacing of ~1000 bases
        let positions = vec![0, 1000, 2100, 3050];
        let encoded = encode_positions_delta(&positions);
        let decoded = decode_positions_delta(&encoded);

        assert_eq!(positions, decoded);
        // 1000 needs 2 bytes, 1100 needs 2 bytes, 950 needs 2 bytes
        assert!(encoded.len() < positions.len() * 4); // Still better than u32
    }

    #[test]
    fn test_delta_encoding_mixed() {
        // Test with mixed small and large deltas
        let positions = vec![0, 5, 15, 1015, 1025, 3000];
        let encoded = encode_positions_delta(&positions);
        let decoded = decode_positions_delta(&encoded);

        assert_eq!(positions, decoded);
    }

    #[test]
    fn test_delta_encoding_very_large() {
        // Test with very large deltas that need full u32 range
        let positions = vec![0, 100_000_000, 200_000_000];
        let encoded = encode_positions_delta(&positions);
        let decoded = decode_positions_delta(&encoded);

        assert_eq!(positions, decoded);
    }

    #[test]
    fn test_delta_encoding_empty() {
        let positions: Vec<u32> = vec![];
        let encoded = encode_positions_delta(&positions);
        let decoded = decode_positions_delta(&encoded);

        assert_eq!(positions, decoded);
        assert_eq!(encoded.len(), 0);
    }

    #[test]
    fn test_delta_encoding_single() {
        let positions = vec![42];
        let encoded = encode_positions_delta(&positions);
        let decoded = decode_positions_delta(&encoded);

        assert_eq!(positions, decoded);
        assert_eq!(encoded.len(), 1);
    }

    // ---- Huffman round-trip tests ----
    // These build local Coding instances to avoid polluting global OnceLock state.

    use minimum_redundancy::{Coding, Code, BitsPerFragment, Frequencies};
    use std::collections::HashMap;

    /// Build a Huffman coding from sample positions, seeding all 256 byte values.
    fn build_test_coding(sample_positions: &[Vec<u32>]) -> (Coding<u8, BitsPerFragment>, [Code; 256]) {
        let mut freq = HashMap::<u8, usize>::new();
        for b in 0..=255u8 { freq.insert(b, 1); }
        for positions in sample_positions {
            let varint_bytes = encode_positions_delta(positions);
            freq.add_occurences_of(varint_bytes.into_iter());
        }
        let coding = Coding::from_frequencies(BitsPerFragment(1), freq);
        let codes = coding.codes_for_values_array();
        (coding, codes)
    }

    #[test]
    fn test_huffman_roundtrip_small_deltas() {
        let positions = vec![0, 10, 20, 35, 45, 55];
        let (coding, codes) = build_test_coding(&[positions.clone()]);

        let encoded = encode_positions_huffman(&positions, &codes);
        let decoded = decode_positions_huffman(&encoded, &coding);
        assert_eq!(positions, decoded);
    }

    #[test]
    fn test_huffman_roundtrip_large_deltas() {
        let positions = vec![0, 1000, 2100, 3050, 10000];
        let (coding, codes) = build_test_coding(&[positions.clone()]);

        let encoded = encode_positions_huffman(&positions, &codes);
        let decoded = decode_positions_huffman(&encoded, &coding);
        assert_eq!(positions, decoded);
    }

    #[test]
    fn test_huffman_roundtrip_very_large_deltas() {
        // Deltas that produce rare varint bytes (multi-byte varints)
        let positions = vec![0, 100_000, 200_000, 300_000_000];
        let (coding, codes) = build_test_coding(&[positions.clone()]);

        let encoded = encode_positions_huffman(&positions, &codes);
        let decoded = decode_positions_huffman(&encoded, &coding);
        assert_eq!(positions, decoded);
    }

    #[test]
    fn test_huffman_roundtrip_empty() {
        let positions: Vec<u32> = vec![];
        let (coding, codes) = build_test_coding(&[positions.clone()]);

        let encoded = encode_positions_huffman(&positions, &codes);
        assert!(encoded.is_empty());
        let decoded = decode_positions_huffman(&encoded, &coding);
        assert_eq!(positions, decoded);
    }

    #[test]
    fn test_huffman_roundtrip_single() {
        let positions = vec![42];
        let (coding, codes) = build_test_coding(&[positions.clone()]);

        let encoded = encode_positions_huffman(&positions, &codes);
        let decoded = decode_positions_huffman(&encoded, &coding);
        assert_eq!(positions, decoded);
    }

    #[test]
    fn test_huffman_count() {
        let positions = vec![0, 10, 20, 35, 45, 55, 100, 200];
        let (_, codes) = build_test_coding(&[positions.clone()]);

        let encoded = encode_positions_huffman(&positions, &codes);
        let count = count_positions_huffman(&encoded);
        assert_eq!(count, positions.len());
    }

    #[test]
    fn test_huffman_unseen_bytes() {
        // Build coding from only small deltas (bytes 0-127 only),
        // then encode positions with large deltas producing bytes 128+ (continuation bytes).
        // This is the exact scenario that caused Bug 1.
        let small_positions = vec![0, 5, 10, 15, 20, 25, 30];
        let (coding, codes) = build_test_coding(&[small_positions]);

        // Now encode positions with large deltas that produce continuation bytes (0x80+)
        let large_positions = vec![0, 1000, 5000, 50000];
        let encoded = encode_positions_huffman(&large_positions, &codes);
        let decoded = decode_positions_huffman(&encoded, &coding);
        assert_eq!(large_positions, decoded, "Unseen byte values should still round-trip correctly");
    }

    #[test]
    fn test_huffman_all_256_codes_valid() {
        // Verify every byte value gets a non-zero-length code
        let positions = vec![0, 10, 20];
        let (_, codes) = build_test_coding(&[positions]);
        for b in 0..=255u8 {
            assert!(codes[b as usize].len > 0,
                "Byte {} has zero-length code -- would corrupt bitstream", b);
        }
    }

    #[test]
    fn test_huffman_roundtrip_many_positions() {
        // Stress test with many positions covering a wide range of delta sizes
        let mut positions = Vec::new();
        let mut pos = 0u32;
        for i in 0..500 {
            positions.push(pos);
            // Mix of small and large deltas
            pos += match i % 5 {
                0 => 3,
                1 => 15,
                2 => 200,
                3 => 5000,
                _ => 1,
            };
        }
        let (coding, codes) = build_test_coding(&[positions.clone()]);
        let encoded = encode_positions_huffman(&positions, &codes);
        let decoded = decode_positions_huffman(&encoded, &coding);
        assert_eq!(positions, decoded);
        assert_eq!(count_positions_huffman(&encoded), positions.len());
    }

    #[test]
    fn test_huffman_coding_trained_on_different_data() {
        // Build coding from one distribution, use it to encode a very different distribution.
        // This tests that all byte values work even when the coding was trained on different data.
        let training = vec![vec![0, 1, 2, 3, 4, 5]]; // tiny deltas only
        let (coding, codes) = build_test_coding(&training);

        // Encode with huge deltas (multi-byte varints with continuation bits)
        let test_positions = vec![0, 100_000_000, 200_000_000, 300_000_000];
        let encoded = encode_positions_huffman(&test_positions, &codes);
        let decoded = decode_positions_huffman(&encoded, &coding);
        assert_eq!(test_positions, decoded);
    }
}
