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

use fxhash::FxHashMap;
use fxhash::FxHashSet;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::{BuildHasherDefault, Hasher};

pub type Kmer64 = u64;
pub type Kmer32 = u32;
pub type KmerHash64 = u64;
pub type KmerHash32 = u32;

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
    pub chain: Vec<(usize, usize, usize)>,
    pub reverse: bool,
    pub score: f64,
}

pub type EdgeIndex = usize;
pub type NodeIndex = usize;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Hash, Eq)]
pub struct ReadData {
    pub index: NodeIndex,
    pub in_edges: Vec<EdgeIndex>,
    pub out_edges: Vec<EdgeIndex>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Hash, Eq)]
pub struct ReadOverlapEdge {
    pub node1: NodeIndex,
    pub node2: NodeIndex,
    pub forward1: bool,
    pub forward2: bool,
    pub overlap_len_bases: usize,
    pub overlap_len_tigs: usize,
    pub shared_tigs: usize,
    pub variable_tigs: usize,
    pub variable_roots: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct OverlapGraph {
    pub reads: Vec<TigRead>,
    pub nodes: FxHashMap<NodeIndex, ReadData>,
    pub edges: Vec<Option<ReadOverlapEdge>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct TigRead {
    pub tig_seq: Vec<u32>,
    pub id: String,
}

// Transitive reduction implementation
impl OverlapGraph {
    pub fn other_node(&self, node: NodeIndex, edge: &ReadOverlapEdge) -> NodeIndex {
        if edge.node1 == node {
            edge.node2
        } else {
            edge.node1
        }
    }
    pub fn transitive_reduction(&mut self) {
        const FUZZ: usize = 10;

        // Initialize reduce array for edges
        let mut reduce = vec![false; self.edges.len()];
        let mut mark: HashMap<NodeIndex, Mark> = HashMap::new();

        // Step 1: Mark all nodes as vacant initially and set reduce to false for all edges
        for (node_id, node_data) in self.nodes.iter() {
            mark.insert(*node_id, Mark::Vacant);
            for &edge_id in &node_data.out_edges {
                reduce[edge_id] = false;
            }
        }

        // Step 2: Iterate over all edges, marking nodes and evaluating paths
        for (node_id, node_data) in self.nodes.iter() {
            let sorted_all_edges = self.get_edges_sorted_by_length(*node_id);

            for &(length, edge_id, outer) in sorted_all_edges.iter() {
                let edge = self.edges[edge_id].as_ref().unwrap();
                let other_node = self.other_node(*node_id, edge);
                mark.insert(other_node, Mark::InPlay);
            }

            let longest = sorted_all_edges.last().unwrap().0 + FUZZ;

            for &(length, edge_id, outer) in sorted_all_edges.iter() {
                let edge = self.edges[edge_id].as_ref().unwrap();
                let other_node = self.other_node(*node_id, edge);
                let outgoing = self.outgoing_edge(other_node, edge_id);
                for &(length2, edge_id2, outer2) in &self.get_edges_sorted_by_length(other_node) {
                    //outgoing edge for other node, requires incoming
                    if outgoing && outer2{
                        continue
                    }
                    //incoming edge for other node, requires outgoing
                    if !outgoing && !outer2{
                        continue
                    }
                    if edge_id2 == edge_id{
                        continue
                    }
                    let next_edge = self.edges[edge_id2].as_ref().unwrap();
                    let third_node = self.other_node(other_node, next_edge);
                    if mark.get(&third_node).unwrap() == &Mark::InPlay {
                        let lensum = if length + length2 < self.reads[other_node].tig_seq.len() {0} else {length + length2 - self.reads[other_node].tig_seq.len()};
                        if lensum <= longest + FUZZ || true
                        {
                            if let Some(true) =
                                mark.get(&third_node).map(|m| *m == Mark::InPlay)
                            {
                                mark.insert(third_node, Mark::Eliminated);
                                println!("Eliminated edge from {} to {}, length1 {}, length2 {}, longets {}, lensum {}, edge_info1 {:?}, edge_info2 {:?}, outgoing_edge {}, direction_out {}", node_id, third_node, length, length2, longest, lensum, &edge, &next_edge, outer2, outgoing);
                            }
                        }
                    }
                }
            }
            // Step 3: Final pass to mark reduced edges
            for &edge_id in node_data.out_edges.iter().chain(node_data.in_edges.iter()){
                let edge = self.edges[edge_id].as_ref().unwrap();
                let other_node = self.other_node(*node_id, edge);
                if let Some(Mark::Eliminated) = mark.get(&other_node) {
                    reduce[edge_id] = true;
                }
            }
        }

        // Apply reduction by removing edges
        for (i, &reduced) in reduce.iter().enumerate() {
            if reduced {
                self.edges[i] = None; // Remove reduced edges
            }
        }

        for (_, node_data) in self.nodes.iter_mut() {
            node_data.out_edges.retain(|&edge_id| self.edges[edge_id].is_some());
            node_data.in_edges.retain(|&edge_id| self.edges[edge_id].is_some());
        }
    }

    fn get_edges_sorted_by_length(&self, node_id: NodeIndex) -> Vec<(usize, EdgeIndex, bool)> {
        let mut sorted_all_edges = vec![];
        let node_data = self.nodes.get(&node_id).unwrap();
        for (l, edge_ind) in node_data
            .out_edges
            .iter()
            .chain(node_data.in_edges.iter())
            .enumerate()
        {
            let outer;
            if l < node_data.out_edges.len() {
                outer = true;
            } else {
                outer = false;
            }
            let edge = self.edges[*edge_ind].as_ref().unwrap();
            let n1 = self.nodes.get(&edge.node1).unwrap();
            let n2 = self.nodes.get(&edge.node2).unwrap();
            let string_length = self.reads[n1.index].tig_seq.len()
                + self.reads[n2.index].tig_seq.len()
                - edge.shared_tigs;
            sorted_all_edges.push((string_length, *edge_ind, outer));
        }
        sorted_all_edges.sort();
        sorted_all_edges
    }

    fn outgoing_edge(&self, node_id: NodeIndex, edge_id: EdgeIndex) -> bool {
        let edge = self.edges[edge_id].as_ref().unwrap();
        if edge.node1 == node_id {
            if edge.forward1{
                true
            }
            else{
                false
            }
        } else {
            if edge.forward2{
                false
            }
            else{
                true
            }
        }
    }
}

// Enum for marking the state of a node during processing
#[derive(PartialEq, Eq, Clone)]
enum Mark {
    Vacant,
    InPlay,
    Eliminated,
}
