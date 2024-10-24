use fxhash::FxHashMap;
use crate::types::*;
use disjoint::DisjointSet;
use fxhash::FxHashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use crate::mapping::*;
use serde::{Deserialize, Serialize};


#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Hash, Eq)]
pub struct ReadData {
    pub index: NodeIndex,
    pub in_edges: Vec<EdgeIndex>,
    pub out_edges: Vec<EdgeIndex>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Hash, Eq)]
pub struct ReadOverlapEdgeTwin {
    pub node1: NodeIndex,
    pub node2: NodeIndex,
    pub forward1: bool,
    pub forward2: bool,
    pub overlap_len_bases: usize,
    pub shared_minimizers: usize,
    pub diff_snpmers: usize
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct OverlapTwinGraph {
    pub reads: Vec<TwinRead>,
    pub nodes: FxHashMap<NodeIndex, ReadData>,
    pub edges: Vec<Option<ReadOverlapEdgeTwin>>,
}


// Transitive reduction implementation
impl OverlapTwinGraph{
    pub fn other_node(&self, node: NodeIndex, edge: &ReadOverlapEdgeTwin) -> NodeIndex {
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
        let mut mark: FxHashMap<NodeIndex, Mark> = FxHashMap::default();

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
                        let lensum = if length + length2 < self.reads[other_node].base_length {0} else {length + length2 - self.reads[other_node].base_length};
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
            let string_length = self.reads[n1.index].base_length
                + self.reads[n2.index].base_length
                - edge.overlap_len_bases;
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



pub fn read_graph_from_overlaps_twin(all_reads_cat: Vec<TwinRead>, overlaps: &Vec<TwinOverlap>, c: usize) -> OverlapTwinGraph
{
    let mut nodes = FxHashMap::default();
    let mut edges = vec![];

    for twlap in overlaps.iter() {
        let i = twlap.i1;
        let j = twlap.i2;
        let read1 = &all_reads_cat[i];
        let read2 = &all_reads_cat[j];

        //check if end-to-end overlap
        let mut forward1 = false;
        let mut forward2 = false;
        let mut r1_r2 = true;
        let mut ol = false;
        let identity = id_est(twlap.shared_minimizers, twlap.diff_snpmers, c as u64);
        let same_strain = if identity > 0.9999 {true} else {false};
        if !same_strain {
            continue;
        }

        let overlap_hang_length = 750;

        if twlap.chain_reverse {
            if twlap.start1 < overlap_hang_length && twlap.start2 + overlap_hang_length < read2.base_length{
                forward1 = false;
                forward2 = true;
                r1_r2 = true;
                ol = true;
            } 
            else if twlap.end1 + overlap_hang_length > read1.base_length && twlap.end2 + overlap_hang_length > read2.base_length 
            {
                forward1 = true;
                forward2 = false;
                r1_r2 = true;
                ol = true;
            }
        } else {
            if twlap.start1 < overlap_hang_length && twlap.end2 + overlap_hang_length > read2.base_length {
                forward1 = true;
                forward2 = true;
                r1_r2 = false;
                ol = true;
            } 
            else if twlap.end1 + overlap_hang_length > read1.base_length  && twlap.start2 < overlap_hang_length {
                forward1 = true;
                forward2 = true;
                r1_r2 = true;
                ol = true;
            }
        }

        let aln_len1 = twlap.end1 - twlap.start1 + 1;
        let aln_len2 = twlap.end2 - twlap.start2 + 1;
        if ol && aln_len1.max(aln_len2) > 1500 {
            log::info!(
                "OVERLAP {} {} {} {} {}-{} {} {}-{}, REVERSE: {}",
                i,
                j,
                identity * 100.,
                read1.base_length,
                twlap.start1,
                twlap.end1,
                read2.base_length,
                twlap.start2,
                twlap.end2,
                twlap.chain_reverse
            );

            let start = if r1_r2 { i } else { j };
            let end = if r1_r2 { j } else { i };
            let new_read_overlap = ReadOverlapEdgeTwin {
                node1: start,
                node2: end,
                forward1,
                forward2,
                overlap_len_bases: aln_len1.max(aln_len2),
                shared_minimizers: twlap.shared_minimizers,
                diff_snpmers: twlap.diff_snpmers
            };

            edges.push(Some(new_read_overlap));
            let ind = edges.len() - 1;
            {
                let rd1 = nodes.entry(i).or_insert(ReadData::default());
                rd1.index = i;
                if r1_r2 && forward1 {
                    rd1.out_edges.push(ind)
                } else if !r1_r2 && !forward2 {
                    rd1.out_edges.push(ind)
                } else {
                    rd1.in_edges.push(ind)
                }
            }
            {
                let rd2 = nodes.entry(j).or_insert(ReadData::default());
                rd2.index = j;
                if r1_r2 && forward2 {
                    rd2.in_edges.push(ind)
                } else if !r1_r2 && !forward1 {
                    rd2.in_edges.push(ind)
                } else {
                    rd2.out_edges.push(ind)
                }
            }
        }
    }

    let graph = OverlapTwinGraph { reads: all_reads_cat, nodes, edges };

    return graph;
}


pub fn print_graph_stdout(graph: &OverlapTwinGraph, file: &str) {
    let mut bufwriter = BufWriter::new(File::create(file).unwrap());
    let all_reads_cat = &graph.reads;
    for edge in graph.edges.iter() {
        if let Some(edge) = edge {
            let i = edge.node1;
            let j = edge.node2;
            let read = &all_reads_cat[i];
            let read2 = &all_reads_cat[j];
            let forward1 = edge.forward1;
            let forward2 = edge.forward2;
            let aln_len = edge.overlap_len_bases;

            if read.id.contains("junk") || read2.id.contains("junk") || read.id.contains("chimera") || read2.id.contains("chimera"){
                continue;
            }

            let res1 = parse_badread(&read.id);
            let res2 = parse_badread(&read2.id);
            if res1.is_none() || res2.is_none() {
                continue;
            }
            let (name1, range1) = res1.unwrap();
            let (name2, range2) = res2.unwrap();
            
            writeln!(bufwriter, 
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i, j, name1, name2, range1, range2, forward1, forward2
            ).unwrap();
        }
    }
}


pub fn get_overlaps_outer_reads_twin(twin_reads: &[TwinRead], outer_read_indices: &[usize],  k: u64, c: u64) -> Vec<TwinOverlap>{

    let mut bufwriter = BufWriter::new(File::create("overlaps.txt").unwrap());
    let vec_format = twin_reads.iter().enumerate().filter(|(i, _)| outer_read_indices.contains(i)).map(|(i, x)| (i,&x.minimizers)).collect::<Vec<_>>();

    let mut inverted_index_hashmap_outer = FxHashMap::default();
    for (i, x) in vec_format.iter(){
        for &y in x.iter(){
            inverted_index_hashmap_outer.entry(y.1).or_insert(vec![]).push(i);
        }
    }
    
    let mut overlaps = vec![];

    let mut compared_reads = FxHashSet::default();
    for i in outer_read_indices.iter() {
        let read = &twin_reads[*i];
        let mut index_count_map = FxHashMap::default();
        for &y in read.minimizers.iter() {
            if let Some(indices) = inverted_index_hashmap_outer.get(&y.1) {
                for &index in indices {
                    if index == i {
                        continue;
                    }
                    *index_count_map.entry(index).or_insert(0) += 1;
                }
            }
        }
        //sort
        let mut top_indices = index_count_map.into_iter().collect::<Vec<_>>();
        top_indices.sort_by(|a, b| b.1.cmp(&a.1));

        for (index, count) in top_indices.iter() {
            let sorted_readpair = if i < *index { (i, *index) } else { (*index, i) };
            if compared_reads.contains(&sorted_readpair) {
                continue;
            }
            let read2 = &twin_reads[**index];
            let twlap = compare_twin_reads(&read, &read2, *i, **index);
            if twlap.is_none(){
                continue;
            }
            let twlap = twlap.unwrap();
            let identity = id_est(twlap.shared_minimizers, twlap.diff_snpmers, c);
            writeln!(bufwriter,
                "intersect tigs: {} i: {} j: {} fsv: {} leni: {} {}-{} lenj: {} {}-{} shared_mini: {} REVERSE:{}, snp_diff: {} snp_share: {}, intersect: {:?}",
                count,
                i,
                *index,
                identity * 100.,
                read.base_length,
                twlap.start1,
                twlap.end1,
                read2.base_length,
                twlap.start2,
                twlap.end2,
                twlap.shared_minimizers,
                twlap.chain_reverse,
                twlap.diff_snpmers,
                twlap.shared_snpmers,
                twlap.intersect
            ).unwrap();
            if identity > 0.9999 {
                overlaps.push(twlap);
            }
            compared_reads.insert(sorted_readpair);
        }
    }

    return overlaps;

}

pub fn remove_contained_reads_twin<'a>(twin_reads: &'a [TwinRead], k: u64, c: u64) -> Vec<usize>{
    let inverted_index_hashmap =
        twin_reads
            .iter()
            .enumerate()
            .fold(FxHashMap::default(), |mut acc, (i, x)| {
                for &y in x.minimizers.iter() {
                    acc.entry(y.1).or_insert(vec![]).push(i);
                }
                acc
            });

    //open file for writing
    let mut bufwriter = BufWriter::new(File::create("histo.txt").unwrap());
    let mut bufwriter_dbg = BufWriter::new(File::create("twin_dbg.txt").unwrap());
    let mut contained_reads = FxHashSet::default();
    let mut outer_reads = vec![];

    for i in 0..twin_reads.len() {
        let mut contained = false;
        let read1 = &twin_reads[i];
        let mut index_count_map = FxHashMap::default();
        for &y in read1.minimizers.iter() {
            if let Some(indices) = inverted_index_hashmap.get(&y.1) {
                for &index in indices {
                    if index == i {
                        continue;
                    }
                    *index_count_map.entry(index).or_insert(0) += 1;
                }
            }
        }

        //look at top 5 indices and do disjoint_set distance
        let mut top_indices = index_count_map.iter().collect::<Vec<_>>();
        top_indices.sort_by(|a, b| b.1.cmp(a.1));
        for (index, count) in top_indices.iter() {
            if **count < read1.minimizers.len() as i32 / 20 {
                break;
            }
            if contained_reads.contains(*index) {
                continue;
            }

            let read2 = &twin_reads[**index];
            let twin_overlap = compare_twin_reads(&read1, &read2, i, **index);
            let res1 = parse_badread(&read1.id);
            let res2 = parse_badread(&read2.id);
            if res1.is_none() || res2.is_none() {
                continue;
            }

            let (name1, _range1) = res1.unwrap();
            let (name2, _range2) = res2.unwrap();
            let within = if { name1 == name2 } { true } else { false };

            if twin_overlap.is_none() {
                continue;
            }

            let twin_overlap = twin_overlap.unwrap();
            let shared_minimizers = twin_overlap.shared_minimizers;
            let diff_snpmers = twin_overlap.diff_snpmers;

            let mut same_strain = false;
            let identity;
            if diff_snpmers == 0 {
                same_strain = true;
                identity = 1.;
            } else {
                identity = id_est(shared_minimizers, diff_snpmers, c);
                if identity > 0.9999 {
                    same_strain = true;
                }
            }

            write!(
                bufwriter,
                "{}\t{}\t{}\t{}\n",
                twin_overlap.shared_snpmers, identity, within, same_strain
            )
            .unwrap();
            let len1 = twin_overlap.end1 - twin_overlap.start1;
            let len2 = twin_overlap.end2 - twin_overlap.start2;
            let ol_len = len1.max(len2);

            writeln!(
                bufwriter_dbg,
                "{} ----- {}   minis: {} shared_snps: {}, diff_snps: {}, identity {}, ol_len {}, read1len: {}, read2len: {}",
                &i,
                &index,
                twin_overlap.shared_minimizers,
                twin_overlap.shared_snpmers,
                twin_overlap.diff_snpmers,
                identity,
                ol_len,
                read1.base_length,
                read2.base_length

            )
            .unwrap();


            if ol_len as f64 > 0.9 * (read1.base_length as f64) && same_strain {
                log::trace!("{} {} CONTAINED", read1.id, read2.id);
                contained = true;
                contained_reads.insert(i);
                break;
            }
        }

        if !contained {
            outer_reads.push(i);
        }
    }
    return outer_reads;
}
