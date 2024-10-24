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


pub fn remove_contained_reads(all_reads_cat: &Vec<TigRead>, disjoint_set: &DisjointSet, used_ids: &FxHashSet<u64>) -> Vec<usize>
{
    let inverted_index_hashmap =
        all_reads_cat
            .iter()
            .enumerate()
            .fold(FxHashMap::default(), |mut acc, (i, x)| {
                for &y in x.tig_seq.iter() {
                    acc.entry(y).or_insert(vec![]).push(i);
                }
                acc
            });

    //open file for writing
    let mut bufwriter = BufWriter::new(File::create("histo.txt").unwrap());
    let mut bufwriter2 = BufWriter::new(File::create("histo.dbg").unwrap());

    //check contained reads
    let mut outer_reads = vec![];
    let mut contained_reads = FxHashSet::default();
    for (i, read) in all_reads_cat.iter().enumerate() {
        let mut contained = false;
        let mut index_count_map = FxHashMap::default();
        for &y in read.tig_seq.iter() {
            if let Some(indices) = inverted_index_hashmap.get(&y) {
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
            if **count < 5 {
                break;
            }
            if contained_reads.contains(*index) {
                continue;
            }
            let read2 = &all_reads_cat[**index];
            dbg!("{} {} {} {}, {}", i, **index, &read.id, &read2.id);
            let tiglap = disjoint_distance(&read.tig_seq, &read2.tig_seq, &disjoint_set, &used_ids, i, **index);
            if read.id.contains("junk") || read2.id.contains("junk") {
                continue;
            }
            let name1 = read.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[0];
            let name2 = read2.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[0];
            let _range1 = read.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[2];
            let _range2 = read2.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[2];
            let aln_len1 = tiglap.tig1_end - tiglap.tig1_start + 1;
            let aln_len2 = tiglap.tig2_end - tiglap.tig2_start + 1;
            let frac_shared = tiglap.shared_tig as f64 / aln_len1.max(aln_len2) as f64;
            let frac_shared_variable = tiglap.variable_tigs as f64 / tiglap.variable_roots as f64;

            let mut within = false;
            if name1 == name2 {
                within = true;
            }

            if !frac_shared_variable.is_nan()
                && aln_len1 as f64 > (read.tig_seq.len() as f64) * 0.4
                && aln_len2 as f64 > (read2.tig_seq.len() as f64) * 0.4
            {
                write!(
                    bufwriter2,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    tiglap.variable_tigs,
                    tiglap.variable_roots,
                    tiglap.shared_tig,
                    aln_len1,
                    within,
                    read.id,
                    read2.id
                )
                .unwrap();
                write!(
                    bufwriter,
                    "{}\t{}\t{}\n",
                    frac_shared, frac_shared_variable, within
                )
                .unwrap();
            }

            let span_read1 = tiglap.tig1_end - tiglap.tig1_start + 1;
            if frac_shared_variable > 0.8 && span_read1 as f64 > (read.tig_seq.len() as f64) * 0.9 {
                contained = true;
                log::trace!("{} {} CONTAINED", read.id, read2.id);
                contained_reads.insert(i);
            }
        }
        if !contained {
            outer_reads.push(i);
        }
    }

    log::info!("OUTER READS: {}", outer_reads.len());
    for i in outer_reads.iter() {
        log::info!("{}", all_reads_cat[*i].id);
    }



    return outer_reads

}

pub fn get_overlaps_outer_reads(
    all_reads_cat: &Vec<TigRead>,
    disjoint_set: &DisjointSet,
    used_ids: &FxHashSet<u64>,
    outer_reads: &Vec<usize>,
) -> Vec<TigdexOverlap> {
    let mut bufwriter = BufWriter::new(File::create("overlaps.txt").unwrap());
    let inverted_index_hashmap_outer = outer_reads
        .iter()
        .map(|&i| {
            all_reads_cat[i]
                .tig_seq
                .iter()
                .map(|&x| (x, i))
                .collect::<Vec<(u32, usize)>>()
        })
        .flatten()
        .fold(FxHashMap::default(), |mut acc, x| {
            acc.entry(x.0).or_insert(vec![]).push(x.1);
            acc
        });
    let mut overlaps = vec![];

    let mut compared_reads = FxHashSet::default();
    for (i, read) in all_reads_cat.iter().enumerate() {
        let mut index_count_map = FxHashMap::default();
        if !outer_reads.contains(&i) {
            continue;
        }
        for &y in read.tig_seq.iter() {
            if let Some(indices) = inverted_index_hashmap_outer.get(&y) {
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
            let read2 = &all_reads_cat[*index];
            let tiglap = disjoint_distance(&read.tig_seq, &read2.tig_seq, &disjoint_set, &used_ids, i, *index);
            let frac_shared_variable = tiglap.variable_tigs as f64 / tiglap.variable_roots as f64;
            writeln!(bufwriter,
                "intersect tigs: {} i: {} j: {} fsv: {} leni: {} {}-{} lenj: {} {}-{} shared_tigs: {} REVERSE:{}",
                count,
                i,
                *index,
                frac_shared_variable,
                read.tig_seq.len(),
                tiglap.tig1_start,
                tiglap.tig1_end,
                read2.tig_seq.len(),
                tiglap.tig2_start,
                tiglap.tig2_end,
                tiglap.shared_tig,
                tiglap.chain_reverse
            ).unwrap();
            if frac_shared_variable > 0.8 {
                overlaps.push(tiglap);
            }
            compared_reads.insert(sorted_readpair);
        }
    }

    return overlaps;
}

pub fn read_graph_from_overlaps(all_reads_cat: Vec<TigRead>, overlaps: &Vec<TigdexOverlap>) -> OverlapGraph
{
    let mut nodes = FxHashMap::default();
    let mut edges = vec![];

    for tiglap in overlaps.iter() {
        let i_tig = tiglap.tig1;
        let j_tig = tiglap.tig2;
        let read = &all_reads_cat[tiglap.tig1].tig_seq;
        let read2 = &all_reads_cat[tiglap.tig2].tig_seq;
        let frac_shared_variable = tiglap.variable_tigs as f64 / tiglap.variable_roots as f64;

        //check if end-to-end overlap
        let mut forward1 = false;
        let mut forward2 = false;
        let mut r1_r2 = true;
        let mut ol = false;

        if tiglap.chain_reverse {
            if tiglap.tig1_start < read.len() / 10 && tiglap.tig2_start < read2.len() / 10 {
                forward1 = false;
                forward2 = true;
                r1_r2 = true;
                ol = true;
            } else if tiglap.tig1_end > read.len() * 9 / 10
                && tiglap.tig2_end > read2.len() * 9 / 10
            {
                forward1 = true;
                forward2 = false;
                r1_r2 = true;
                ol = true;
            }
        } else {
            if tiglap.tig1_start < read.len() / 10 && tiglap.tig2_end > read2.len() * 9 / 10 {
                forward1 = true;
                forward2 = true;
                r1_r2 = false;
                ol = true;
            } else if tiglap.tig2_start < read2.len() / 10
                && tiglap.tig1_end > read.len() * 9 / 10
            {
                forward1 = true;
                forward2 = true;
                r1_r2 = true;
                ol = true;
            }
        }
        let aln_len1 = tiglap.tig1_end - tiglap.tig1_start + 1;
        let aln_len2 = tiglap.tig2_end - tiglap.tig2_start + 1;
        if ol && frac_shared_variable > 0.8 && aln_len1.max(aln_len2) > 20 {
            log::info!(
                "OVERLAP {} {} {} {} {}-{} {} {}-{}, REVERSE: {}",
                i_tig,
                j_tig,
                frac_shared_variable,
                read.len(),
                tiglap.tig1_start,
                tiglap.tig1_end,
                read2.len(),
                tiglap.tig2_start,
                tiglap.tig2_end,
                tiglap.chain_reverse
            );

            let start = if r1_r2 { i_tig } else { j_tig };
            let end = if r1_r2 { j_tig } else { i_tig };
            let new_read_overlap = ReadOverlapEdge {
                node1: start,
                node2: end,
                forward1,
                forward2,
                overlap_len_bases: 0,
                overlap_len_tigs: aln_len1.max(aln_len2),
                shared_tigs: tiglap.shared_tig,
                variable_tigs: tiglap.variable_tigs,
                variable_roots: tiglap.variable_roots,
            };

            edges.push(Some(new_read_overlap));
            let ind = edges.len() - 1;
            {
                let rd1 = nodes.entry(i_tig).or_insert(ReadData::default());
                rd1.index = i_tig;
                if r1_r2 && forward1 {
                    rd1.out_edges.push(ind)
                } else if !r1_r2 && !forward2 {
                    rd1.out_edges.push(ind)
                } else {
                    rd1.in_edges.push(ind)
                }
            }
            {
                let rd2 = nodes.entry(j_tig).or_insert(ReadData::default());
                rd2.index = j_tig;
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

    let graph = OverlapGraph { reads: all_reads_cat, nodes, edges };

    return graph;
}


pub fn print_graph_stdout(graph: &OverlapGraph, file: &str) {
    let mut bufwriter = BufWriter::new(File::create(file).unwrap());
    let all_reads_cat = &graph.reads;
    for edge in graph.edges.iter() {
        if let Some(edge) = edge {
            let i = edge.node1;
            let j = edge.node2;
            let read = &all_reads_cat[i];
            let read2 = &all_reads_cat[j];
            let frac_shared_variable = edge.variable_tigs as f64 / edge.variable_roots as f64;
            let forward1 = edge.forward1;
            let forward2 = edge.forward2;
            let aln_len = edge.overlap_len_tigs;

            if read.id.contains("junk") || read2.id.contains("junk") || read.id.contains("chimera") || read2.id.contains("chimera"){
                continue;
            }
            let name1 = read.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[0];
            let name2 = read2.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[0];
            let range1 = read.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[2];
            let range2 = read2.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[2];

            writeln!(bufwriter, 
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i, j, name1, name2, range1, range2, forward1, forward2
            );
        }
    }
}
