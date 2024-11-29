use std::path::Path;
use fxhash::hash64;
use crate::cli::Cli;
use fxhash::FxHashMap;
use crate::graph::*;
use statrs::distribution::{Binomial, DiscreteCDF};
use rayon::prelude::*;
use std::sync::Mutex;
use std::sync::RwLock;
use crate::types::*;
use fxhash::FxHashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use crate::mapping::*;

#[derive(Debug, Clone, PartialEq, Default, Hash, Eq)]
pub struct ReadData {
    pub index: NodeIndex,
    pub in_edges: Vec<EdgeIndex>,
    pub out_edges: Vec<EdgeIndex>,
    pub base_length: usize,
    pub read_id: String
}

impl GraphNode for ReadData {
    fn in_edges(&self) -> &[EdgeIndex] { &self.in_edges }
    fn out_edges(&self) -> &[EdgeIndex] { &self.out_edges }

    fn both_edges(&self) -> impl Iterator<Item = &EdgeIndex> {
        self.in_edges.iter().chain(self.out_edges.iter())
    }

    fn in_edges_mut(&mut self) -> &mut Vec<EdgeIndex> {
        &mut self.in_edges
    }

    fn out_edges_mut(&mut self) -> &mut Vec<EdgeIndex> {
        &mut self.out_edges
    }
}

#[derive(Debug, Clone,PartialEq, Default, Hash, Eq)]
pub struct ReadOverlapEdgeTwin {
    pub node1: NodeIndex,
    pub node2: NodeIndex,
    pub hang1: usize,
    pub hang2: usize,
    pub forward1: bool,
    pub forward2: bool,
    pub overlap1_len: usize,
    pub overlap2_len: usize,
    pub overlap_len_bases: usize,
    pub shared_minimizers: usize,
    pub diff_snpmers: usize,
    pub shared_snpmers: usize
}

impl GraphEdge for ReadOverlapEdgeTwin {
    fn node1(&self) -> NodeIndex { self.node1 }
    fn node2(&self) -> NodeIndex { self.node2 }
    fn orientation1(&self) -> bool {
        self.forward1
    }
    fn orientation2(&self) -> bool {
        self.forward2
    }
}

impl ReadOverlapEdgeTwin{
    pub fn get_orientation(&self, node1: NodeIndex, node2: NodeIndex) -> (bool, bool){
        if self.node1 == node1 && self.node2 == node2{
            (self.forward1, self.forward2)
        }
        else if self.node1 == node2 && self.node2 == node1{
            (!self.forward2, !self.forward1)
        }
        else{
            dbg!(self, node1, node2);
            panic!("Nodes do not match edge")
        }
    }
}

pub type OverlapTwinGraph = BidirectedGraph<ReadData, ReadOverlapEdgeTwin>;

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

            for &(_, edge_id, _) in sorted_all_edges.iter() {
                let edge = self.edges[edge_id].as_ref().unwrap();
                let other_node = self.other_node(*node_id, edge);
                mark.insert(other_node, Mark::InPlay);
            }

            let mut mark_direction_map = FxHashMap::default();
            let longest = sorted_all_edges.last().unwrap().0 + FUZZ;

            // Iterating over the present node's edges
            for &(length, edge_id, _) in sorted_all_edges.iter() {
                let edge = self.edges[edge_id].as_ref().unwrap();
                let other_node = self.other_node(*node_id, edge);
                let outgoing = self.outgoing_edge(other_node, edge_id);
                let direction_out_of_initial = self.edges[edge_id].as_ref().unwrap().node_edge_direction(node_id);

                // Iterating over the other node's edges
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
                        let lensum = if length + length2 < self.nodes[&other_node].base_length {0} else {length + length2 - self.nodes[&other_node].base_length};
                        if lensum <= longest + FUZZ || true
                        {
                            if let Some(true) =
                                mark.get(&third_node).map(|m| *m == Mark::InPlay)
                            {
                                mark.insert(third_node, Mark::Eliminated);
                                let valid_directions = mark_direction_map.entry(third_node).or_insert(FxHashSet::default());
                                valid_directions.insert(direction_out_of_initial);
                                log::trace!("Potential reduction from {} to {}, length1 {}, length2 {}, longets {}, lensum {}, edge_info1 {:?}, edge_info2 {:?}, outgoing_edge {}, direction_out {}", node_id, third_node, length, length2, longest, lensum, &edge, &next_edge, outer2, outgoing);
                            }
                        }
                    }
                }
            }
            // Step 3: Final pass to mark reduced edges
            for &edge_id in node_data.out_edges.iter().chain(node_data.in_edges.iter()){
                let edge = self.edges[edge_id].as_ref().unwrap();
                let other_node = self.other_node(*node_id, edge);
                let direction_out_of_initial = self.edges[edge_id].as_ref().unwrap().node_edge_direction(node_id);
                if let Some(Mark::Eliminated) = mark.get(&other_node) {
                    if let Some(valid_directions) = mark_direction_map.get(&other_node){
                        if valid_directions.contains(&direction_out_of_initial){
                            reduce[edge_id] = true;
                            log::trace!("Reduced from {} to {}. INFO:{:?}", node_id, other_node, &edge);
                        }
                    }
                }
                mark.insert(other_node, Mark::Vacant);
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

            if node_data.out_edges.is_empty() && node_data.in_edges.is_empty() {
                log::trace!("Node {} has no edges", node_data.index);
            }
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
            let string_length = self.nodes[&n1.index].base_length
                + self.nodes[&n2.index].base_length
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

pub fn read_graph_from_overlaps_twin(all_reads_cat: &Vec<TwinRead>, overlaps: &Vec<TwinOverlap>, args: &Cli) -> OverlapTwinGraph
{
    let raw_edges_file = Path::new(args.output_dir.as_str()).join("edges_raw.txt");
    let mut bufwriter = BufWriter::new(File::create(raw_edges_file).unwrap());
    let mut nodes = FxHashMap::default();
    let mut seen_overlaps = FxHashSet::default();
    let mut edges = vec![];

    for twlap in overlaps.iter() {
        let i = twlap.i1;
        let j = twlap.i2;
        let sorted_ij = if i < j { (i, j) } else { (j, i) };
        if seen_overlaps.contains(&sorted_ij) {
            continue;
        }
        else{
            seen_overlaps.insert(sorted_ij);
        }
        let read1 = &all_reads_cat[i];
        let read2 = &all_reads_cat[j];

        //check if end-to-end overlap
        let mut forward1 = false;
        let mut forward2 = false;
        let mut r1_r2 = true;
        let mut ol = false;

        let mut hang1 = 0;
        let mut hang2 = 0;

        let identity = id_est(twlap.shared_minimizers, twlap.diff_snpmers, args.c as u64);
        let same_strain = same_strain(twlap.shared_minimizers, twlap.diff_snpmers, twlap.shared_snpmers, args.c as u64, args.snpmer_threshold, args.snpmer_error_rate);

        if !same_strain {
            continue;
        }

        let overlap_hang_length = 750;

        if twlap.chain_reverse {
            if twlap.start1 < overlap_hang_length && twlap.start2 < overlap_hang_length {
                hang1 = twlap.start1;
                hang2 = twlap.start2;
                forward1 = false;
                forward2 = true;
                r1_r2 = true;
                ol = true;
            } 
            else if twlap.end1 + overlap_hang_length > read1.base_length && twlap.end2 + overlap_hang_length > read2.base_length 
            {
                hang1 = read1.base_length - twlap.end1;
                hang2 = read2.base_length - twlap.end2;
                forward1 = true;
                forward2 = false;
                r1_r2 = true;
                ol = true;
            }
        } else {
            if twlap.start1 < overlap_hang_length && twlap.end2 + overlap_hang_length > read2.base_length {
                hang1 = twlap.start1;
                hang2 = read2.base_length - twlap.end2;
                forward1 = true;
                forward2 = true;
                r1_r2 = false;
                ol = true;
            } 
            else if twlap.end1 + overlap_hang_length > read1.base_length  && twlap.start2 < overlap_hang_length {
                hang1 = read1.base_length - twlap.end1;
                hang2 = twlap.start2;
                forward1 = true;
                forward2 = true;
                r1_r2 = true;
                ol = true;
            }
        }

        let aln_len1 = twlap.end1 - twlap.start1 + 1;
        let aln_len2 = twlap.end2 - twlap.start2 + 1;
        if ol && aln_len1.max(aln_len2) > 1000 {
            writeln!(bufwriter,
                "READ EDGE: {}:{} {}:{} ID: {} SNP_SHARE:{}, SNP_DIFF:{} OL_INFO:{} {}-{} {} {}-{}, REVERSE: {}",
                &read1.id,
                i,
                &read2.id,
                j,
                identity * 100.,
                twlap.shared_snpmers,
                twlap.diff_snpmers,
                read1.base_length,
                twlap.start1,
                twlap.end1,
                read2.base_length,
                twlap.start2,
                twlap.end2,
                twlap.chain_reverse
            ).unwrap();

            let start = if r1_r2 { i } else { j };
            let end = if r1_r2 { j } else { i };
            let new_read_overlap = ReadOverlapEdgeTwin {
                node1: start,
                node2: end,
                hang1: hang1,
                hang2: hang2,
                overlap1_len: aln_len1,
                overlap2_len: aln_len2,
                forward1,
                forward2,
                overlap_len_bases: aln_len1.max(aln_len2).min(read1.base_length.min(read2.base_length)),
                shared_minimizers: twlap.shared_minimizers,
                diff_snpmers: twlap.diff_snpmers,
                shared_snpmers: twlap.shared_snpmers
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
                rd1.base_length = read1.base_length;
                rd1.read_id = read1.id.clone();
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
                rd2.base_length = read2.base_length;
                rd2.read_id = read2.id.clone();
            }
        }
    }

    let graph = OverlapTwinGraph { nodes, edges };

    return graph;
}


pub fn print_graph_stdout<T>(graph: &OverlapTwinGraph, file: T) where T: AsRef<std::path::Path> {
    let mut bufwriter = BufWriter::new(File::create(file).unwrap());
    for (l, edge) in graph.edges.iter().enumerate() {
        if let Some(edge) = edge {
            let i = edge.node1;
            let j = edge.node2;
            let read = &graph.nodes[&i];
            let read2 = &graph.nodes[&j];
            let forward1 = edge.forward1;
            let forward2 = edge.forward2;
            let _aln_len = edge.overlap_len_bases;

            if read.read_id.contains("junk") || read2.read_id.contains("junk"){
                log::trace!("{} {} edge {} CHIMERA OR JUNK", i,j,l);
                continue;
            }

            let res1 = parse_badread(&read.read_id);
            let res2 = parse_badread(&read2.read_id);
            let (name1, range1) = res1.unwrap();
            let (name2, range2) = res2.unwrap();
            
            writeln!(bufwriter, 
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i, j, name1, name2, range1, range2, forward1, forward2
            ).unwrap();
        }
    }
}


pub fn get_overlaps_outer_reads_twin(twin_reads: &[TwinRead], outer_read_indices: &[usize], args: &Cli) -> Vec<TwinOverlap>{

    let overlaps_file = Path::new(args.output_dir.as_str()).join("overlaps.txt");
    let bufwriter = Mutex::new(BufWriter::new(File::create(overlaps_file).unwrap()));
    let vec_format = twin_reads.iter().enumerate().filter(|(i, _)| outer_read_indices.contains(i)).map(|(i, x)| (i,&x.minimizers)).collect::<Vec<_>>();

    let mut inverted_index_hashmap_outer = FxHashMap::default();
    for (i, x) in vec_format.iter(){
        for &y in x.iter(){
            inverted_index_hashmap_outer.entry(y.1).or_insert(vec![]).push(i);
        }
    }
    
    let overlaps = Mutex::new(vec![]);

    outer_read_indices.into_par_iter().for_each(|&i| { 
        let read = &twin_reads[i];
        let mut index_count_map = FxHashMap::default();
        for &y in read.minimizers.iter() {
            if let Some(indices) = inverted_index_hashmap_outer.get(&y.1) {
                for &&index in indices {

                    //No self comparisons
                    if index == i {
                        continue;
                    }
                    *index_count_map.entry(index).or_insert(0) += 1;
                }
            }
        }
        //sort
        let mut top_indices = index_count_map.into_iter().filter(|(_, count)| *count > 500 / args.c).collect::<Vec<_>>();
        top_indices.sort_by(|a, b| b.1.cmp(&a.1));

        for &(index, _) in top_indices.iter(){
            let sorted_readpair = if i < index { (i, index) } else { (index, i) };

            //Only compare once. I think we get slightly different results if we
            //compare in both directinos, but this forces consistency. 
            if i < index{
                continue;
            }

            let read2 = &twin_reads[index];
            let twlaps = compare_twin_reads(&read, &read2, None, None, i, index);
            for twlap in twlaps{
                let identity = id_est(twlap.shared_minimizers, twlap.diff_snpmers, args.c as u64);
                writeln!(bufwriter.lock().unwrap(),
                    "{} {} fsv:{} leni:{} {}-{} lenj:{} {}-{} shared_mini:{} REVERSE:{}, snp_diff:{} snp_share:{}, intersect:{:?}, snp_both:{:?}, r1 {} r2 {}",
                    i,
                    index,
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
                    twlap.intersect,
                    twlap.snpmers_in_both,
                    &read.id,
                    &read2.id


                ).unwrap();
                if same_strain(twlap.shared_minimizers, twlap.diff_snpmers, twlap.shared_snpmers, args.c as u64, args.snpmer_threshold, args.snpmer_error_rate) {
                    overlaps.lock().unwrap().push(twlap);
                }
            }
        }
    });

    return overlaps.into_inner().unwrap();

}

pub fn remove_contained_reads_twin<'a>(indices: Option<Vec<usize>>, twin_reads: &'a [TwinRead],  args: &Cli) -> Vec<usize>{
    //let start = std::time::Instant::now();
    let inverted_index_hashmap =
        twin_reads
            .iter()
            .enumerate()
            .filter(|x| x.1.est_id.is_none() 
            || x.1.est_id.unwrap() > args.quality_value_cutoff)
            .fold(FxHashMap::default(), |mut acc, (i, x)| {
                for &y in x.minimizers.iter() {
                    if hash64(&y.1) > u64::MAX/3 {
                        continue;
                    }
                    acc.entry(y.1).or_insert(vec![]).push(i);
                }
                acc
            });
    //println!("Time to build inverted index hashmap: {:?}", start.elapsed());

    //open file for writing
    let name = if indices.is_none() { "all-cont.txt" } else { "subset-cont.txt" };
    let bufwriter_dbg = Mutex::new(BufWriter::new(File::create(name).unwrap()));
    let contained_reads = Mutex::new(FxHashSet::default());
    let outer_reads = Mutex::new(vec![]);

    let range = if let Some(indices) = indices {
        indices
    } else {
        (0..twin_reads.len()).into_iter().collect::<Vec<_>>()
    };
    range.into_par_iter().for_each(|i| {
        let mut contained = false;
        let read1 = &twin_reads[i];
        if read1.est_id.is_some() && read1.est_id.unwrap() < args.quality_value_cutoff {
            return;
        }
        //let start = std::time::Instant::now();
        let mut index_count_map = FxHashMap::default();
        for &y in read1.minimizers.iter() {
            if hash64(&y.1) > u64::MAX/3 {
                continue;
            }
            if let Some(indices) = inverted_index_hashmap.get(&y.1) {
                for &index in indices {
                    if index == i {
                        continue;
                    }
                    *index_count_map.entry(index).or_insert(0) += 1;
                }
            }
        }
        //println!("Querying took: {:?}", start.elapsed());

        //let start = std::time::Instant::now();
        let mut top_indices = index_count_map.iter().collect::<Vec<_>>();
        top_indices.sort_by(|a, b| b.1.cmp(a.1));
        let top_indices = top_indices.into_iter()
        .filter(|(_,count)| **count > read1.minimizers.len() as i32 / 20 && **count > 5)
        .filter(|(index,_)| twin_reads[**index].base_length > read1.base_length)
        .collect::<Vec<_>>();
        //println!("Sorting took: {:?}", start.elapsed());


        //let start = std::time::Instant::now();
        for (index, _) in top_indices.into_iter() {
            if contained{
                break;
            }
            let read2 = &twin_reads[*index];
            let twin_overlaps = compare_twin_reads(&read1, &read2, None, None, i, *index);
            for twin_overlap in twin_overlaps{
                let shared_minimizers = twin_overlap.shared_minimizers;
                let diff_snpmers = twin_overlap.diff_snpmers;

                let identity = id_est(shared_minimizers, diff_snpmers, args.c as u64);
                let same_strain = same_strain(shared_minimizers, diff_snpmers, twin_overlap.shared_snpmers, args.c as u64, args.snpmer_threshold, args.snpmer_error_rate);

                let len1 = twin_overlap.end1 - twin_overlap.start1;
                let len2 = twin_overlap.end2 - twin_overlap.start2;
                let ol_len = len1.max(len2);

                // only do this when log is config to trace
                if log::log_enabled!(log::Level::Trace) {
                    writeln!(
                        bufwriter_dbg.lock().unwrap(),
                        "{} {}:{}-{} ----- {} {}:{}-{}   minis: {} shared_snps: {}, diff_snps: {}, identity {}, ol_len {}, read1len: {}, read2len: {}",
                        read1.id.split_ascii_whitespace().next().unwrap(),
                        &i,
                        twin_overlap.start1,
                        twin_overlap.end1,
                        read2.id.split_ascii_whitespace().next().unwrap(),
                        &index,
                        twin_overlap.start2,
                        twin_overlap.end2,
                        twin_overlap.shared_minimizers,
                        twin_overlap.shared_snpmers,
                        twin_overlap.diff_snpmers,
                        identity,
                        ol_len,
                        read1.base_length,
                        read2.base_length

                    )
                    .unwrap();
                }


                //Can't be reflexive.
                if ol_len as f64 + (30. * args.c as f64) > 0.95 * (read1.base_length as f64) && same_strain && read1.base_length < read2.base_length {
                    writeln!(bufwriter_dbg.lock().unwrap(), "{} {} CONTAINED", read1.id, read2.id).unwrap();
                    contained = true;
                    contained_reads.lock().unwrap().insert(i);
                    break;
                }
            }
        }
        //println!("Comparing took: {:?}", start.elapsed());

        if !contained{
            outer_reads.lock().unwrap().push(i);
        }
    });
    return outer_reads.into_inner().unwrap();
}

pub fn binomial_test(n: u64, k: u64, p: f64) -> f64 {
    // n: number of trials
    // k: number of successes
    // p: probability of success

    // Create a binomial distribution
    let binomial = Binomial::new(p, n).unwrap();

    // Calculate the probability of observing k or more successes
    let p_value = 1.0 - binomial.cdf(k);

    p_value
}

pub fn same_strain(minimizers: usize, snp_diff: usize, snp_shared: usize,  c: u64, snpmer_threshold: f64, snpmer_error_rate: f64) -> bool {
    let identity = id_est(minimizers, snp_diff, c);
    let high_id;
    if identity >= snpmer_threshold / 100.{
        high_id = true;
    } else {
        high_id = false;
    }
    let p_val = binomial_test((snp_diff + snp_shared) as u64, snp_diff as u64, snpmer_error_rate);
    let miscalled_snpmers;
    if p_val > 0.05{
        miscalled_snpmers = true;
    } else {
        miscalled_snpmers = false; 
    }

    return high_id || miscalled_snpmers;
}
