use std::path::Path;
use fxhash::hash64;
use crate::cli::Cli;
use fxhash::FxHashMap;
use crate::graph::*;
use statrs::distribution::{Binomial, DiscreteCDF};
use rayon::prelude::*;
use std::sync::Mutex;
use crate::types::*;
use serde::{Serialize, Deserialize};
use fxhash::FxHashSet;
use std::fs::File;
use std::io::{BufWriter, Write};
use crate::mapping::*;

pub struct OverlapConfig{
    pub hang1 : usize,
    pub hang2 : usize,
    pub forward1 : bool,
    pub forward2 : bool,
    pub overlap1_len : usize,
    pub overlap2_len : usize,
}

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

#[derive(Debug, Clone,PartialEq, Default, Hash, Eq, Serialize, Deserialize)]
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

impl OverlapTwinGraph{
    pub fn other_node(&self, node: NodeIndex, edge: &ReadOverlapEdgeTwin) -> NodeIndex {
        if edge.node1 == node {
            edge.node2
        } else {
            edge.node1
        }
    }

    //Initial construction allowed lax overlaps, now prune if 
    //it doesn't cause tipping. Can look into more sophisticated heuristics later.
    pub fn prune_lax_overlaps(&mut self){
        let mut edges_to_remove = FxHashSet::default();
        for i in 0..self.edges.len(){
            let edge_opt = &self.edges[i];
            if edge_opt.is_none(){
                continue;
            }

            let edge = edge_opt.as_ref().unwrap();
            // Good overlap
            if edge.diff_snpmers == 0{
                continue;
            }

            // Check that node1 and node2 have non-lax overlaps
            let node1 = self.nodes.get(&edge.node1).unwrap();
            let node2 = self.nodes.get(&edge.node2).unwrap();
            let dir1 = edge.node_edge_direction(&edge.node1);
            let dir2 = edge.node_edge_direction(&edge.node2);

            let edge_nd1 = node1.edges_direction(&dir1);
            let edge_nd2 = node2.edges_direction(&dir2);

            let mut node1_good_found = false;
            for edge_id in edge_nd1{
                if let Some(edge) = &self.edges[*edge_id]{
                    if edge.diff_snpmers == 0{
                        node1_good_found = true;
                        break;
                    }
                }
            }

            let mut node2_good_found = false;
            for edge_id in edge_nd2{
                if let Some(edge) = &self.edges[*edge_id]{
                    if edge.diff_snpmers == 0{
                        node2_good_found = true;
                        break;
                    }
                }
            }

            if node1_good_found || node2_good_found{
                edges_to_remove.insert(i);
                self.edges[i] = None;
            }
        }

        log::debug!("Pruning {} lax overlaps", edges_to_remove.len());
        self.remove_edges(edges_to_remove);

        let mut singleton_nodes_to_remove = vec![];
        for (i, node) in self.nodes.iter(){
            if node.both_edges().collect::<Vec<_>>().is_empty(){
                log::trace!("Node {} has no edges after lax cuts", node.index);
                singleton_nodes_to_remove.push(*i);
            }
        }

        self.remove_nodes(&singleton_nodes_to_remove, false);
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

            //Need this for singleton nodes; happens usually due to lax cuts
            if sorted_all_edges.is_empty() {
                continue;
            }

            for &(_, edge_id, _) in sorted_all_edges.iter() {
                let edge = self.edges[edge_id].as_ref().unwrap();
                let other_node = self.other_node(*node_id, edge);
                mark.insert(other_node, Mark::InPlay);
            }

            let mut mark_direction_map = FxHashMap::default();
            let longest = sorted_all_edges.last().unwrap().0 + FUZZ;

            // Iterating over the present node's edges
            //
            //            edge_id
            // node_id: o>----->o : other_node
            for &(length, edge_id, first_direction) in sorted_all_edges.iter() {
                let edge = self.edges[edge_id].as_ref().unwrap();
                let other_node = self.other_node(*node_id, edge);
                let present_edge_direction_other = edge.node_edge_direction(&other_node);
                
                //Allow self loops. This doesn't happen in practice for the string graph (yet),
                //but it's a valid operation that should not be reduced. 
                if other_node == *node_id {
                    continue;
                }

                // Iterating over the other node's edges
                //   edge_id
                //o>------>o : other_node
                //            edge_id2
                //         o >------>o : third_node
                for &(length2, edge_id2, new_edge_direction_other) in &self.get_edges_sorted_by_length(other_node) {
                    //outgoing edge for other node, requires incoming
                    if new_edge_direction_other == present_edge_direction_other {
                        continue;
                    }
                    
                    let next_edge = self.edges[edge_id2].as_ref().unwrap();
                    let third_node = self.other_node(other_node, next_edge);
                    let next_third_direction = next_edge.node_edge_direction(&third_node);

                    if mark.get(&third_node).unwrap() == &Mark::InPlay {
                        let lensum = if length + length2 < self.nodes[&other_node].base_length {0} else {length + length2 - self.nodes[&other_node].base_length};
                        if lensum <= longest + FUZZ || true
                        {
                            if let Some(true) =
                                mark.get(&third_node).map(|m| *m == Mark::InPlay)
                            {
                                mark.insert(third_node, Mark::Eliminated);
                                let valid_directions = mark_direction_map.entry(third_node).or_insert(FxHashSet::default());
                                valid_directions.insert((first_direction, next_third_direction));
                                log::trace!("Potential reduction from {} to {}, length1 {}, length2 {}, longets {}, lensum {}, edge_info1 {:?}, edge_info2 {:?}", node_id, third_node, length, length2, longest, lensum, &edge, &next_edge);
                            }
                        }
                    }
                }
            }
            // Step 3: Final pass to mark reduced edges

            // Possible reduction:
            //o>------>o 
            //         o >------>o 

            // Require direction concordance
            //o<---------------->o : NOT OK; 
            //o<----------------<o : NOT OK; 
            //o>---------------->o : OK
            for &edge_id in node_data.out_edges.iter().chain(node_data.in_edges.iter()){
                let edge = self.edges[edge_id].as_ref().unwrap();
                let direction_out_of_initial = self.edges[edge_id].as_ref().unwrap().node_edge_direction(node_id);

                let other_node = self.other_node(*node_id, edge);
                let direction_into_other = edge.node_edge_direction(&other_node);
                if let Some(Mark::Eliminated) = mark.get(&other_node) {
                    if let Some(valid_directions) = mark_direction_map.get(&other_node){
                        if valid_directions.contains(&(direction_out_of_initial, direction_into_other)){
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

    fn get_edges_sorted_by_length(&self, node_id: NodeIndex) -> Vec<(usize, EdgeIndex, Direction)> {
        let mut sorted_all_edges = vec![];
        let node_data = self.nodes.get(&node_id).unwrap();
        for (l, edge_ind) in node_data
            .out_edges
            .iter()
            .chain(node_data.in_edges.iter())
            .enumerate()
        {
            let direction;
            if l < node_data.out_edges.len() {
                direction = Direction::Outgoing; 
            } else {
                direction = Direction::Incoming;
            }
            let edge = self.edges[*edge_ind].as_ref().unwrap();
            let n1 = self.nodes.get(&edge.node1).unwrap();
            let n2 = self.nodes.get(&edge.node2).unwrap();
            let string_length = self.nodes[&n1.index].base_length
                + self.nodes[&n2.index].base_length
                - edge.overlap_len_bases;
            sorted_all_edges.push((string_length, *edge_ind, direction));
        }
        sorted_all_edges.sort_by(|a,b| a.0.cmp(&b.0));
        sorted_all_edges
    }
}

// Enum for marking the state of a node during processing
#[derive(PartialEq, Eq, Clone)]
enum Mark {
    Vacant,
    InPlay,
    Eliminated,
}

pub fn read_graph_from_overlaps_twin(all_reads_cat: &[TwinRead], overlaps: &[TwinOverlap], args: &Cli) -> OverlapTwinGraph
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
        let identity = id_est(twlap.shared_minimizers, twlap.diff_snpmers, args.c as u64);

        let aln_len1 = twlap.end1 - twlap.start1;
        let aln_len2 = twlap.end2 - twlap.start2;
        let overlap_hang_length = 750;
        if aln_len1.max(aln_len2) < 1000 {
            continue;
        }

        let mut overlap_possibilites_out = vec![];
        let mut overlap_possibilites_in = vec![];

        if twlap.chain_reverse {
            if twlap.start1 < overlap_hang_length && twlap.start2 < overlap_hang_length {
                let ol_config = OverlapConfig {
                    hang1: twlap.start1,
                    hang2: twlap.start2,
                    forward1: false,
                    forward2: true,
                    overlap1_len: aln_len1,
                    overlap2_len: aln_len2,
                };
                overlap_possibilites_in.push(ol_config);
            } 
            else if twlap.end1 + overlap_hang_length > read1.base_length && twlap.end2 + overlap_hang_length > read2.base_length 
            {
                let ol_config = OverlapConfig {
                    hang1: read1.base_length - twlap.end1 - 1,
                    hang2: read2.base_length - twlap.end2 - 1,
                    forward1: true,
                    forward2: false,
                    overlap1_len: aln_len1,
                    overlap2_len: aln_len2,
                };
                overlap_possibilites_out.push(ol_config);
            }
        } else {
            if twlap.start1 < overlap_hang_length && twlap.end2 + overlap_hang_length > read2.base_length {
                let ol_config = OverlapConfig {
                    hang1: twlap.start1,
                    hang2: read2.base_length - twlap.end2 - 1,
                    forward1: false,
                    forward2: false,
                    overlap1_len: aln_len1,
                    overlap2_len: aln_len2,
                };
                overlap_possibilites_in.push(ol_config);
            } 
            else if twlap.end1 + overlap_hang_length > read1.base_length  && twlap.start2 < overlap_hang_length {
                let ol_config = OverlapConfig {
                    hang1: read1.base_length - twlap.end1 - 1,
                    hang2: twlap.start2,
                    forward1: true,
                    forward2: true,
                    overlap1_len: aln_len1,
                    overlap2_len: aln_len2,
                };
                overlap_possibilites_out.push(ol_config);
            }
        }

        // Only allow one overlap per pair of reads for now. 
        //    o
        //   | |
        //    x
        // Can have two edges, x y + + and y x + +. Add at most 1 in edge and one out edge. 

        let best_overlap_forward = overlap_possibilites_out.iter().max_by_key(|x| (x.overlap1_len + x.overlap2_len) - (x.hang1 + x.hang2));
        let best_overlap_backward = overlap_possibilites_in.iter().max_by_key(|x| (x.overlap1_len + x.overlap2_len) - (x.hang1 + x.hang2));

        for best_overlap in [best_overlap_forward, best_overlap_backward]{
            if let Some(best_overlap) = best_overlap {
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

                let new_read_overlap = ReadOverlapEdgeTwin {
                    node1: i,
                    node2: j,
                    hang1: best_overlap.hang1,
                    hang2: best_overlap.hang2,
                    overlap1_len: best_overlap.overlap1_len,
                    overlap2_len: best_overlap.overlap2_len,
                    forward1: best_overlap.forward1,
                    forward2: best_overlap.forward2,
                    overlap_len_bases: best_overlap.overlap1_len.max(best_overlap.overlap2_len),
                    shared_minimizers: twlap.shared_minimizers,
                    diff_snpmers: twlap.diff_snpmers,
                    shared_snpmers: twlap.shared_snpmers
                };

                edges.push(Some(new_read_overlap));
                let ind = edges.len() - 1;
                {
                    let rd1 = nodes.entry(i).or_insert(ReadData::default());
                    rd1.index = i;
                    if best_overlap.forward1 {
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
                    if best_overlap.forward2 {
                        rd2.in_edges.push(ind)
                    } else {
                        rd2.out_edges.push(ind)
                    }
                    rd2.base_length = read2.base_length;
                    rd2.read_id = read2.id.clone();
                }
            }
        }
    }

    let mut graph = OverlapTwinGraph { nodes, edges };

    graph.prune_lax_overlaps();
    graph.transitive_reduction();

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
    let _file = Path::new(args.output_dir.as_str()).join("comparisons.txt");
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
            let _sorted_readpair = if i < index { (i, index) } else { (index, i) };

            //Only compare once. I think we get slightly different results if we
            //compare in both directinos, but this forces consistency. 
            if i < index{
                continue;
            }

            let read2 = &twin_reads[index];
            let twlaps = compare_twin_reads(&read, &read2, None, None, i, index, true);
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

                let same_strain_lax = same_strain(twlap.shared_minimizers, twlap.diff_snpmers, twlap.shared_snpmers, args.c as u64, args.snpmer_threshold, args.snpmer_error_rate);
                if same_strain_lax{
                    overlaps.lock().unwrap().push(twlap);
                }
            }
        }
    });

    return overlaps.into_inner().unwrap();

}

pub fn remove_contained_reads_twin<'a>(indices: Option<Vec<usize>>, twin_reads: &'a [TwinRead],  args: &Cli) -> Vec<usize>{
    let start = std::time::Instant::now();
    let downsample_factor = (args.contain_subsample_rate / args.c).max(1) as u64;
    log::info!("Building inverted index hashmap for all reads...");
    let inverted_index_hashmap =
        twin_reads
            .iter()
            .enumerate()
            .filter(|x| x.1.est_id.is_none() 
            || x.1.est_id.unwrap() > args.quality_value_cutoff)
            .fold(FxHashMap::default(), |mut acc, (i, x)| {
                for &y in x.minimizers.iter() {
                    if hash64(&y.1) > u64::MAX / downsample_factor {
                        continue;
                    }
                    acc.entry(y.1).or_insert(FxHashSet::default()).insert(i);
                }
                acc
            });
    //println!("Time to build inverted index hashmap: {:?}", start.elapsed());

    //open file for writing
    let name = if indices.is_none() { "all-cont.txt" } else { "subset-cont.txt" };
    let output_path = Path::new(args.output_dir.as_str()).join(name);
    let bufwriter_dbg = Mutex::new(BufWriter::new(File::create(output_path).unwrap()));
    let contained_reads = Mutex::new(FxHashSet::default());
    let outer_reads = Mutex::new(vec![]);

    let range = if let Some(indices) = indices {
        indices
    } else {
        (0..twin_reads.len()).into_iter().collect::<Vec<_>>()
    };

    log::info!("Removing contained reads...");
    range.into_par_iter().for_each(|i| {
        let mut contained = false;
        let read1 = &twin_reads[i];
        if read1.est_id.is_some() && read1.est_id.unwrap() < args.quality_value_cutoff {
            return;
        }

        let start = std::time::Instant::now();
        let mut index_count_map = FxHashMap::default();
        let mut index_range_map = FxHashMap::default();
        for &y in read1.minimizers.iter() {
            //Downsample to 100 compression factor
            if hash64(&y.1) > u64::MAX / downsample_factor{
                continue;
            }
            if let Some(indices) = inverted_index_hashmap.get(&y.1) {
                for &index in indices {
                    if index == i {
                        continue;
                    }
                    *index_count_map.entry(index).or_insert(0) += 1;
                    let range = index_range_map.entry(index).or_insert([u32::MAX,0]);
                    if (y.0 as u32) < range[0]{
                        range[0] = y.0 as u32;
                    }
                    if y.0 as u32 > range[1]{
                        range[1] = y.0 as u32;
                    }
                }
            }
        }

        let inverted_indexing_time = start.elapsed();
        let start = std::time::Instant::now();
        let mut top_indices = index_count_map.iter().collect::<Vec<_>>();

        top_indices.retain(|(index,_)| twin_reads[**index].base_length > read1.base_length);
        top_indices.retain(|(_,count)| **count > read1.minimizers.len() as u32 / 10 && **count > 5);
        top_indices.retain(|(index,_)| {
            let range = index_range_map.get(&index).unwrap();
            range[1] > range[0] && (range[1] - range[0]) as f64 + (60. * args.c as f64) > read1.base_length as f64 * 0.80
        });

        top_indices.sort_by(|a, b| {
            let lr_a = index_range_map.get(&a.0).unwrap();
            let length_a;
            if lr_a[0] > lr_a[1]{
                length_a = 0;
            }
            else{
                length_a = lr_a[1] - lr_a[0];
            }
            let lr_b = index_range_map.get(&b.0).unwrap();
            let length_b;
            if lr_b[0] > lr_b[1]{
                length_b = 0;
            }
            else{
                length_b = lr_b[1] - lr_b[0];
            }

            (length_b * b.1).cmp(&(length_a * a.1))
        });

        let sorting_filtering_time = start.elapsed();
        writeln!(bufwriter_dbg.lock().unwrap(), "CONTAIN: {} NUMBER OF HITS {}", &read1.id, top_indices.len()).unwrap();

        let start = std::time::Instant::now();
        let num_tries = 50;
        let mut num_fails = 0;
        let mut max_ol = 0;
        for (index, order_count) in top_indices.into_iter() {
            if contained{
                break;
            }
            let read2 = &twin_reads[*index];
            let twin_overlaps = compare_twin_reads(&read1, &read2, None, None, i, *index, true);
            if twin_overlaps.is_empty(){
                continue;
            }
            let twin_overlap = twin_overlaps.first().unwrap();
            let shared_minimizers = twin_overlap.shared_minimizers;
            let diff_snpmers = twin_overlap.diff_snpmers;

            let identity = id_est(shared_minimizers, diff_snpmers, args.c as u64);
            let same_strain = same_strain(shared_minimizers, diff_snpmers, twin_overlap.shared_snpmers, args.c as u64, args.snpmer_threshold_contain, 0.);

            let len1 = twin_overlap.end1 - twin_overlap.start1;
            let len2 = twin_overlap.end2 - twin_overlap.start2;
            let ol_len = len1;
            
            if ol_len > max_ol{
                max_ol = ol_len;
                num_fails = 0;
            }
            else{
                num_fails += 1;
                if num_fails > num_tries{
                    contained = false;
                    break;
                }
            }

            // only do this when log is config to trace
            if log::log_enabled!(log::Level::Debug) {
                writeln!(
                    bufwriter_dbg.lock().unwrap(),
                    "{} {}:{}-{} ----- {} {}:{}-{}   minis: {} shared_snps: {}, diff_snps: {}, identity {}, ol_len {}, read1len: {}, read2len: {}, order_count: {}",
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
                    len1.max(len2),
                    read1.base_length,
                    read2.base_length,
                    order_count

                )
                .unwrap();
            }


            //Can't be reflexive.
            if r1_contained_r2(&twin_overlap, read1, read2, same_strain, args.c){
                let comparing_time = start.elapsed();
                writeln!(bufwriter_dbg.lock().unwrap(), "{} {} CONTAINED. Times {} {} {}", read1.id, read2.id, inverted_indexing_time.as_micros(), sorting_filtering_time.as_micros(), comparing_time.as_micros()).unwrap();
                contained = true;
                contained_reads.lock().unwrap().insert(i);
                break;
            }
            
        }
        //println!("Comparing took: {:?}", start.elapsed());

        if !contained{
            outer_reads.lock().unwrap().push(i);
        }
    });
    log::info!("{} reads are contained", contained_reads.lock().unwrap().len());
    log::info!("Time elapsed for removing contained reads is: {:?}", start.elapsed());
    return outer_reads.into_inner().unwrap();
}

#[inline]
pub fn r1_contained_r2(twin_overlap: &TwinOverlap, read1: &TwinRead, read2: &TwinRead, same_strain: bool, c: usize) -> bool {
    let ol_len = twin_overlap.end1 - twin_overlap.start1;
    if ol_len as f64 + (30. * c as f64) > 0.95 * (read1.base_length as f64) && same_strain && read1.base_length < read2.base_length {
        return true;
    }
    return false;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Direction::{Incoming, Outgoing};

    struct MockEdge {
        pub i: usize,
        pub j: usize,
        pub d1: Direction,
        pub d2: Direction,
        pub diff_snpmers: usize,
    }
    
    impl MockEdge{
        pub fn new(i: usize, j: usize, d1: Direction, d2: Direction) -> Self{
            Self{
                i,
                j,
                d1,
                d2,
                diff_snpmers: 0
            }
        }
    }

    fn mock_graph_from_edges(edge_list: Vec<MockEdge>) -> OverlapTwinGraph {
        let mut nodes = FxHashMap::default();
        let mut edges = vec![];

        for edge in edge_list {
            let new_edge = ReadOverlapEdgeTwin {
                node1: edge.i,
                node2: edge.j,
                hang1: 0,
                hang2: 0,
                overlap1_len: 2000,
                overlap2_len: 2000,
                forward1: edge.d1 == Direction::Outgoing,
                forward2: edge.d2 == Direction::Incoming,
                overlap_len_bases: 2000,
                shared_minimizers: 100,
                diff_snpmers: edge.diff_snpmers,
                shared_snpmers: 10
            };

            edges.push(Some(new_edge));
            let ind = edges.len() - 1;
            {
                let nodes_len = nodes.len();
                let rd1 = nodes.entry(edge.i).or_insert(ReadData::default());
                rd1.base_length = 3000;
                rd1.read_id = nodes_len.to_string();
                rd1.index = edge.i;
                if edge.d1 == Direction::Outgoing {
                    rd1.out_edges.push(ind)
                } else {
                    rd1.in_edges.push(ind)
                }
            }
            {
                let nodes_len = nodes.len();
                let rd2 = nodes.entry(edge.j).or_insert(ReadData::default());
                rd2.index = edge.j;
                rd2.read_id = nodes_len.to_string();
                rd2.base_length = 3000;
                if edge.d2 == Direction::Incoming{
                    rd2.in_edges.push(ind)
                } else {
                    rd2.out_edges.push(ind)
                }
            }
        }

        OverlapTwinGraph { nodes, edges }
    }

    #[test]
    fn test_same_strain() {
        let minimizers = 100;
        let snp_diff = 10;
        let snp_shared = 10;
        let c = 10;
        let snpmer_threshold = 99.9;
        let snpmer_error_rate = 0.025;
        assert_eq!(same_strain(minimizers, snp_diff, snp_shared, c, snpmer_threshold, snpmer_error_rate), false);

        let minimizers = 100;
        let snp_diff = 1;
        let snp_shared = 0;
        assert_eq!(same_strain(minimizers, snp_diff, snp_shared, c, snpmer_threshold, snpmer_error_rate), true);

        let minimizers = 10000;
        let snp_diff = 5;
        let snp_shared = 0;
        assert_eq!(same_strain(minimizers, snp_diff, snp_shared, c, snpmer_threshold, snpmer_error_rate), true);

        let minimizers = 10;
        let snp_diff = 1;
        let snp_shared = 0;
        assert_eq!(same_strain(minimizers, snp_diff, snp_shared, c, snpmer_threshold, snpmer_error_rate), false);

        let minimizers = 100;
        let snp_diff = 5;
        let snp_shared = 10;
        let snpmer_error_rate = 0.50;
        assert_eq!(same_strain(minimizers, snp_diff, snp_shared, c, snpmer_threshold, snpmer_error_rate), true);
    }

    #[test]
    fn test_transitive_reduction_1() {
        //  >-------------> GOOD
        // o>----->o>----->o
        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Incoming),
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(0, 2, Outgoing, Incoming),
        ];
        let mut graph = mock_graph_from_edges(mock_edges1);
        assert_eq!(graph.edges.len(), 3);

        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(some_edges, 2);
    }

    #[test]
    fn test_transitive_reduction_2() {
        //  >-------------< BAD
        // o>----->o>----->o
        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Incoming),
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(0, 2, Outgoing, Outgoing),
        ];

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(some_edges, 3);
    }

    #[test]
    fn test_transitive_reduction_3() {
        //  <-------< BAD
        // o>-->o>-->o
        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Incoming),
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(0, 2, Incoming, Outgoing),
        ];

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(some_edges, 3);
    }

    #[test]
    fn test_transitive_reduction_4() {
        //  >-------> BAD
        // o>--<o>-->o
        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Outgoing),
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(0, 2, Outgoing, Incoming),
        ];

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(some_edges, 3);
    }

    #[test]
    fn test_transitive_reduction_5() {
        //  GOOD
        // 0<--<1>-->2
        //  <-------<
        // 1 to 0 and 2 to 0
        let mock_edges1 = vec![
            MockEdge::new(0, 1, Incoming, Outgoing),
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(0, 2, Incoming, Outgoing),
        ];

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(some_edges, 2);
    }

    #[test]
    fn test_transitive_reduction_cycle() {
        // BAD CYCLE
        // 0>--->1>-->2
        //  ^ ------|

        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Incoming),
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(2, 0, Outgoing, Incoming),
        ];

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(some_edges, 3);
    }

    #[test]
    fn test_transitive_reduction_out_triangle(){
        //      o 
        //     v  v
        //    /   \ 
        //   ^     ^
        //  o > - < o
        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Outgoing),
            MockEdge::new(1, 2, Outgoing, Outgoing),
            MockEdge::new(2, 0, Outgoing, Outgoing),
        ];

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(some_edges, 3);
    }

    #[test]
    fn test_transitive_reduction_in_triangle(){
        //      o 
        //     ^  ^
        //    /   \ 
        //   v     v
        //  o < - > o
        let mock_edges1 = vec![
            MockEdge::new(0, 1, Incoming, Incoming),
            MockEdge::new(1, 2, Incoming, Incoming),
            MockEdge::new(2, 0, Incoming, Incoming), 
        ];

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(some_edges, 3);
    }


    #[test]
    fn test_transitive_reduction_loops(){
        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Incoming),
            MockEdge::new(1, 0, Outgoing, Incoming),
        ];

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(some_edges, 2);

        let mock_edges2 = vec![
            MockEdge::new(0, 0, Outgoing, Incoming)
        ];

        let mut graph = mock_graph_from_edges(mock_edges2);
        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(some_edges, 1);
    }

    #[test]
    fn test_transitive_reduction_many_nodes(){
        //
        //0    1    2    3    4    5    6    7    8 
        //o>-->o>-->o>-->o>-->o>-->o>-->o>-->o>-->o
        // >------->o    o<--<o    o>------->o  : 1,2
        //     o>------->o : 3
        // >------------>o : 4

        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Incoming),
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(2, 3, Outgoing, Incoming),
            MockEdge::new(3, 4, Outgoing, Incoming),
            MockEdge::new(4, 5, Outgoing, Incoming),
            MockEdge::new(5, 6, Outgoing, Incoming),
            MockEdge::new(6, 7, Outgoing, Incoming),
            MockEdge::new(7, 8, Outgoing, Incoming),
            MockEdge::new(0, 2, Outgoing, Incoming),
            MockEdge::new(1, 3, Outgoing, Incoming),
            MockEdge::new(0, 3, Outgoing, Incoming),
            MockEdge::new(4, 3, Outgoing, Incoming),
            MockEdge::new(5, 7, Outgoing, Incoming),
        ];

        let mut graph = mock_graph_from_edges(mock_edges1);
        assert_eq!(graph.edges.len(), 13);
        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(some_edges, 9);
    }

    //For unitigging
    #[test]
    fn test_get_nonbranching_paths_1(){
        //  >-------------> GOOD
        // o>----->o>----->o
        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Incoming),
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(0, 2, Outgoing, Incoming),
        ];

        let graph = mock_graph_from_edges(mock_edges1);
        let nbp = graph.find_non_branching_paths();
        assert_eq!(nbp.len(), 3);

        let nbp_set = nbp.into_iter().map(|(x,y)| (x, y)).collect::<FxHashSet<_>>();

        let true_vec = vec![
            (vec![0], vec![]),
            (vec![1], vec![]),
            (vec![2], vec![])
        ];
        let true_set = true_vec.into_iter().collect::<FxHashSet<_>>();
        assert_eq!(nbp_set, true_set);
    }

    #[test]
    fn test_get_nonbranching_paths_cycle(){
        // o>----->o>----->o
        // ^               |
        //  <-------------<
        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Incoming),
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(2, 0, Outgoing, Incoming),
        ];

        let graph = mock_graph_from_edges(mock_edges1);
        let nbp = graph.find_non_branching_paths();
        assert_eq!(nbp.len(), 1);

        //Order nondeterministic, just check sizes
        //let true_vec = vec![
        //    (vec![0,1,2], vec![0,1]),
        //];
        assert_eq!(nbp[0].0.len(), 3);
        assert_eq!(nbp[0].1.len(), 2);


        // o>------->o
        // ^       |
        //  <-----<
        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Incoming),
            MockEdge::new(1, 0, Outgoing, Incoming),
        ];

        let graph = mock_graph_from_edges(mock_edges1);
        let nbp = graph.find_non_branching_paths();
        assert_eq!(nbp.len(), 1);

        //Order nondeterministic, just check sizes
        //let true_vec = vec![
        //    (vec![0,1,2], vec![0,1]),
        //];
        assert_eq!(nbp[0].0.len(), 2);
        assert_eq!(nbp[0].1.len(), 1);
    }

    #[test]
    fn test_get_nonbranching_paths_self_cycle(){
        // o>--
        // ^  |
        //  
        let mock_edges1 = vec![
            MockEdge::new(0, 0, Outgoing, Incoming),
        ];

        let graph = mock_graph_from_edges(mock_edges1);
        let nbp = graph.find_non_branching_paths();
        assert_eq!(nbp.len(), 1);

        //Order nondeterministic, just check sizes
        //let true_vec = vec![
        //    (vec![0,1,2], vec![0,1]),
        //];
        assert_eq!(nbp[0].0.len(), 1);
        assert_eq!(nbp[0].1.len(), 0);
    }

    #[test]
    fn test_prune_lax_overlaps_nofp1(){
        // o>--x-->o>----->o
        let mut mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Incoming),
            MockEdge::new(1, 2, Outgoing, Incoming),
        ];

        mock_edges1[0].diff_snpmers = 2;

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.prune_lax_overlaps();
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 2);

    }

    #[test]
    fn test_prune_lax_overlaps_nofp2(){
        // o<--x-->o>---x-<o
        let mut mock_edges1 = vec![
            MockEdge::new(0, 1, Incoming, Incoming),
            MockEdge::new(1, 2, Incoming, Outgoing),
        ];

        mock_edges1[0].diff_snpmers = 2;
        mock_edges1[1].diff_snpmers = 2;

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.prune_lax_overlaps();
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 2);
    }

    #[test]
    fn test_prune_lax_overlaps_simple(){
        // 1 -- > 2
        //   ↘ (x)
        // 3 --> 4
        let mut mock_edges1 = vec![
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(1, 4, Outgoing, Incoming),
            MockEdge::new(3, 4, Outgoing, Incoming),
        ];

        mock_edges1[1].diff_snpmers = 5;

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.prune_lax_overlaps();
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 2);

        // Change orientations
        // 1 >--< 2
        //    ↘ (x)
        // 3 <--> 4
        let mut mock_edges1 = vec![
            MockEdge::new(1, 2, Outgoing, Outgoing),
            MockEdge::new(1, 4, Outgoing, Incoming),
            MockEdge::new(3, 4, Incoming, Incoming),
        ];

        mock_edges1[1].diff_snpmers = 5;

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.prune_lax_overlaps();
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 2);
    }

    #[test]
    fn test_prune_lax_overlaps_cross(){
        // 1 -x- > 2
        //   ↗↘ 
        // 3 -x-> 4
        let mut mock_edges1 = vec![
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(1, 4, Outgoing, Incoming),
            MockEdge::new(3, 2, Outgoing, Incoming),
            MockEdge::new(3, 4, Outgoing, Incoming),
        ];

        mock_edges1[0].diff_snpmers = 5;
        mock_edges1[3].diff_snpmers = 5;

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.prune_lax_overlaps();
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 2);

        // Change orientations
        // 1 >--< 2
        //    ↘ 
        // 3 <--> 4
        let mut mock_edges1 = vec![
            MockEdge::new(1, 2, Outgoing, Outgoing),
            MockEdge::new(1, 4, Outgoing, Incoming),
            MockEdge::new(3, 2, Incoming, Outgoing),
            MockEdge::new(3, 4, Incoming, Incoming),
        ];

        mock_edges1[0].diff_snpmers = 5;
        mock_edges1[3].diff_snpmers = 5;

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.prune_lax_overlaps();
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 2);
    }
}
