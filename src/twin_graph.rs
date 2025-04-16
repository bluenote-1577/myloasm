use fxhash::hash64;
use crate::cli::Cli;
use crate::constants::{MINIMIZER_END_NTH_OVERLAP, OVERLAP_HANG_LENGTH};
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
use crate::map_processing;
use std::path::PathBuf;
use flate2::write::GzEncoder;
use flate2::Compression;


#[derive(Debug, Clone, PartialEq, Default, Hash, Eq, Serialize, Deserialize)]
pub struct OverlapConfig{
    pub hang1 : usize,
    pub hang2 : usize,
    pub forward1 : bool,
    pub forward2 : bool,
    pub overlap1_len : usize,
    pub overlap2_len : usize,
    pub read_i: usize,
    pub read_j: usize,
    pub shared_mini: usize,
    pub shared_snpmer: usize,
    pub diff_snpmer: usize,
    pub contained: bool,
    pub large_indel: bool
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

#[derive(Debug, Clone, PartialEq, Default, Eq)]
pub struct SnpmerOverlapConfig{
    pub snpmer_overlap_set: FxHashSet<(u32, bool)>,
    pub edge: EdgeIndex,
    pub length: usize
}

impl SnpmerOverlapConfig{
    pub fn greater_than(&self, other: &SnpmerOverlapConfig) -> bool {
        if self.length > other.length {
            for pos_match in other.snpmer_overlap_set.iter(){
                if !self.snpmer_overlap_set.contains(pos_match){
                    return false;
                }
            }
            return true;
        }
        return false;
    }

    pub fn soft_better_than(&self, other: &SnpmerOverlapConfig) -> bool {
        if self.length > other.length {
            for pos_match in self.snpmer_overlap_set.iter(){
                // Is not better than if a mismatch in self is a match in other
                if !pos_match.1{
                    let pos_match_ideal = (pos_match.0, true);
                    if other.snpmer_overlap_set.contains(&pos_match_ideal){
                        return false;
                    }
                }
            }
            return true;
        }
        return false;
    }

    pub fn sort(to_sort: &mut Vec<SnpmerOverlapConfig>){
        to_sort.sort_by(|a,b| (a.snpmer_overlap_set.len(), a.length).cmp(&(b.snpmer_overlap_set.len(), b.length)));
    }

    pub fn get_maximal(overlaps: &[SnpmerOverlapConfig]) -> Vec<usize> {
        let mut maximal = vec![];
        // look at largest first
        for i in (0..overlaps.len()).rev(){
            let mut is_maximal = true;
            for j in (i..overlaps.len()).rev(){
                if i == j{
                    continue;
                }
                if overlaps[j].greater_than(&overlaps[i]){
                    is_maximal = false;
                    break;
                }
            }
            if is_maximal{
                maximal.push(i);
            }
        }
        return maximal;
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
    pub shared_snpmers: usize,
    pub large_indel: bool
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
    fn edge_id_est(&self, c: usize) -> f64{
        id_est(self.shared_minimizers, self.diff_snpmers, c as u64, self.large_indel)
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

    // fn get_snpmer_overlap_set_vec(&self, node: &ReadData, direction: &Direction) -> Vec<SnpmerOverlapConfig>{
    //     let mut snpmer_overlap_set_vec = vec![];
    //     let edges = node.edges_direction(&direction);
    //     //Get snpmer_overlap_set_vec
    //     for edge_id in edges.iter(){
    //         let mut snpmer_overlap_set = FxHashSet::default();
    //         let edge = self.edges[*edge_id].as_ref().unwrap();
    //         let mut min_pos = u32::MAX;
    //         let mut max_pos = 0;
    //         for snpmer_hit in edge.snpmer_hits.iter(){
    //             let equal = snpmer_hit.bases.0 == snpmer_hit.bases.1;
    //             let pos;
    //             if node.index == edge.node1(){
    //                 pos = snpmer_hit.pos1;
    //             }
    //             else{
    //                 pos = snpmer_hit.pos2;
    //             }
    //             snpmer_overlap_set.insert((pos, equal));
    //             if pos < min_pos{
    //                 min_pos = pos;
    //             }
    //             if pos > max_pos{
    //                 max_pos = pos;
    //             }
    //         }
            
    //         let length = if min_pos == u32::MAX {0} else { max_pos - min_pos};
    //         let config = SnpmerOverlapConfig{
    //             snpmer_overlap_set,
    //             edge: *edge_id,
    //             length: length as usize
    //         };
    //         snpmer_overlap_set_vec.push(config);
    //     }
    //     SnpmerOverlapConfig::sort(&mut snpmer_overlap_set_vec);
    //     return snpmer_overlap_set_vec;
    // }

    // pub fn get_maximal_overlaps(&mut self, args: &Cli){
    //     log::info!("Getting maximal overlaps...");
    //     let edges_to_remove = Mutex::new(FxHashSet::default());
    //     self.nodes.iter().for_each(|(_, node)| {
    //         for direction in [Direction::Incoming, Direction::Outgoing]{
    //             let snpmer_overlap_set_vec = self.get_snpmer_overlap_set_vec(node, &direction);
    //             let maximal = SnpmerOverlapConfig::get_maximal(&snpmer_overlap_set_vec);
    //             let mut optimal_edge_ids = vec![];
    //             let edges = node.edges_direction(&direction);
    //             log::trace!("Read {} ({}), direction: {:?}", node.read_id, node.index, direction);
    //             for i in 0..edges.len(){
    //                 let edge_id = snpmer_overlap_set_vec[i].edge;
    //                 let edge = &self.edges[edge_id].as_ref().unwrap();
    //                 if maximal.contains(&i){
    //                     log::trace!("MAXIMAL ({}-{}) snp_diff:{} snp_share:{} ol:{}", 
    //                     &edge.node1(), edge.node2(), edge.diff_snpmers, edge.shared_snpmers, edge.overlap_len_bases);
    //                 }
    //                 else{
    //                     log::trace!("MINIMAL ({}-{}) snp_diff:{} snp_share:{} ol:{}", 
    //                     &edge.node1(), edge.node2(), edge.diff_snpmers, edge.shared_snpmers, edge.overlap_len_bases);
    //                     edges_to_remove.lock().unwrap().insert(i);
    //                 }
    //             }
    //             for index in maximal.iter(){
    //                 let edge_id = snpmer_overlap_set_vec[*index].edge;
    //                 let mut optimal = true;
    //                 for index2 in maximal.iter(){
    //                     if index == index2{
    //                         continue;
    //                     }
    //                     if snpmer_overlap_set_vec[*index2].soft_better_than(&snpmer_overlap_set_vec[*index]){
    //                         optimal = false;

    //                         edges_to_remove.lock().unwrap().insert(edge_id);
    //                         break;
    //                     }
    //                 }
    //                 if optimal{
    //                     let edge = &self.edges[edge_id].as_ref().unwrap();
    //                     optimal_edge_ids.push(edge_id);
    //                     let mut print_snpvec = snpmer_overlap_set_vec[*index].snpmer_overlap_set.iter().collect::<Vec<_>>();
    //                     print_snpvec.sort();
    //                     log::trace!("OPTIMAL ({}-{}) snp_diff:{} snp_share:{} ol:{} vec:{:?}", edge.node1(), edge.node2(), edge.diff_snpmers, edge.shared_snpmers, edge.overlap_len_bases, print_snpvec);
    //                 }
    //             }

    //             if optimal_edge_ids.len() == 0{
    //                 continue
    //             }

    //             let mut optimal_edges_with_fsv = optimal_edge_ids.iter().map(|x|{
    //                  let edge = &self.edges[*x].as_ref().unwrap();
    //                  let fsv = id_est(edge.shared_minimizers, edge.diff_snpmers, args.c as u64);
    //                  (*x, fsv)
    //             }).collect::<Vec<_>>();
                
    //             optimal_edges_with_fsv.sort_by(|a,b| b.1.partial_cmp(&a.1).unwrap());
    //             let highest_fsv = optimal_edges_with_fsv[0].1;
    //             for (edge_id, fsv) in optimal_edges_with_fsv{
    //                 let edge = &self.edges[edge_id].as_ref().unwrap();
    //                 if highest_fsv - fsv > 0.002{
    //                     edges_to_remove.lock().unwrap().insert(edge_id);
    //                     log::trace!("NON OPTIMAL ({}-{}) snp_diff:{} snp_share:{} ol:{}", edge.node1(), edge.node2(), edge.diff_snpmers, edge.shared_snpmers, edge.overlap_len_bases);
    //                 }
    //             }
    //         }
    //     });

    //     self.remove_edges(edges_to_remove.into_inner().unwrap());
    // }

    pub fn other_node(&self, node: NodeIndex, edge: &ReadOverlapEdgeTwin) -> NodeIndex {
        if edge.node1 == node {
            edge.node2
        } else {
            edge.node1
        }
    }

    pub fn prune_low_minimizer_overlaps(&mut self, c: usize, absolute_ratio: f64, relative_ratio: f64){
        let mut edges_to_remove = FxHashSet::default();
        for i in 0..self.edges.len(){
            let edge_opt = &self.edges[i];
            if edge_opt.is_none(){
                continue;
            }

            let edge = edge_opt.as_ref().unwrap();
            let num_mini = edge.shared_minimizers;
            let mini_overlap_fraction = edge.overlap_len_bases as f64 / num_mini as f64;
            if mini_overlap_fraction / c as f64 > absolute_ratio{
                edges_to_remove.insert(i);
                continue;
            }

            let node1 = self.nodes.get(&edge.node1).unwrap();
            let node2 = self.nodes.get(&edge.node2).unwrap();
            let dir1 = edge.node_edge_direction(&edge.node1);
            let dir2 = edge.node_edge_direction(&edge.node2);
            let n1_edges = node1.edges_direction(&dir1);
            let n2_edges = node2.edges_direction(&dir2);

            let best_overlap_fraction_n1 = n1_edges.iter()
                .filter_map(|edge_id| self.edges[*edge_id].as_ref())
                .map(|edge| edge.overlap_len_bases as f64 / edge.shared_minimizers as f64)
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap_or(0.0);
            let best_overlap_fraction_n2 = n2_edges.iter()
                .filter_map(|edge_id| self.edges[*edge_id].as_ref())
                .map(|edge| edge.overlap_len_bases as f64 / edge.shared_minimizers as f64)
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap_or(0.0);

            let best_overlap_fraction = best_overlap_fraction_n1.max(best_overlap_fraction_n2);
            if  mini_overlap_fraction as f64 / best_overlap_fraction > relative_ratio{
                edges_to_remove.insert(i);
                continue;
            }
        }
        log::trace!("Pruning {} low minimizer overlaps", edges_to_remove.len());

        for edge in edges_to_remove.iter(){
            let edge = self.edges[*edge].as_ref().unwrap();
            log::trace!("Removing edge {}-{} overlap len {} shared mini {}", edge.node1, edge.node2, edge.overlap_len_bases, edge.shared_minimizers);
        }

        self.remove_edges(edges_to_remove);
    }

    //Initial construction allowed lax overlaps, now prune if 
    //it doesn't cause tipping. Can look into more sophisticated heuristics later.
    pub fn prune_lax_overlaps(&mut self, c: usize, twin_reads: Option<&[TwinRead]>, snpmer_threshold_strict: f64, snpmer_error_rate_strict: f64, disable_rescue: bool){
        let rescue = true;
        let mut edges_to_remove = FxHashSet::default();
        for i in 0..self.edges.len(){
            let edge_opt = &self.edges[i];
            if edge_opt.is_none(){
                continue;
            }

            let edge = edge_opt.as_ref().unwrap();
            let mut snpmer_threshold = snpmer_threshold_strict;
            if let Some(twin_reads) = twin_reads{
                let read1_snpmer_id_thresh = twin_reads[edge.node1].snpmer_id_threshold.unwrap();
                let read2_snpmer_id_thresh = twin_reads[edge.node2].snpmer_id_threshold.unwrap();
                //snpmer_threshold = read1_snpmer_id_thresh.min(read2_snpmer_id_thresh);
                // TODO try more selective... 
                snpmer_threshold = read1_snpmer_id_thresh.max(read2_snpmer_id_thresh);
            }
            // Good overlap
            if same_strain_edge(edge, c, snpmer_threshold, snpmer_error_rate_strict){
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
            let mut short_relative_to_good1= false;
            let mut short_relative_to_good2= false;

            // TODO Get max overlap both sides, if this edge has length > 1.5x the max overlap then retain if snpmer_threshold > 99.9, diff <= 2. Only for nanopore? 

            for edge_id in edge_nd1 {
                if let Some(edge_new) = &self.edges[*edge_id]{
                    if same_strain_edge(edge_new, c, snpmer_threshold, snpmer_error_rate_strict){
                        if edge.overlap_len_bases < edge_new.overlap_len_bases * 3 / 2{
                            short_relative_to_good1 = true;
                        }
                        node1_good_found = true;
                    }
                }
            }

            let mut node2_good_found = false;
            for edge_id in edge_nd2 {
                if let Some(edge_new) = &self.edges[*edge_id]{
                    if same_strain_edge(edge_new, c, snpmer_threshold, snpmer_error_rate_strict){
                        if edge.overlap_len_bases < edge_new.overlap_len_bases * 3 / 2 {
                            short_relative_to_good2 = true;
                        }
                        node2_good_found = true;
                    }
                }
            }

            //SWITCH TODO
            // Require both nodes to have good overlaps in order to remove the edge.
            // In other words, require

            //   -->x          
            //   o ----> z     o ----> z
            //       y -->         y -->

            // Remove o->z in first case. Do not remove o->z in second. 

            if (node1_good_found && node2_good_found) || !rescue {
                if disable_rescue{
                    edges_to_remove.insert(i);
                    self.edges[i] = None;
                }
                else{
                    if short_relative_to_good1 && short_relative_to_good2{
                        edges_to_remove.insert(i);
                        self.edges[i] = None;
                    }
                    else if !same_strain_edge(edge, c, snpmer_threshold - 0.05, snpmer_error_rate_strict) || edge.diff_snpmers > 2{
                        edges_to_remove.insert(i);
                        self.edges[i] = None;
                    }
                    else{
                        log::trace!("Rescuing edge {}-{} diff snpmers:{} snpmer threshold:{} fsv:{} ol_len:{}", edge.node1, edge.node2, edge.diff_snpmers, snpmer_threshold, edge.edge_id_est(c), edge.overlap_len_bases);
                    }
                }
            }
        }

        log::trace!("Pruning {} lax overlaps", edges_to_remove.len());
        self.remove_edges(edges_to_remove);

        // let mut singleton_nodes_to_remove = vec![];
        // for (i, node) in self.nodes.iter(){
        //     if node.both_edges().collect::<Vec<_>>().is_empty(){
        //         log::trace!("Node {} has no edges after lax cuts", node.index);
        //         //singleton_nodes_to_remove.push(*i);
        //     }
        // }

        // self.remove_nodes(&singleton_nodes_to_remove, false);
    }   

    pub fn transitive_reduction(&mut self) {
        const FUZZ: usize = 750;

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
            let mut mark_direction_map = FxHashMap::default();
            for direction in [Direction::Incoming, Direction::Outgoing]{
                let sorted_all_edges = self.get_edges_sorted_by_length(*node_id, &direction);

                //Need this for singleton nodes; happens usually due to lax cuts
                if sorted_all_edges.is_empty() {
                    continue;
                }

                for &(_, edge_id) in sorted_all_edges.iter() {
                    let edge = self.edges[edge_id].as_ref().unwrap();
                    let other_node = self.other_node(*node_id, edge);
                    mark.insert(other_node, Mark::InPlay);
                }

                let longest_direction = sorted_all_edges.last().unwrap().0;

                // Iterating over the present node's edges
                //
                //            edge_id
                // node_id: o>----->o : other_node
                for &(length, edge_id) in sorted_all_edges.iter() {
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
                    for &(length2, edge_id2) in &self.get_edges_sorted_by_length(other_node, &present_edge_direction_other.reverse()) {

                        let next_edge = self.edges[edge_id2].as_ref().unwrap();
                        let third_node = self.other_node(other_node, next_edge);
                        let next_third_direction = next_edge.node_edge_direction(&third_node);

                        if mark.get(&third_node).unwrap() == &Mark::InPlay {
                            let lensum = if length + length2 < self.nodes[&other_node].base_length {0} else {length + length2 - self.nodes[&other_node].base_length};
                            //let lensum = length + length2 - self.nodes[&other_node].base_length;
                            if lensum <= longest_direction + FUZZ
                            {
                                if let Some(true) =
                                    mark.get(&third_node).map(|m| *m == Mark::InPlay)
                                {
                                    mark.insert(third_node, Mark::Eliminated);
                                    let valid_directions = mark_direction_map.entry(third_node).or_insert(FxHashSet::default());
                                    valid_directions.insert((direction, next_third_direction));
                                    log::trace!("Potential reduction from {} to {}, length1 {}, length2 {}, longets {}, lensum {}, edge_info1 {:?}, edge_info2 {:?}", node_id, third_node, length, length2, longest_direction, lensum, &edge, &next_edge);
                                }
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

    fn get_edges_sorted_by_length(&self, node_id: NodeIndex, direction: &Direction) -> Vec<(usize, EdgeIndex)> {
        let mut sorted_all_edges = vec![];
        let node_data = self.nodes.get(&node_id).unwrap();
        for (_, edge_ind) in node_data.edges_direction(&direction).into_iter().enumerate()
        {
            let edge = self.edges[*edge_ind].as_ref().unwrap();
            let n1 = self.nodes.get(&edge.node1).unwrap();
            let n2 = self.nodes.get(&edge.node2).unwrap();
            let string_length = self.nodes[&n1.index].base_length
                + self.nodes[&n2.index].base_length
                - edge.overlap_len_bases - edge.hang1 - edge.hang2;
            sorted_all_edges.push((string_length, *edge_ind));
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

pub fn get_overlaps_outer_reads_twin(twin_reads: &[TwinRead], outer_read_indices: &[usize], args: &Cli, overlap_file_path: Option<&PathBuf>) -> Vec<OverlapConfig>{

    let bufwriter = if let Some(overlap_file_path) = overlap_file_path{
        let file = File::create(overlap_file_path).unwrap();
        Mutex::new(Box::new(GzEncoder::new(BufWriter::new(file), Compression::default())) as Box<dyn Write + Send>)
        //Mutex::new(BufWriter::new(std::fs::File::create(overlaps_file).unwrap()))
    }
    else{
        Mutex::new(Box::new(std::io::sink()) as Box<dyn Write + Send>)
    };

    let outer_twin_reads = outer_read_indices.iter().map(|x| (*x, &twin_reads[*x])).collect::<FxHashMap<_,_>>();

    // Contained reads may still lurk in the outer reads, remove again for sure.
    let contained_reads_again = Mutex::new(FxHashSet::default());

    let inverted_index = get_minimizer_index(None, Some(&outer_twin_reads));
    
    let overlaps = Mutex::new(vec![]);


    outer_read_indices.into_par_iter().for_each(|&i| { 
        let read = &twin_reads[i];
        let mini_anchors = find_exact_matches_with_full_index(&read.minimizers_vec(), &inverted_index, None, Some(&outer_twin_reads));

        let comparison_options = CompareTwinReadOptions{
            compare_snpmers: true,
            retain_chain: false,
            force_one_to_one_alignments: true,
            supplementary_threshold_score: Some(500.),
            supplementary_threshold_ratio: Some(0.25),
            secondary_threshold: None,
            read1_mininimizers: None, // indexed anchors
            read1_snpmers: Some(read.snpmers_vec()),
        };

        for (outer_ref_id, anchors) in mini_anchors.into_iter(){

            //Only compare once. I think we get slightly different results if we
            //compare in both directinos, but this forces consistency. 
            if i <= outer_ref_id as usize {
                continue;
            }

            let read2 = &twin_reads[outer_ref_id as usize];

            if !dovetail_possibility(&anchors, &read, &read2){
                continue;
            }

            let twlaps = compare_twin_reads(&read, &read2, Some(&anchors), None, i, outer_ref_id as usize, &comparison_options, args);

            if twlaps.len() > 1{
                log::trace!("Multiple overlaps for {}:{} and {}:{}", &read.id, i,  &read2.id, outer_ref_id);
                for twlap in twlaps.iter(){
                    log::trace!("{}-{} Overlap: {}-{} {}-{}, reverse {}", i, outer_ref_id, twlap.start1, twlap.end1, twlap.start2, twlap.end2, twlap.chain_reverse);
                }
            }

            let mut possible_containment = false;
            //Check for contained read
            for twlap in twlaps.iter(){
                let mut twlap_contain = false;
                let mut smaller_read_index = i;
                let mut snpmer_threshold = args.snpmer_threshold_strict;
                if r1_contained_r2(&twlap, &read, &read2, true, args.c){
                    twlap_contain = true;
                    possible_containment = true;
                    snpmer_threshold = read.snpmer_id_threshold.unwrap_or(100.);
                }
                else if r1_contained_r2(&twlap, &read2, &read, true, args.c){
                    twlap_contain = true;
                    possible_containment = true;
                    smaller_read_index = outer_ref_id as usize;
                    snpmer_threshold = read2.snpmer_id_threshold.unwrap_or(100.);
                }
                if twlap_contain {
                    if same_strain(twlap.shared_minimizers, twlap.diff_snpmers, twlap.shared_snpmers, args.c as u64, snpmer_threshold, args.snpmer_error_rate_strict, twlap.large_indel){
                        contained_reads_again.lock().unwrap().insert(smaller_read_index);
                        log::trace!("Contained read {} in read {} STRICT", i, outer_ref_id);
                    }
                }
            }

            if possible_containment {
                if twlaps.len() == 1 {
                    let twlap = &twlaps[0];
                    let contained_overlap_config = OverlapConfig{
                        hang1: 0,
                        hang2: 0,
                        //Forward doesn't make sense for contained reads
                        forward1: true,
                        forward2: true,
                        overlap1_len: twlap.end1 - twlap.start1,
                        overlap2_len: twlap.end2 - twlap.start2,
                        read_i: twlap.i1,
                        read_j: twlap.i2,
                        shared_mini: twlap.shared_minimizers,
                        shared_snpmer: twlap.shared_snpmers,
                        diff_snpmer: twlap.diff_snpmers,
                        contained: true,
                        large_indel: twlap.large_indel
                    };
                    log::trace!("Contained read {} in read {}", twlap.i2, twlap.i1);
                    overlaps.lock().unwrap().push(contained_overlap_config);
                }
                continue
            }

            let best_overlaps = comparison_to_overlap(twlaps, &twin_reads, args, &bufwriter);
            for best_ol in best_overlaps{
                overlaps.lock().unwrap().push(best_ol);
            }
        }
    });

    let contained_reads = contained_reads_again.into_inner().unwrap();
    let mut ol = overlaps.into_inner().unwrap();
    ol.retain(|x| !contained_reads.contains(&x.read_i) && !contained_reads.contains(&x.read_j));

    return ol;
}

fn overlap_hang_length(twin_read: &TwinRead) -> (usize, usize){
    let kmer_error_est = 1./(twin_read.est_id.unwrap_or(100.0) / 100.).powf(twin_read.k as f64);
    let nth_read = ((MINIMIZER_END_NTH_OVERLAP as f64 * kmer_error_est) as usize).min(100);
    let (start_hang_cutoff, end_hang_cutoff) = map_processing::first_last_mini_in_range(0, twin_read.base_length, twin_read.k  as usize, nth_read,  &twin_read.minimizer_positions);
    let return_start = (start_hang_cutoff).min(OVERLAP_HANG_LENGTH);
    let return_end = (twin_read.base_length - end_hang_cutoff).min(OVERLAP_HANG_LENGTH);
    return (return_start, return_end);
}

pub fn comparison_to_overlap<T>(twlaps: Vec<TwinOverlap>, twin_reads: &[TwinRead], args: &Cli, writer: &Mutex<T>) -> Vec<OverlapConfig>
where T : Write + Send
{

    let mut overlap_possib_out = vec![];
    let mut overlap_possib_in = vec![];

    for twlap in twlaps{
        let overlap_possibility = overlap_config_from_twlap(&twlap, twin_reads, args, writer);
        if let Some(ol) = overlap_possibility{
            let first_index = ol.read_i == twlap.i1;
            if first_index && ol.forward1{
                overlap_possib_out.push(ol);
            } 
            else if first_index && !ol.forward1{
                overlap_possib_in.push(ol);
            }
            else if !first_index && ol.forward2{
                overlap_possib_in.push(ol);
            }
            else if !first_index && !ol.forward2{
                overlap_possib_out.push(ol);
            }
        }
    }
    
    let best_overlap_forward = overlap_possib_out.into_iter().max_by_key(|x| (x.overlap1_len + x.overlap2_len) as i32 - (x.hang1 + x.hang2) as i32);
    let best_overlap_backward = overlap_possib_in.into_iter().max_by_key(|x| (x.overlap1_len + x.overlap2_len) as i32 - (x.hang1 + x.hang2) as i32);
    let mut best_overlaps = vec![];

    if best_overlap_forward.is_some() && best_overlap_backward.is_none(){
        best_overlaps.push(best_overlap_forward.unwrap());
    }
    else if best_overlap_forward.is_none() && best_overlap_backward.is_some(){
        best_overlaps.push(best_overlap_backward.unwrap());
    }
    else if best_overlap_forward.is_some() && best_overlap_backward.is_some(){

        // Check if circular overlap is concordant.
        let ol_f = best_overlap_forward.unwrap();
        let ol_b = best_overlap_backward.unwrap();

        let reverse_f = ol_f.forward1 != ol_f.forward2;
        let reverse_b = ol_b.forward1 != ol_b.forward2;

        //Consistent
        if reverse_f == reverse_b{
            best_overlaps.push(ol_f);
            best_overlaps.push(ol_b);
        }
        //Inconsistent; take best
        else{
            let best_ol = [ol_f, ol_b].into_iter().max_by_key(|x| (x.overlap1_len + x.overlap2_len) as i32 - (x.hang1 + x.hang2) as i32).unwrap();
            best_overlaps.push(best_ol);
        }
    }

    return best_overlaps;
}

fn overlap_config_from_twlap<T>(twlap: &TwinOverlap, twin_reads: &[TwinRead], args: &Cli, writer: &Mutex<T>) -> Option<OverlapConfig>
where T: Write + Send
{

    let mut overlap_possibilities = vec![];

    let i = twlap.i1;
    let j = twlap.i2;
    let mut exist_overlap = false;

    let read1 = &twin_reads[i];
    let read2 = &twin_reads[j];

    //check if end-to-end overlap
    let identity = id_est(twlap.shared_minimizers, twlap.diff_snpmers, args.c as u64, twlap.large_indel);
    let same_strain_lax = same_strain(twlap.shared_minimizers, twlap.diff_snpmers, twlap.shared_snpmers, args.c as u64, args.snpmer_threshold_lax, args.snpmer_error_rate_lax, twlap.large_indel);

    let aln_len1 = twlap.end1 - twlap.start1;
    let aln_len2 = twlap.end2 - twlap.start2;

    if aln_len1.max(aln_len2) < args.min_ol{
        return None;
    }

    let (hang1_start, hang1_end) = overlap_hang_length(read1);
    let (hang2_start, hang2_end) = overlap_hang_length(read2);

    let hang_start = hang1_start.max(hang2_start);
    let hang_end = hang1_end.max(hang2_end);

    //let (ext_s1, ext_e1, ext_s2, ext_e2) = alignment::extend_ends_chain(&twin_reads[twlap.i1].dna_seq, &twin_reads[twlap.i2].dna_seq, twlap, args);

    //let aln_len1 = ext_e1 - ext_s1 + 1;
    //let aln_len2 = ext_e2 - ext_s2 + 1;

    //let mini_chain = twlap.minimizer_chain.as_ref().unwrap();

    // println!("EXT {}, {}, {}, {} READ {}", ext_s1, ext_e1, ext_s2, ext_e2, &read1.id);
    // println!("OLRANGE {}, {}, {}, {} READ {}", twlap.start1, twlap.end1, twlap.start2, twlap.end2, &read1.id);
    // println!("HANG {}, {}, {}, {} READ {}", hang1_start, hang1_end, hang2_start, hang2_end, &read1.id);

    // let (hang1_start, hang1_end) = (OVERLAP_HANG_LENGTH, OVERLAP_HANG_LENGTH);
    // let (hang2_start, hang2_end) = (OVERLAP_HANG_LENGTH, OVERLAP_HANG_LENGTH);
    
    if twlap.chain_reverse {
        if twlap.start1 < hang_start && twlap.start2 < hang_start {
            let ol_config = OverlapConfig {
                hang1: twlap.start1,
                hang2: twlap.start2,
                forward1: false,
                forward2: true,
                overlap1_len: aln_len1,
                overlap2_len: aln_len2,
                read_i: i,
                read_j: j, 
                shared_mini: twlap.shared_minimizers,
                shared_snpmer: twlap.shared_snpmers,
                diff_snpmer: twlap.diff_snpmers,
                contained: false,
                large_indel: twlap.large_indel
            };
            if same_strain_lax{
                overlap_possibilities.push(ol_config);
            }
            exist_overlap = true;
        } 
        else if twlap.end1 + hang_end > read1.base_length && twlap.end2 + hang_end > read2.base_length 
        {
            let ol_config = OverlapConfig {
                //hang1: read1.base_length - twlap.end1 - 1,
                //hang2: read2.base_length - twlap.end2 - 1,
                hang1: read1.base_length - twlap.end1 - 1,
                hang2: read2.base_length - twlap.end2 - 1,
                forward1: true,
                forward2: false,
                overlap1_len: aln_len1,
                overlap2_len: aln_len2,
                read_i: i,
                read_j: j,
                shared_mini: twlap.shared_minimizers,
                shared_snpmer: twlap.shared_snpmers,
                diff_snpmer: twlap.diff_snpmers,
                contained: false,
                large_indel: twlap.large_indel
            };
            if same_strain_lax{
                overlap_possibilities.push(ol_config);
            }
            exist_overlap = true;
        }
    } else {
        if twlap.start1 < hang_start && twlap.end2 + hang_end > read2.base_length {
            let ol_config = OverlapConfig {
                //hang1: twlap.start1,
                //hang2: read2.base_length - twlap.end2 - 1,
                hang1: twlap.start1,
                hang2: read2.base_length - twlap.end2 - 1,
                forward1: false,
                forward2: false,
                overlap1_len: aln_len1,
                overlap2_len: aln_len2,
                read_i: i,
                read_j: j,
                shared_mini: twlap.shared_minimizers,
                shared_snpmer: twlap.shared_snpmers,
                diff_snpmer: twlap.diff_snpmers,
                contained: false,
                large_indel: twlap.large_indel
            };
            if same_strain_lax{
                overlap_possibilities.push(ol_config);
            }
            exist_overlap = true;
        } 
        else if twlap.end1 + hang_end > read1.base_length  && twlap.start2 < hang_start{
            let ol_config = OverlapConfig {
                //hang1: read1.base_length - twlap.end1 - 1,
                //hang2: twlap.start2,
                hang1: read1.base_length - twlap.end1 - 1,
                hang2: twlap.start2,
                forward1: true,
                forward2: true,
                overlap1_len: aln_len1,
                overlap2_len: aln_len2,
                read_i: i,
                read_j: j,
                shared_mini: twlap.shared_minimizers,
                shared_snpmer: twlap.shared_snpmers,
                diff_snpmer: twlap.diff_snpmers,
                contained: false,
                large_indel: twlap.large_indel
            };
            if same_strain_lax{
                overlap_possibilities.push(ol_config);
            }
            exist_overlap = true;
        }
    }

    if exist_overlap {
        let mut bufwriter = writer.lock().unwrap();
        let mut possibilties_string = String::new();
        for possib in overlap_possibilities.iter(){
            possibilties_string.push_str(&format!("{}:{}-{}:{} HANG {} {}", possib.read_i, possib.forward1, possib.read_j, possib.forward2, possib.hang1, possib.hang2));
        }
        writeln!(bufwriter,
            "{} {} {} {} fsv:{} SHARE:{} DIFF:{} MINI: {}, LEN1:{} {}-{} LEN2:{} {}-{}, REVERSE: {}, Possibilties: {}",
            i,
            j,
            &read1.id.split_ascii_whitespace().next().unwrap(),
            &read2.id.split_ascii_whitespace().next().unwrap(),
            identity * 100.,
            twlap.shared_snpmers,
            twlap.diff_snpmers,
            twlap.shared_minimizers,
            read1.base_length,
            twlap.start1,
            twlap.end1,
            read2.base_length,
            twlap.start2,
            twlap.end2,
            twlap.chain_reverse,
            possibilties_string
        ).unwrap();
    }

    let best_overlap = overlap_possibilities.into_iter().max_by_key(|x| (x.overlap1_len as i64 + x.overlap2_len as i64) - (x.hang1 as i64 + x.hang2 as i64) - (x.hang1 as i64 - x.hang2 as i64).abs());
    return best_overlap;
}

pub fn read_graph_from_overlaps_twin(overlaps: Vec<OverlapConfig>, twin_reads: &[TwinRead], outer_reads: Option<&Vec<usize>>, args: &Cli) -> OverlapTwinGraph
{
    let mut nodes = FxHashMap::default();
    let mut edges = vec![];

    for overlap in overlaps{
        if overlap.contained{
            continue
        }

        let i = overlap.read_i;
        let j = overlap.read_j;
        let new_read_overlap = ReadOverlapEdgeTwin {
            node1: overlap.read_i,
            node2: overlap.read_j,
            hang1: overlap.hang1,
            hang2: overlap.hang2,
            overlap1_len: overlap.overlap1_len,
            overlap2_len: overlap.overlap2_len,
            forward1: overlap.forward1,
            forward2: overlap.forward2,
            overlap_len_bases: overlap.overlap1_len.max(overlap.overlap2_len),
            shared_minimizers: overlap.shared_mini,
            diff_snpmers: overlap.diff_snpmer,
            shared_snpmers: overlap.shared_snpmer,
            large_indel: false
        };

        edges.push(Some(new_read_overlap));
        let ind = edges.len() - 1;
        {
            let rd1 = nodes.entry(i).or_insert(ReadData::default());
            rd1.index = i;
            if overlap.forward1 {
                rd1.out_edges.push(ind)
            } else {
                rd1.in_edges.push(ind)
            }
            rd1.base_length = twin_reads[i].base_length;
            rd1.read_id = twin_reads[i].id.clone();
        }
        {
            let rd2 = nodes.entry(j).or_insert(ReadData::default());
            rd2.index = j;
            if overlap.forward2 {
                rd2.in_edges.push(ind)
            } else {
                rd2.out_edges.push(ind)
            }
            rd2.base_length = twin_reads[j].base_length;
            rd2.read_id = twin_reads[j].id.clone();
        }
    }

    if let Some(outer_read_indices) = outer_reads{
        let mut num_singletons = 0;
        let mut used_reads = FxHashSet::default();
        for read in nodes.values(){
            used_reads.insert(read.index);
        }

        for &i in outer_read_indices.iter(){
            if !used_reads.contains(&i){
                // Small reads are likely to be singletons cause they can't overlap anyways
                if twin_reads[i].base_length < args.min_ol * 3 / 2{
                    continue;
                }
                let rd = nodes.entry(i).or_insert(ReadData::default());
                rd.index = i;
                rd.base_length = twin_reads[i].base_length;
                rd.read_id = twin_reads[i].id.clone();
                log::trace!("Read {} {} has no overlaps. Adding singleton to unitig graph", i, rd.read_id);
                num_singletons += 1;
            }
        }

        log::debug!("Number of singleton reads: {}", num_singletons);
    }

    let mut graph = OverlapTwinGraph { nodes, edges };

    //graph.get_maximal_overlaps(args);
    graph.prune_low_minimizer_overlaps(args.c, args.absolute_minimizer_cut_ratio, args.relative_minimizer_cut_ratio);
    graph.prune_lax_overlaps(args.c, Some(twin_reads), args.snpmer_threshold_strict, args.snpmer_error_rate_strict, args.disable_error_overlap_rescue);
    graph.transitive_reduction();

    log::trace!("Number of nodes in initial read overlap graph: {}", graph.nodes.len());

    return graph;
}

pub fn remove_contained_reads_twin(query_indices: Option<Vec<usize>>, ref_indices: Option<Vec<usize>>, twin_reads: &[TwinRead], first_iteration: bool, temp_dir: &PathBuf,  args: &Cli) -> Vec<usize>{
    let start = std::time::Instant::now();
    let downsample_factor = (args.contain_subsample_rate / args.c).max(1) as u64;
    let reads_to_index;
    if let Some(ref_indices) = ref_indices{
        reads_to_index = ref_indices.iter().map(|x| *x).collect::<FxHashSet<_>>();
    }
    else{
        reads_to_index = (0..twin_reads.len()).collect::<FxHashSet<_>>();
    }
    let mut inverted_index_hashmap =
        twin_reads
            .iter()
            .enumerate()
            .filter(
                |x| 
                reads_to_index.contains(&x.0) &&
                (x.1.est_id.is_none()
                 || x.1.est_id.unwrap() > args.quality_value_cutoff)
            )
            .fold(FxHashMap::default(), |mut acc, (i, x)| {
                for y in x.minimizer_kmers().iter(){
                    if hash64(y) > u64::MAX / downsample_factor {
                        continue;
                    }
                    acc.entry(*y).or_insert(FxHashSet::default()).insert(i as u32);
                }
                acc
            });
    
    let mut kmer_to_count = inverted_index_hashmap.iter().map(|(_,v)| v.len()).collect::<Vec<_>>();
    kmer_to_count.sort_by(|a,b| b.cmp(&a));
    let threshold = kmer_to_count[kmer_to_count.len() / 100_000];
    inverted_index_hashmap.retain(|_,v| v.len() < threshold);
    log::debug!("Number of kmer indices in inverted index: {}. Threshold: {}", inverted_index_hashmap.len(), threshold);
    //println!("Time to build inverted index hashmap: {:?}", start.elapsed());

    //open file for writing
    let name = if query_indices.is_none() { "all-cont.txt.gz" } else { "subset-cont.txt.gz" };
    //let output_path = Path::new(args.output_dir.as_str()).join(name);
    let output_path = temp_dir.join(name);
    let bufwriter_dbg;
    if first_iteration && !log::log_enabled!(log::Level::Trace) {
        //write to null
        bufwriter_dbg = Mutex::new(Box::new(std::io::sink()) as Box<dyn Write + Send>);
    }
    else{
        let writer = BufWriter::new(File::create(&output_path).unwrap());
        bufwriter_dbg = Mutex::new(Box::new(GzEncoder::new(writer, Compression::default())) as Box<dyn Write + Send>);
    }
    let contained_reads = Mutex::new(FxHashSet::default());
    let outer_reads = Mutex::new(vec![]);

    let query_range = if let Some(special_indices) = query_indices{
        special_indices
    } else {
        // ALl reads get queried
        (0..twin_reads.len()).into_iter().collect::<Vec<_>>()
    };

    parallel_remove_contained(query_range, twin_reads, &inverted_index_hashmap, &bufwriter_dbg, &contained_reads, &outer_reads, downsample_factor, args);
    let num_contained_reads = contained_reads.lock().unwrap().len();
    log::info!("{} reads are contained; {} outer reads", num_contained_reads, outer_reads.lock().unwrap().len());

    log::info!("Time elapsed for removing contained reads is: {:?}", start.elapsed());
    return outer_reads.into_inner().unwrap();
}

fn parallel_remove_contained<T>(
    range: Vec<usize>,
    twin_reads: &[TwinRead],
    inverted_index_hashmap: &FxHashMap<Kmer48, FxHashSet<u32>>,
    bufwriter_dbg: &Mutex<T>,
    contained_reads: &Mutex<FxHashSet<usize>>,
    outer_reads: &Mutex<Vec<usize>>,
    downsample_factor: u64,
    args: &Cli,
) where T : Write + Send
{
    range.into_par_iter().for_each(|i| {
        let mut contained = false;
        let read1 = &twin_reads[i];
        if read1.est_id.is_some() && read1.est_id.unwrap() < args.quality_value_cutoff {
            return;
        }

        let mut comparison_options = CompareTwinReadOptions::default();
        comparison_options.read1_mininimizers = Some(read1.minimizers_vec());
        comparison_options.read1_snpmers = Some(read1.snpmers_vec());

        let start = std::time::Instant::now();
        let mut index_count_map = FxHashMap::default();
        let mut index_range_map = FxHashMap::default();

        for y in comparison_options.read1_mininimizers.as_ref().unwrap().iter() {
            //Downsample to 100 compression factor
            if hash64(&y.1) > u64::MAX / downsample_factor{
                continue;
            }
            if let Some(indices) = inverted_index_hashmap.get(&y.1) {
                for &index in indices {
                    if index == i as u32 {
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

        top_indices.retain(|(index,_)| twin_reads[(**index) as usize].base_length > read1.base_length);
        //top_indices.retain(|(_,count)| **count > read1.minimizers.len() as u32 / 10 && **count > 5);
        top_indices.retain(|(index,_)| {
            let range = index_range_map.get(&index).unwrap();
            range[1] > range[0] && (range[1] - range[0]) as f64 + (20. * args.contain_subsample_rate as f64) > read1.base_length as f64 * 0.90
        });

        top_indices.truncate(500);

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
        if log::log_enabled!(log::Level::Trace) {
            writeln!(bufwriter_dbg.lock().unwrap(), "CONTAIN: {} NUMBER OF HITS {}. THRESHOLD: {}", &read1.id, top_indices.len(), read1.snpmer_id_threshold.unwrap_or(100.)).unwrap();
        }

        let start = std::time::Instant::now();
        let num_tries = 50;
        let mut num_fails = 0;
        let mut max_ol = 0;
        for (index, order_count) in top_indices.into_iter() {
            
            if contained{
                break;
            }

            let read2 = &twin_reads[(*index) as usize];
            let twin_overlaps = compare_twin_reads(&read1, &read2, None, None, i, (*index) as usize, &comparison_options, args);
            if twin_overlaps.is_empty(){
                continue;
            }
            let twin_overlap = twin_overlaps.first().unwrap();
            let shared_minimizers = twin_overlap.shared_minimizers;
            let diff_snpmers = twin_overlap.diff_snpmers;

            let identity = id_est(shared_minimizers, diff_snpmers, args.c as u64, twin_overlap.large_indel);
            let snpmer_threshold = read1.snpmer_id_threshold.unwrap_or(args.snpmer_threshold_strict);
            let same_strain = same_strain(shared_minimizers, diff_snpmers, twin_overlap.shared_snpmers, args.c as u64, snpmer_threshold, args.snpmer_error_rate_strict, twin_overlap.large_indel);

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
            if log::log_enabled!(log::Level::Trace) {
                writeln!(
                    bufwriter_dbg.lock().unwrap(),
                    "{} {}:{}-{} ----- {} {}:{}-{}   minis: {} shared_snps: {}, diff_snps: {}, identity {}, ol_len {}, read1len: {}, read2len: {}, order_count: {}, snpmer_thresh: {}",
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
                    order_count,
                    snpmer_threshold

                )
                .unwrap();
            }


            //Can't be reflexive.
            if r1_contained_r2(&twin_overlap, read1, read2, same_strain, args.c){
                let comparing_time = start.elapsed();
                if log::log_enabled!(log::Level::Trace) {
                    writeln!(bufwriter_dbg.lock().unwrap(), "{} {} CONTAINED. Times {} {} {}", read1.id, read2.id, inverted_indexing_time.as_micros(), sorting_filtering_time.as_micros(), comparing_time.as_micros()).unwrap();
                }
                else{
                    writeln!(bufwriter_dbg.lock().unwrap(), "{} -> {}", read1.id, read2.id).unwrap();
                }
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

pub fn same_strain_edge(edge: &ReadOverlapEdgeTwin, c: usize, snpmer_threshold: f64, snpmer_error_rate: f64) -> bool {
    return same_strain(edge.shared_minimizers, edge.diff_snpmers, edge.shared_snpmers, c as u64, snpmer_threshold, snpmer_error_rate, edge.large_indel);
}

pub fn same_strain(minimizers: usize, snp_diff: usize, snp_shared: usize,  c: u64, snpmer_threshold: f64, snpmer_error_rate: f64, large_indel: bool) -> bool {
    assert!(snpmer_threshold > 1.0); // Some ambiguity with percentages vs fractions...
    let identity = id_est(minimizers, snp_diff, c, large_indel);
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

#[inline]
fn dovetail_possibility(anchors: &Anchors, read1: &TwinRead, read2: &TwinRead) -> bool {

    let read1_length = read1.base_length;
    let read2_length = read2.base_length;

    let mut read1_possible = false;
    let mut read2_possible = false;

    //TODO change 750 to an adaptive threshold based on solid k-mers and error rates?
    for anchor in anchors.anchors.iter(){
        if anchor.pos1 < OVERLAP_HANG_LENGTH as u32 || anchor.pos1 + OVERLAP_HANG_LENGTH as u32 > read1_length as u32{
            read1_possible = true;
        }
        if anchor.pos2 < OVERLAP_HANG_LENGTH as u32 || anchor.pos2 + OVERLAP_HANG_LENGTH as u32 > read2_length as u32{
            read2_possible = true;
        }

        if read2_possible && read1_possible{
            return true;
        }
    }

    return false;
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
        pub overlap_len_bases: usize,
        pub shared_minimizers: usize,
    }
    
    impl MockEdge{
        pub fn new(i: usize, j: usize, d1: Direction, d2: Direction) -> Self{
            Self{
                i,
                j,
                d1,
                d2,
                diff_snpmers: 0,
                overlap_len_bases: 2000,
                shared_minimizers: 200,
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
                overlap1_len: edge.overlap_len_bases,
                overlap2_len: edge.overlap_len_bases,
                forward1: edge.d1 == Direction::Outgoing,
                forward2: edge.d2 == Direction::Incoming,
                overlap_len_bases: edge.overlap_len_bases,
                shared_minimizers: edge.shared_minimizers,
                diff_snpmers: edge.diff_snpmers,
                shared_snpmers: 10,
                large_indel: false,
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

    fn overlap_lengths(graph: &mut OverlapTwinGraph, lens : &[usize]){
        for (i,edge) in graph.edges.iter_mut().enumerate(){
            if let Some(edge) = edge{
                edge.overlap1_len = lens[i];
                edge.overlap2_len = lens[i];
                edge.overlap_len_bases = lens[i];
            }
        }
    }

    #[test]
    fn test_same_strain() {
        let minimizers = 100;
        let snp_diff = 10;
        let snp_shared = 10;
        let c = 10;
        let snpmer_threshold = 99.9;
        let snpmer_error_rate = 0.025;
        assert_eq!(same_strain(minimizers, snp_diff, snp_shared, c, snpmer_threshold, snpmer_error_rate, false), false);

        let minimizers = 100;
        let snp_diff = 1;
        let snp_shared = 0;
        assert_eq!(same_strain(minimizers, snp_diff, snp_shared, c, snpmer_threshold, snpmer_error_rate, false), true);

        let minimizers = 10000;
        let snp_diff = 5;
        let snp_shared = 0;
        assert_eq!(same_strain(minimizers, snp_diff, snp_shared, c, snpmer_threshold, snpmer_error_rate, false), true);

        let minimizers = 10;
        let snp_diff = 1;
        let snp_shared = 0;
        assert_eq!(same_strain(minimizers, snp_diff, snp_shared, c, snpmer_threshold, snpmer_error_rate, false), false);

        let minimizers = 100;
        let snp_diff = 5;
        let snp_shared = 10;
        let snpmer_error_rate = 0.50;
        assert_eq!(same_strain(minimizers, snp_diff, snp_shared, c, snpmer_threshold, snpmer_error_rate, false), true);
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

        let lens = [2000, 2000, 1000];

        let mut graph = mock_graph_from_edges(mock_edges1);
        overlap_lengths(&mut graph, &lens);
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


        let lens = [2000, 2000, 1000];

        let mut graph = mock_graph_from_edges(mock_edges1);
        overlap_lengths(&mut graph, &lens);
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

        let lens = [2000, 2000, 1000];
        let mut graph = mock_graph_from_edges(mock_edges1);
        overlap_lengths(&mut graph, &lens);
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
        let lens = [2000, 2000, 1000];

        let mut graph = mock_graph_from_edges(mock_edges1);
        overlap_lengths(&mut graph, &lens);
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

        let lens = [1000, 2000, 2000];

        let mut graph = mock_graph_from_edges(mock_edges1);
        overlap_lengths(&mut graph, &lens);
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
    fn test_transitive_reduction_small_cycle() {

        // 3 overlapping reads on a small plasmid, multiple double edges
        // -------------------> 0
        // xxxx     |||||||||||
        // ---->   ------------ 1
        // ||||     xxx   ||||||
        // ------------>  ----- 2

        let mock_edges1 = vec![
            MockEdge::new(0, 1, Outgoing, Incoming), // 2000
            MockEdge::new(1, 0, Outgoing, Incoming), // 1000
            MockEdge::new(0, 2, Outgoing, Incoming), // 1000
            MockEdge::new(2, 0, Outgoing, Incoming), // 2000
            MockEdge::new(1, 2, Outgoing, Incoming), // 2000
            MockEdge::new(2, 1, Outgoing, Incoming), // 1000
        ];

        let mut overlap_lengths = vec![2000, 1000, 1000, 2000, 2000, 1000];

        let mut graph = mock_graph_from_edges(mock_edges1);
        for edge in graph.edges.iter_mut(){
            if let Some(edge) = edge{
                edge.overlap1_len = overlap_lengths.pop().unwrap();
                edge.overlap2_len = edge.overlap1_len;
                edge.overlap_len_bases = edge.overlap1_len;
            }
        }

        dbg!(&graph);

        graph.transitive_reduction();
        let some_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        // We want exactly a triangle.
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

        let lens = 
        [
            2000,
            2000,
            2000,
            2000,
            2000,
            2000,
            2000,
            2000,
            1100,
            900,
            100,
            2000,
            800,
        ];

        let mut graph = mock_graph_from_edges(mock_edges1);
        overlap_lengths(&mut graph, &lens);
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
        graph.prune_lax_overlaps(8, None, 100., 0., true);
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
        graph.prune_lax_overlaps(8, None, 100., 0., true);
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
        graph.prune_lax_overlaps(8, None, 100., 0., true);
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
        graph.prune_lax_overlaps(8, None, 100., 0., true);
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
        graph.prune_lax_overlaps(8, None, 100., 0., true);
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
        graph.prune_lax_overlaps(8, None, 100., 0., true);
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 2);
    }

    fn create_snpmer_config(positions: Vec<(u32, bool)>, edge: usize, length: usize) -> SnpmerOverlapConfig {
        let mut config = SnpmerOverlapConfig::default();
        config.snpmer_overlap_set = positions.into_iter().collect();
        config.edge = edge;
        config.length = length;
        config
    }

    #[test]
    fn test_greater_than_basic() {
        // Test case 1: Clear superset with greater length
        let config1 = create_snpmer_config(
            vec![(1, true), (2, true), (3, false)],
            0,
            100
        );
        let config2 = create_snpmer_config(
            vec![(1, true), (2, true)],
            1,
            50
        );
        assert!(config1.greater_than(&config2));
        assert!(!config2.greater_than(&config1));
    }

    #[test]
    fn test_greater_than_same_positions() {
        // Test case where sets have same positions but different lengths
        let config1 = create_snpmer_config(
            vec![(1, true), (2, false)],
            0,
            100
        );
        let config2 = create_snpmer_config(
            vec![(1, true), (2, false)],
            1,
            50
        );
        assert!(config1.greater_than(&config2));
        assert!(!config2.greater_than(&config1));
    }

    #[test]
    fn test_greater_than_length_requirement() {
        // Test that greater_than requires strictly greater length
        let config1 = create_snpmer_config(
            vec![(1, true), (2, true), (3, false)],
            0,
            100
        );
        let config2 = create_snpmer_config(
            vec![(1, true), (2, true)],
            1,
            100
        );
        assert!(!config1.greater_than(&config2));
    }

    #[test]
    fn test_get_maximal_simple() {
        let mut configs = vec![
            create_snpmer_config(vec![(1, true), (2, true)], 0, 100),
            create_snpmer_config(vec![(1, true)], 1, 50),
            create_snpmer_config(vec![(2, true)], 2, 75)
        ];


        SnpmerOverlapConfig::sort(&mut configs);
        
        let maximal = SnpmerOverlapConfig::get_maximal(&configs);
        assert_eq!(maximal, vec![2]);
    }

    #[test]
    fn test_get_maximal_incomparable() {
        // Test case where sets are incomparable
        let mut configs = vec![
            create_snpmer_config(vec![(1, true), (2, false)], 0, 100),
            create_snpmer_config(vec![(2, true), (3, false)], 1, 100),
            create_snpmer_config(vec![(1, true), (3, true)], 2, 100)
        ];
        SnpmerOverlapConfig::sort(&mut configs);
        
        let maximal = SnpmerOverlapConfig::get_maximal(&configs);
        // All should be maximal as they're incomparable
        assert_eq!(maximal.len(), 3);
    }

    #[test]
    fn test_get_maximal_chain() {
        // Test case with a chain of strict inclusions
        let mut configs = vec![
            create_snpmer_config(vec![(1, true)], 0, 50),
            create_snpmer_config(vec![(1, true), (2, true)], 1, 75),
            create_snpmer_config(vec![(1, true), (2, true), (3, true)], 2, 100)
        ];

        SnpmerOverlapConfig::sort(&mut configs);
        
        let maximal = SnpmerOverlapConfig::get_maximal(&configs);
        assert_eq!(maximal, vec![2]);
    }

    #[test]
    fn test_get_maximal_empty() {
        let configs: Vec<SnpmerOverlapConfig> = vec![];
        let maximal = SnpmerOverlapConfig::get_maximal(&configs);
        assert!(maximal.is_empty());
    }

    #[test]
    fn test_get_maximal_single() {
        let configs = vec![
            create_snpmer_config(vec![(1, true)], 0, 50)
        ];
        
        let maximal = SnpmerOverlapConfig::get_maximal(&configs);
        assert_eq!(maximal, vec![0]);
    }

    #[test]
    fn test_get_maximal_mixed_matches() {
        // Test with mixed matching and non-matching positions
        let mut configs = vec![
            create_snpmer_config(vec![(1, true), (2, false), (3, true)], 0, 100),
            create_snpmer_config(vec![(1, true), (2, false)], 1, 75),
            create_snpmer_config(vec![(1, true), (2, false), (3, false)], 2, 90)
        ];

        SnpmerOverlapConfig::sort(&mut configs);
        dbg!(&configs);
        
        let maximal = SnpmerOverlapConfig::get_maximal(&configs);
        assert_eq!(maximal, vec![2, 1]);
    }

    #[test]
    fn test_prune_lax_overlaps_simple_rescue_heuristic_negative(){
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
        graph.prune_lax_overlaps(8, None, 100., 0., false);
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
        graph.prune_lax_overlaps(8, None, 100., 0., false);
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 2);
    }

    #[test]
    fn test_prune_lax_overlaps_simple_rescue_heuristic_positive(){
        // 1 -- > 2
        //   ↘ (x)
        // 3 --> 4
        let mut mock_edges1 = vec![
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(1, 4, Outgoing, Incoming),
            MockEdge::new(3, 4, Outgoing, Incoming),
        ];

        mock_edges1[1].diff_snpmers = 1;
        mock_edges1[1].overlap_len_bases = 5000;
        mock_edges1[1].shared_minimizers = 500;

        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.prune_lax_overlaps(8, None, 100., 0., false);
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 3);

        // Change orientations
        // 1 >--< 2
        //    ↘ (x)
        // 3 <--> 4
        let mut mock_edges1 = vec![
            MockEdge::new(1, 2, Outgoing, Outgoing),
            MockEdge::new(1, 4, Outgoing, Incoming),
            MockEdge::new(3, 4, Incoming, Incoming),
        ];

        mock_edges1[1].diff_snpmers = 1;
        mock_edges1[1].overlap_len_bases = 3000;
        mock_edges1[1].shared_minimizers = 400;


        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.prune_lax_overlaps(8, None, 100., 0., false);
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 3);
    }

    #[test]
    fn test_prune_lax_overlaps_cross_rescue(){
        // 1 -- > 2
        //   x
        // 3 --> 4
        let mut mock_edges1 = vec![
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(1, 4, Outgoing, Incoming),
            MockEdge::new(3, 4, Outgoing, Incoming),
            MockEdge::new(3, 2, Outgoing, Incoming),
        ];

        mock_edges1[1].diff_snpmers = 1;
        mock_edges1[1].overlap_len_bases = 5000;
        mock_edges1[1].shared_minimizers = 500;

        mock_edges1[3].diff_snpmers = 3;
        mock_edges1[3].overlap_len_bases = 5000;
        mock_edges1[3].shared_minimizers = 1000;


        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.prune_lax_overlaps(8, None, 100., 0., false);
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 3);
    }

    #[test]
    fn test_prune_lax_overlaps_multiple(){
        // 1 ---------> 5
        // 1 -- > 2
        //   ↘ (x)
        // 3 --> 4 <-------6
        let mut mock_edges1 = vec![
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(1, 4, Outgoing, Incoming),
            MockEdge::new(3, 4, Outgoing, Incoming),
            MockEdge::new(1, 5, Outgoing, Incoming),
            MockEdge::new(6, 4, Outgoing, Incoming),
        ];

        mock_edges1[1].diff_snpmers = 1;
        mock_edges1[1].overlap_len_bases = 5000;
        mock_edges1[1].shared_minimizers = 500;

        mock_edges1[3].overlap_len_bases = 10000;
        mock_edges1[3].shared_minimizers = 1000;

        mock_edges1[4].overlap_len_bases = 10000;
        mock_edges1[4].shared_minimizers = 1000;


        let mut graph = mock_graph_from_edges(mock_edges1);
        graph.prune_lax_overlaps(8, None, 100., 0., false);
        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();

        assert_eq!(good_edges, 4);
    }

    #[test]
    fn test_prune_low_minimizer(){
        let mut mock_edges1 = vec![
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(1, 4, Outgoing, Incoming),
            MockEdge::new(3, 4, Outgoing, Incoming),
        ];

        mock_edges1[1].shared_minimizers = 5;
        let mut graph = mock_graph_from_edges(mock_edges1);

        graph.prune_low_minimizer_overlaps(10, 8., 5.);

        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(good_edges, 2);
    }

    #[test]
    fn test_prune_low_minimizer_relative(){
        let mut mock_edges1 = vec![
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(1, 4, Outgoing, Incoming),
            MockEdge::new(3, 4, Outgoing, Incoming),
        ];

        mock_edges1[0].shared_minimizers = 50;
        mock_edges1[1].shared_minimizers = 10;
        mock_edges1[2].shared_minimizers = 50;
        let mut graph = mock_graph_from_edges(mock_edges1);

        graph.prune_low_minimizer_overlaps(10, 100., 4.);

        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(good_edges, 2);
    }

    #[test]
    fn test_prune_low_minimizer_negative(){
        let mock_edges1 = vec![
            MockEdge::new(1, 2, Outgoing, Incoming),
            MockEdge::new(1, 4, Outgoing, Incoming),
            MockEdge::new(3, 4, Outgoing, Incoming),
        ];

        let mut graph = mock_graph_from_edges(mock_edges1);

        graph.prune_low_minimizer_overlaps(10, 8., 4.);

        let good_edges = graph.edges.iter().filter(|x| x.is_some()).count();
        assert_eq!(good_edges, 3);
    }
}

