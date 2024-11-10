use crate::twin_graph::*;
use crate::graph::*;
use crate::types::*;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use std::io::BufWriter;
use std::io::Write;
use bio_seq::prelude::*;

#[derive(Debug, Clone, PartialEq, Default)]
pub struct UnitigNode {
    // Vector of read indices that make up this unitig in order
    pub read_indices_ori: Vec<(NodeIndex, bool)>,
    // Internal overlaps between consecutive reads in the unitig
    pub internal_overlaps: Vec<ReadOverlapEdgeTwin>,
    // Length of the unitig in bases
    pub length: usize,
    // Coverage depth
    pub depth: Option<f64>,
    // Incoming and outgoing edges
    pub in_edges: Vec<EdgeIndex>,
    pub out_edges: Vec<EdgeIndex>,

    pub right_cut: usize,
    pub left_cut: usize,
}

impl GraphNode for UnitigNode {
    fn in_edges(&self) -> &[EdgeIndex] { &self.in_edges }
    fn out_edges(&self) -> &[EdgeIndex] { &self.out_edges }
    fn in_edges_mut(&mut self) -> &mut Vec<EdgeIndex> {
        &mut self.in_edges
    }
    fn out_edges_mut(&mut self) -> &mut Vec<EdgeIndex> {
        &mut self.out_edges
    }
}

impl UnitigNode {
    // Calculate the consensus sequence of the unitig
    pub fn raw_consensus(&self, reads: &[TwinRead]) -> (Seq<Dna>, Vec<(usize,usize)>){
        let mut base_seq = Seq::new();
        let mut ranges = vec![];
        let mut carryover = 0;
        let overlap_lengths = self.internal_overlaps.iter().map(|overlap| overlap.overlap_len_bases).collect::<Vec<_>>();
        for (i, indori) in self.read_indices_ori.iter().enumerate(){
            let ind = indori.0;
            let ori = indori.1;
            let mut range;
            if self.read_indices_ori.len() == 1{
                range = (self.left_cut, reads[ind].base_length - self.right_cut);
            }
            else if i == 0{
                range = (self.left_cut, reads[ind].base_length - overlap_lengths[0]);
            }
            else if i == self.read_indices_ori.len() - 1{
                range = (carryover, reads[ind].base_length - self.right_cut);
            }
            else{
                range = (carryover, reads[ind].base_length
                     - overlap_lengths[i]);
            }
            if range.0 >= range.1{
                range = (0,0);
                carryover += range.0 - range.1;
            }
            else{
                carryover = 0;
            }
            ranges.push(range);
            if ori{
                base_seq.append(&reads[ind].dna_seq[range.0..range.1]);
            }
            else{
                base_seq.append(&reads[ind].dna_seq.revcomp()[range.0..range.1]);
            }
        }
        return (base_seq, ranges);
    }
}

// UnitigEdge is essentially a wrapper around ReadOverlapEdgeTwin with unitig context
#[derive(Debug, Clone, PartialEq)]
pub struct UnitigEdge {
    // The actual overlap information
    pub overlap: ReadOverlapEdgeTwin,
    // Which reads from the unitigs are involved in this overlap
    pub from_read_idx: usize, // Index into the from_unitig's read_indices
    pub to_read_idx: usize,   // Index into the to_unitig's read_indices
    // The unitigs being connected
    pub from_unitig: NodeIndex,
    pub to_unitig: NodeIndex,

    pub f1: bool,
    pub f2: bool,
}

impl GraphEdge for UnitigEdge {
    fn node1(&self) -> NodeIndex { self.from_unitig }
    fn node2(&self) -> NodeIndex { self.to_unitig }
    fn orientation1(&self) -> bool { self.f1 }
    fn orientation2(&self) -> bool { self.f2 }
}

impl <'a> GraphEdge for &'a UnitigEdge {
    fn node1(&self) -> NodeIndex { self.from_unitig }
    fn node2(&self) -> NodeIndex { self.to_unitig }
    fn orientation1(&self) -> bool { self.f1 }
    fn orientation2(&self) -> bool { self.f2 }
}

pub type UnitigGraph = BidirectedGraph<UnitigNode,UnitigEdge>;

impl UnitigGraph {

    pub fn re_unitig(&mut self) {
        let mut new_unitig_graph = UnitigGraph {
            nodes: FxHashMap::default(),
            edges: Vec::new(),
        };

        let mut old_nodes_to_new_nodes = FxHashMap::default();
        let mut old_terminal_edges: FxHashSet<NodeIndex> = FxHashSet::default();
        let mut old_nodes_new_ori: FxHashMap<NodeIndex, bool> = FxHashMap::default();
        let unitig_nbp = self.find_non_branching_paths();
        for (nodepath, edgepath) in unitig_nbp{
            let mut new_read_indices_ori = vec![];
            let mut new_internal_overlaps = vec![];
            let overlaps = edgepath.iter().map(|&edge_idx| self.edges[edge_idx].as_ref().unwrap()).collect::<Vec<_>>();
            let mut new_unitig = UnitigNode::default();
            let unitig_path_ori = orientation_list(&nodepath,&overlaps);
            for (i,(&ori, &node_ind)) in unitig_path_ori.iter().zip(nodepath.iter()).enumerate(){
                old_nodes_to_new_nodes.insert(node_ind, new_unitig_graph.nodes.len());
                old_nodes_new_ori.insert(node_ind, ori);
                let unitig = self.nodes.get_mut(&node_ind).unwrap();
                let mut rdio = std::mem::take(&mut unitig.read_indices_ori);
                let internal_overlaps = std::mem::take(&mut unitig.internal_overlaps);
                if ori{
                    new_read_indices_ori.extend(rdio);
                    new_internal_overlaps.extend(internal_overlaps);
                }
                else{
                    rdio = rdio.into_iter().rev().map(|(ind,ori)| (ind,!ori)).collect();
                    new_read_indices_ori.extend(rdio);
                    new_internal_overlaps.extend(internal_overlaps.into_iter().rev());
                }
                if i != nodepath.len() - 1{
                    new_internal_overlaps.push(self.edges[edgepath[i]].as_ref().unwrap().overlap.clone());
                }
                new_unitig.length += unitig.length;
            }
            new_unitig.read_indices_ori = new_read_indices_ori;
            new_unitig.internal_overlaps = new_internal_overlaps;
            new_unitig_graph.nodes.insert(new_unitig_graph.nodes.len(), new_unitig);

            let last_node = &self.nodes[&nodepath[nodepath.len()-1]];
            if last_node.in_edges().len() > 1{
                old_terminal_edges.extend(last_node.in_edges());
            }
            if last_node.out_edges().len() > 1{
                old_terminal_edges.extend(last_node.out_edges());
            }
            let first_node = &self.nodes[&nodepath[0]];
            if first_node.in_edges().len() > 1{
                old_terminal_edges.extend(first_node.in_edges());
            }
            if first_node.out_edges().len() > 1{
                old_terminal_edges.extend(first_node.out_edges());
            }
        }

        //Re-edge the graph using the old unitig edges but mapped to the new nodes
        for old_edge_id in old_terminal_edges{
            let old_edge = self.edges[old_edge_id].as_ref().unwrap();
            let new_node1_ind = old_nodes_to_new_nodes[&old_edge.node1()];
            let new_node2_ind = old_nodes_to_new_nodes[&old_edge.node2()];
            let old_orientations = self.edges[old_edge_id].as_ref().unwrap().get_orientation(old_edge.node1(), old_edge.node2());
            let new_orientations = (old_nodes_new_ori[&old_edge.node1()], old_nodes_new_ori[&old_edge.node2()]);
            
            let new_unitig_edge = UnitigEdge{
                overlap: self.edges[old_edge_id].as_ref().unwrap().overlap.clone(),
                from_read_idx: self.edges[old_edge_id].as_ref().unwrap().from_read_idx,
                to_read_idx: self.edges[old_edge_id].as_ref().unwrap().to_read_idx,
                from_unitig: new_node1_ind,
                to_unitig: new_node2_ind,
                f1: new_orientations.0 == old_orientations.0,
                f2: new_orientations.1 == old_orientations.1,
            };

            let new_n1 = new_unitig_graph.nodes.get_mut(&new_node1_ind).unwrap();

            if new_unitig_edge.f1{
                new_n1.out_edges.push(new_unitig_graph.edges.len());
            }
            else{
                new_n1.in_edges.push(new_unitig_graph.edges.len());
            }

            let new_n2 = new_unitig_graph.nodes.get_mut(&new_node2_ind).unwrap();
            if new_unitig_edge.f2{
                new_n2.in_edges.push(new_unitig_graph.edges.len());
            }
            else{
                new_n2.out_edges.push(new_unitig_graph.edges.len());
            }

            new_unitig_graph.edges.push(Some(new_unitig_edge));
        }

        *self = new_unitig_graph;

    }
    // Constructor that creates unitig graph from overlap graph
    pub fn from_overlap_graph(overlap_graph: &OverlapTwinGraph, reads: &[TwinRead]) -> Self {
        // Start with empty graph
        let mut unitig_graph = UnitigGraph {
            nodes: FxHashMap::default(),
            edges: Vec::new(),
        };

        let nbp = overlap_graph.find_non_branching_paths();
        for (nodevec, edgevec) in nbp {
            let overlaps: Vec<ReadOverlapEdgeTwin> = edgevec
                .iter()
                .map(|&edge_idx| overlap_graph.edges[edge_idx].as_ref().unwrap().clone())
                .collect();
            let read_length: usize = nodevec
                .iter()
                .map(|&node_idx| reads[node_idx].base_length)
                .sum();
            //DOESN"T WORK
            let overlap_lengths: usize = overlaps
                .iter()
                .map(|overlap| overlap.overlap_len_bases)
                .sum();
            let unitig_length = read_length - overlap_lengths;

            let oris = orientation_list(&nodevec, &overlaps);
            
            let read_indices_ori = (0..nodevec.len()).map(|i| (nodevec[i], oris[i])).collect();

            let unitig = UnitigNode {
                read_indices_ori,
                internal_overlaps: overlaps,
                length: unitig_length,
                depth: None,
                in_edges: Vec::new(),
                out_edges: Vec::new(),
                right_cut: 0,
                left_cut: 0,
            };
            unitig_graph
                .nodes
                .insert(unitig_graph.nodes.len(), unitig);
        }

        let mut terminal_edge_to_unitigs = FxHashMap::default();
        let mut unitigs_to_terminal_edges = FxHashMap::default();
        let mut edge_assignments = vec![];
        let mut seen_edges = FxHashSet::default();
        for (id, unitig) in unitig_graph.nodes.iter() {
            let terminal_reads;
            let mut terminal_edges = vec![];
            if unitig.read_indices_ori.len() == 1 {
                terminal_reads = vec![unitig.read_indices_ori[0]];
            } else {
                terminal_reads = vec![
                    unitig.read_indices_ori[0],
                    unitig.read_indices_ori[unitig.read_indices_ori.len() - 1],
                ];
            }
            //<- [x-x-x-x-x]-> Get the terminal edges
            for (read_idx, read_ori) in terminal_reads {
                let read_in_graph = overlap_graph.nodes.get(&read_idx).unwrap();
                for edge_idx in read_in_graph
                    .in_edges
                    .iter()
                    .chain(read_in_graph.out_edges.iter())
                {
                    let edge = overlap_graph.edges[*edge_idx].as_ref().unwrap();
                    if unitig.internal_overlaps.len() == 0 {
                        terminal_edges.push((edge_idx, read_idx, read_ori));
                        terminal_edge_to_unitigs
                            .entry(*edge_idx)
                            .or_insert_with(Vec::new)
                            .push((id, read_idx, read_ori));
                    } else if edge != &unitig.internal_overlaps[0]
                        && edge != &unitig.internal_overlaps[unitig.internal_overlaps.len() - 1]
                    {
                        terminal_edges.push((edge_idx, read_idx, read_ori));
                        terminal_edge_to_unitigs
                            .entry(*edge_idx)
                            .or_insert_with(Vec::new)
                            .push((id, read_idx, read_ori));
                    }
                }
            }
            unitigs_to_terminal_edges.insert(id, terminal_edges);
        }

        //Check for overlaps between terminal reads of unitigs
        for (unitig_id, term) in unitigs_to_terminal_edges {
            for (edge_id, read_id, read_ori) in term {
                let other_unitigs = terminal_edge_to_unitigs.get(&edge_id).unwrap();
                for (unitig_id2, read_id2, read_ori2) in other_unitigs {
                    if unitig_id == *unitig_id2 {
                        continue;
                    }

                    let sorted_uni_ids = if unitig_id < *unitig_id2 {
                        (*unitig_id, **unitig_id2)
                    } else {
                        (**unitig_id2, *unitig_id)
                    };

                    if seen_edges.contains(&sorted_uni_ids) {
                        continue;
                    }
                    else{
                        seen_edges.insert(sorted_uni_ids);
                    }

                    let overlap_edge = overlap_graph.edges[*edge_id].as_ref().unwrap();
                    let orientation = overlap_edge.get_orientation(read_id, *read_id2);

                    //If the terminal orientation is -, the read will be reverse complemented in
                    //the unitig. So the link, which assumes the original (non-RC) read, must be
                    //flipped.
                    let ori1 = read_ori == orientation.0;
                    let ori2 = *read_ori2 == orientation.1;

                    let unitig_edge = UnitigEdge {
                        overlap: overlap_edge.clone(),
                        from_read_idx: read_id,
                        to_read_idx: *read_id2,
                        from_unitig: *unitig_id,
                        to_unitig: **unitig_id2,
                        f1: ori1,
                        f2: ori2,
                    };
                    unitig_graph.edges.push(Some(unitig_edge));
                    edge_assignments.push((
                        *unitig_id,
                        **unitig_id2,
                        unitig_graph.edges.len() - 1,
                        ori1,
                        ori2,
                    ));
                }
            }
        }

        for (from_id, to_id, edge_id, ori1, ori2) in edge_assignments {
            let uni1 = unitig_graph.nodes.get_mut(&from_id).unwrap();

            if ori1 {
                uni1.out_edges.push(edge_id);
            } else {
                uni1.in_edges.push(edge_id);
            }

            let uni2 = unitig_graph.nodes.get_mut(&to_id).unwrap();
            if ori2 {
                uni2.in_edges.push(edge_id);
            } else {
                uni2.out_edges.push(edge_id);
            }
        }
        unitig_graph
    }

    // Convert to GFA format
    pub fn to_gfa(&self, filename: &str, output_readgroups: bool, reads: &[TwinRead]) {
        let mut bufwriter = BufWriter::new(std::fs::File::create(filename).unwrap());
        let mut gfa = String::new();

        // Header
        gfa.push_str("H\tVN:Z:1.0\n");

        // Segments
        for (_id, unitig) in &self.nodes {
            let (base_seq, ranges) = unitig.raw_consensus(reads);
            gfa.push_str(&format!(
                "S\tu{}\t{}\tLN:i:{}\tDP:f:{:.1}\n",
                unitig.read_indices_ori[0].0,
                //String::from_utf8(unitig.raw_consensus()).unwrap(),
                base_seq,
                unitig.length,
                unitig.depth.unwrap_or(0.)
            ));

            if output_readgroups {
                let mut curr_pos = 0;
                let mut count = 0;
                for (read_idx, read_ori) in unitig.read_indices_ori.iter() {
                    let ori_string = if *read_ori { "+" } else { "-" };
                    let read = &reads[*read_idx];
                    let range = ranges[count];
                    let length = range.1 - range.0;
                    gfa.push_str(&format!(
                        "a\tu{},{}-{}\t{}\t{}\t{}\t{}\t{}\n",
                        unitig.read_indices_ori[0].0,
                        range.0, range.1, read_idx, curr_pos, read.id, ori_string, length
                    ));
                    curr_pos += length;
                    count += 1;
                }
            }
        }

        // Links
        for edge in self.edges.iter().flatten() {
            let from_orient = if edge.f1 { "+" } else { "-" };
            let to_orient = if edge.f2 { "+" } else { "-" };
            let id1 = self.nodes[&edge.from_unitig].read_indices_ori[0].0;
            let id2 = self.nodes[&edge.to_unitig].read_indices_ori[0].0;

            gfa.push_str(&format!(
                "L\tu{}\t{}\tu{}\t{}\t{}M\n",
                id1,
                from_orient,
                id2,
                to_orient,
                edge.overlap.overlap_len_bases
            ));
        }

        write!(bufwriter, "{}", gfa).unwrap();
    }

    pub fn test_consistent_left_right_edges(&self) {
        for (_, unitig) in self.nodes.iter() {
            let mut in_readset = FxHashSet::default();
            for edge in unitig.in_edges.iter() {
                let e = self.edges[*edge].as_ref().unwrap();
                // If the edge's from is the left end, make sure
                //the from_orientation is false. 
                if e.from_read_idx == unitig.read_indices_ori[0].0 {
                    in_readset.insert(e.from_read_idx);
                    debug_assert!(
                        !e.f1,
                        "Edge: {:?} and node {:?}",
                        e, unitig.read_indices_ori
                    );
                } else if e.to_read_idx == unitig.read_indices_ori[0].0 {
                    in_readset.insert(e.to_read_idx);
                    debug_assert!(e.f2, "Edge: {:?} and node {:?}", e, unitig.read_indices_ori);
                } else if e.from_read_idx
                    == unitig.read_indices_ori[unitig.read_indices_ori.len() - 1].0
                {
                    in_readset.insert(e.from_read_idx);
                    debug_assert!(!e.f1, "Edge: {:?} and node{:?}", e, unitig.read_indices_ori);
                } else if e.to_read_idx
                    == unitig.read_indices_ori[unitig.read_indices_ori.len() - 1].0
                {
                    in_readset.insert(e.to_read_idx);
                    debug_assert!(e.f2, "Edge: {:?} and node{:?}", e, unitig.read_indices_ori);
                } else {
                    panic!("Edge does not connect to unitig");
                }
            }
            if in_readset.len() > 1 && unitig.read_indices_ori.len() != 1 {
                dbg!(&unitig.read_indices_ori);
                dbg!(in_readset);
                let edges = unitig
                    .in_edges
                    .iter()
                    .map(|e| self.edges[*e].as_ref().unwrap())
                    .collect::<Vec<_>>();
                dbg!(edges);
                panic!("Incoming edges link to left and right end of unitig");
            }

            let mut out_readset = FxHashSet::default();
            for edge in unitig.out_edges.iter() {
                let e = self.edges[*edge].as_ref().unwrap();
                if e.from_read_idx == unitig.read_indices_ori[0].0 {
                    out_readset.insert(e.from_read_idx);
                    debug_assert!(e.f1, "Edge: {:?} and node {:?}", e, unitig.read_indices_ori);
                } else if e.to_read_idx == unitig.read_indices_ori[0].0 {
                    out_readset.insert(e.to_read_idx);
                    debug_assert!(
                        !e.f2,
                        "Edge: {:?} and node {:?}",
                        e, unitig.read_indices_ori
                    );
                } else if e.from_read_idx
                    == unitig.read_indices_ori[unitig.read_indices_ori.len() - 1].0
                {
                    out_readset.insert(e.from_read_idx);
                    debug_assert!(e.f1, "Edge: {:?} and node{:?}", e, unitig.read_indices_ori);
                } else if e.to_read_idx
                    == unitig.read_indices_ori[unitig.read_indices_ori.len() - 1].0
                {
                    out_readset.insert(e.to_read_idx);
                    debug_assert!(!e.f2, "Edge: {:?} and node{:?}", e, unitig.read_indices_ori);
                } else {
                    panic!("Edge does not connect to unitig");
                }
            }
            if out_readset.len() > 1 && unitig.read_indices_ori.len() != 1 {
                dbg!(&unitig.read_indices_ori);
                dbg!(out_readset);
                let edges = unitig
                    .out_edges
                    .iter()
                    .map(|e| self.edges[*e].as_ref().unwrap())
                    .collect::<Vec<_>>();
                dbg!(edges);
                panic!("Incoming edges link to left and right end of unitig: OUT");
            }
        }
    }

    pub fn cut_overlap_boundaries(&mut self) {
        let mut visited_nodes = FxHashSet::default();
        let mut changes = vec![];
        for (id, _) in self.nodes.iter() {
            if visited_nodes.contains(id) {
                continue;
            }
            let mut explore_vec = vec![(*id, None)];
            loop {
                if explore_vec.is_empty() {
                    break;
                }
                let (node, inc_edge) = explore_vec.pop().unwrap();
                let (left_cut, right_cut) = self.cut_node_and_inc_edge(&node, inc_edge);
                changes.push((node, left_cut, right_cut));
                visited_nodes.insert(node);
                let unitig = self.nodes.get(&node).unwrap();
                for edge in unitig.in_edges.iter().chain(unitig.out_edges.iter()) {
                    let edge = self.edges[*edge].as_ref().unwrap();
                    let other_node = if edge.from_unitig == node {
                        edge.to_unitig
                    } else {
                        edge.from_unitig
                    };
                    if !visited_nodes.contains(&other_node) {
                        explore_vec.push((other_node, Some(edge)));
                    }
                }
            }
        }
        for (node, left_cut, right_cut) in changes {
            let unitig = self.nodes.get_mut(&node).unwrap();
            unitig.left_cut = left_cut;
            unitig.right_cut = right_cut;
            if unitig.length > left_cut + right_cut {
                unitig.length -= left_cut + right_cut;
            } else {
//                dbg!(unitig.length);
//                dbg!(left_cut);
//                dbg!(right_cut);
//                dbg!(&unitig);
                unitig.length = 0;
            }
        }
    }

    pub fn remove_tips(&mut self, length: usize, num_reads: usize) {

        let node_to_sizeread_map = self.get_all_connected_components();
        let mut unitigs_to_remove = Vec::new();
        // First pass: find all dead ends
        let dead_ends: Vec<NodeIndex> = self.nodes.keys()
            .filter(|&node_idx| {
                let unitig = &self.nodes[node_idx];
                unitig.in_edges.is_empty() || unitig.out_edges.is_empty()
            })
            .copied()
            .collect();

        for dead_end_ind in dead_ends {
            if let Some(unitig) = self.nodes.get(&dead_end_ind) {
                let (bp_size_cc, reads_in_cc) = node_to_sizeread_map[&dead_end_ind];
                if unitig.length < length.min(bp_size_cc/10)
                    || unitig.read_indices_ori.len() < num_reads.min(reads_in_cc/10) {
                    unitigs_to_remove.push(dead_end_ind);
                }
            }
        }

        log::info!("Removing {} tips", unitigs_to_remove.len());
        self.remove_nodes(&unitigs_to_remove);
    }


    pub fn pop_bubbles(&mut self, _max_length: usize, _max_divergence: f64) {
        // Implementation
    }

    fn cut_node_and_inc_edge(&self, node: &NodeIndex, inc: Option<&UnitigEdge>) -> (usize, usize) {
        let cut_dir_edges;
        let mut direction = Direction::Incoming;
        let unitig = self.nodes.get(node).unwrap();
        if unitig.in_edges.len() + unitig.out_edges.len() == 0{
            return (0, 0);
        }
        if let Some(edge) = inc {
            //Incoming direction
            direction = edge.node_edge_direction(node);
            if direction == Direction::Incoming {
                cut_dir_edges = &unitig.in_edges;
            } else {
                cut_dir_edges = &unitig.out_edges;
            }
        } else {
            if unitig.in_edges.len() != 0{
                cut_dir_edges = &unitig.in_edges;
            }
            else{
                cut_dir_edges = &unitig.out_edges;
            }
        }
        let mut max_overlap = 0;
        for e_ind in cut_dir_edges.iter() {
            let edge = self.edges[*e_ind].as_ref().unwrap();
            if edge.overlap.overlap_len_bases > max_overlap {
                max_overlap = edge.overlap.overlap_len_bases;
            }
        }
        //Cut left or right depending on if cut_dir_edges corresponds to leftmost or rightmost read
        if unitig.read_indices_ori.len() == 1 {
            if direction == Direction::Incoming {
                return (max_overlap, 0);
            } else {
                return (0, max_overlap);
            }
        } else {
            for e_id in cut_dir_edges.iter() {
                let edge = self.edges[*e_id].as_ref().unwrap();
                let first_read = unitig.read_indices_ori[0].0;
                let last_read = unitig.read_indices_ori
                    [unitig.read_indices_ori.len() - 1]
                    .0;
                if edge.from_read_idx == first_read {
                    return (edge.overlap.overlap_len_bases, 0);
                } else if edge.to_read_idx == first_read {
                    return (edge.overlap.overlap_len_bases, 0);
                } else if edge.from_read_idx == last_read {
                    return (0, edge.overlap.overlap_len_bases);
                } else if edge.to_read_idx == last_read {
                    return (0, edge.overlap.overlap_len_bases);
                } else {
                    panic!("Edge does not connect to unitig");
                }
            }
            dbg!(&self.nodes[node]);
            dbg!(&cut_dir_edges);
            dbg!(direction);
            dbg!(inc);
            dbg!(node);
            panic!();
        }
    }


    fn get_all_connected_components(&self) -> FxHashMap<NodeIndex, (usize, usize)> {
        use std::collections::VecDeque;
        let mut visited = FxHashSet::default();
        let mut component_sizes = FxHashMap::default();

        // Process each node if not already visited
        for &start_node in self.nodes.keys() {
            if visited.contains(&start_node) {
                continue;
            }

            let mut current_visits = vec![];

            // BFS on this component
            let mut queue = VecDeque::new();
            queue.push_back(start_node);
            visited.insert(start_node);
            current_visits.push(start_node);

            let mut component_length = 0;
            let mut component_reads = 0;

            while let Some(curr_node) = queue.pop_front() {
                let unitig = &self.nodes[&curr_node];
                component_length += unitig.length;
                component_reads += unitig.read_indices_ori.len();

                // Process all neighbors
                for &edge_idx in unitig.in_edges.iter().chain(unitig.out_edges.iter()) {
                    let edge = self.edges[edge_idx].as_ref().unwrap();
                    let next_node = if edge.from_unitig == curr_node {
                        edge.to_unitig
                    } else {
                        edge.from_unitig
                    };

                    if !visited.contains(&next_node) {
                        visited.insert(next_node);
                        queue.push_back(next_node);
                        current_visits.push(next_node);
                    }
                }
            }

            // Set the component size for all nodes in this component
            for &node in current_visits.iter() {
                component_sizes.insert(node, (component_length, component_reads));
            }
        }

        component_sizes
    }

    

}





