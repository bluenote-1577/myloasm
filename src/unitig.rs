use crate::twin_graph::*;
use crate::types::*;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use std::io::BufWriter;
use std::io::Write;
use bio_seq::prelude::*;

#[derive(Debug, Clone, PartialEq)]
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
                range = (carryover, reads[ind].base_length - overlap_lengths[i]);
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

#[derive(Debug, Clone, PartialEq)]
pub struct UnitigGraph {
    pub unitigs: FxHashMap<NodeIndex, UnitigNode>,
    pub edges: Vec<Option<UnitigEdge>>,
}

impl UnitigGraph {
    // Constructor that creates unitig graph from overlap graph
    pub fn from_overlap_graph(overlap_graph: &OverlapTwinGraph) -> Self {
        // Start with empty graph
        let mut unitig_graph = UnitigGraph {
            unitigs: FxHashMap::default(),
            edges: Vec::new(),
        };

        let nbp = find_non_branching_paths(overlap_graph);
        for (nodevec, edgevec) in nbp {
            let overlaps: Vec<ReadOverlapEdgeTwin> = edgevec
                .iter()
                .map(|&edge_idx| overlap_graph.edges[edge_idx].as_ref().unwrap().clone())
                .collect();
            let read_length: usize = nodevec
                .iter()
                .map(|&node_idx| overlap_graph.reads[node_idx].base_length)
                .sum();
            let overlap_lengths: usize = overlaps
                .iter()
                .map(|overlap| overlap.overlap_len_bases)
                .sum();
            let unitig_length = read_length - overlap_lengths;

            let mut oris = vec![];
            for i in 0..nodevec.len() - 1 {
                let read1 = nodevec[i];
                let read2 = nodevec[i + 1];
                let edge = &overlaps[i];
                let orientation = edge.get_orientation(read1, read2);
                if i == 0 {
                    oris.push(orientation.0);
                }
                oris.push(orientation.1);
            }
            if nodevec.len() == 1 {
                oris.push(true);
            }

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
                .unitigs
                .insert(unitig_graph.unitigs.len(), unitig);
        }

        let mut terminal_edge_to_unitigs = FxHashMap::default();
        let mut unitigs_to_terminal_edges = FxHashMap::default();
        let mut edge_assignments = vec![];
        for (id, unitig) in unitig_graph.unitigs.iter() {
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
            let uni1 = unitig_graph.unitigs.get_mut(&from_id).unwrap();

            if ori1 {
                uni1.out_edges.push(edge_id);
            } else {
                uni1.in_edges.push(edge_id);
            }

            let uni2 = unitig_graph.unitigs.get_mut(&to_id).unwrap();
            if ori2 {
                uni2.in_edges.push(edge_id);
            } else {
                uni2.out_edges.push(edge_id);
            }
        }
        unitig_graph
    }

    // Convert to GFA format
    pub fn to_gfa(&self, filename: &str, output_readgroups: bool, graph: &OverlapTwinGraph) {
        let mut bufwriter = BufWriter::new(std::fs::File::create(filename).unwrap());
        let mut gfa = String::new();

        // Header
        gfa.push_str("H\tVN:Z:1.0\n");

        // Segments
        for (id, unitig) in &self.unitigs {
            let (base_seq, ranges) = unitig.raw_consensus(&graph.reads);
            gfa.push_str(&format!(
                "S\t{}\t{}\tLN:i:{}\tDP:f:{:.1}\n",
                id,
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
                    let read = &graph.reads[*read_idx];
                    let range = ranges[count];
                    let length = range.1 - range.0;
                    gfa.push_str(&format!(
                        "a\t{},{}-{}\t{}\t{}\t{}\t{}\t{}\n",
                        id, range.0, range.1, read_idx, curr_pos, read.id, ori_string, length
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

            gfa.push_str(&format!(
                "L\t{}\t{}\t{}\t{}\t{}M\n",
                edge.from_unitig,
                edge.to_unitig,
                from_orient,
                to_orient,
                edge.overlap.overlap_len_bases
            ));
        }

        write!(bufwriter, "{}", gfa).unwrap();
    }

    pub fn test_consistent_left_right_edges(&self) {
        for (_, unitig) in self.unitigs.iter() {
            let mut in_readset = FxHashSet::default();
            for edge in unitig.in_edges.iter() {
                let e = self.edges[*edge].as_ref().unwrap();
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
        for (id, _) in self.unitigs.iter() {
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
                let unitig = self.unitigs.get(&node).unwrap();
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
            let unitig = self.unitigs.get_mut(&node).unwrap();
            unitig.left_cut = left_cut;
            unitig.right_cut = right_cut;
        }
    }

    // Graph manipulation methods
    pub fn remove_tips(&mut self, min_length: usize) {
        // Implementation
    }

    pub fn pop_bubbles(&mut self, max_length: usize, max_divergence: f64) {
        // Implementation
    }

    fn cut_node_and_inc_edge(&self, node: &NodeIndex, inc: Option<&UnitigEdge>) -> (usize, usize) {
        let cut_dir_edges;
        let mut direction = Direction::Incoming;
        let unitig = self.unitigs.get(node).unwrap();
        if unitig.in_edges.len() + unitig.out_edges.len() == 0{
            return (0, 0);
        }
        if let Some(edge) = inc {
            //Incoming direction
            direction = unitig_node_edge_direction(node, &edge);
            if direction == Direction::Incoming {
                cut_dir_edges = &unitig.in_edges;
            } else {
                cut_dir_edges = &unitig.out_edges;
            }
        } else {
            cut_dir_edges = &unitig.in_edges;
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
            dbg!(&self.unitigs[node]);
            dbg!(&cut_dir_edges);
            dbg!(direction);
            dbg!(inc);
            dbg!(node);
            panic!();
        }
    }
}

fn unitig_node_edge_direction(node: &NodeIndex, edge: &UnitigEdge) -> Direction {
    if edge.from_unitig == *node {
        if edge.f1 {
            Direction::Outgoing
        } else {
            Direction::Incoming
        }
    } else if edge.to_unitig == *node {
        if edge.f2 {
            Direction::Incoming
        } else {
            Direction::Outgoing
        }
    } else {
        panic!("Edge does not connect to node");
    }
}

// Helper to check if a node is a branch point considering balanced degrees
fn is_branch_point(graph: &OverlapTwinGraph, node_idx: NodeIndex) -> bool {
    if let Some(node) = graph.nodes.get(&node_idx) {
        // A non-branching internal node must have in_deg = out_deg = 1
        let in_deg = node.in_edges.len();
        let out_deg = node.out_edges.len();

        in_deg != 1 || out_deg != 1
    } else {
        true
    }
}

//Traverse directions based on previous node
//return left and right traversals (in/out)
fn traverse_graph_for_unitigs(
    graph: &OverlapTwinGraph,
    node_idx: NodeIndex,
    previous_node: Option<NodeIndex>,
) -> Vec<EdgeIndex> {
    let node = graph.nodes.get(&node_idx).unwrap();
    let in_deg = node.in_edges.len();
    let out_deg = node.out_edges.len();
    if in_deg == 1 && out_deg == 1 {
        if previous_node.is_none() {
            return vec![node.in_edges[0], node.out_edges[0]];
        } else {
            for edge_idx in node.in_edges.iter().chain(node.out_edges.iter()) {
                let edge = graph.edges.get(*edge_idx).unwrap().as_ref().unwrap();
                if edge.node1 != previous_node.unwrap() && edge.node2 != previous_node.unwrap() {
                    return vec![*edge_idx];
                }
            }
            return vec![];
        }
    } else if in_deg == 1 && out_deg != 1 {
        return vec![node.in_edges[0]];
    } else if in_deg != 1 && out_deg == 1 {
        return vec![node.out_edges[0]];
    } else {
        return vec![];
    }
}

// Helper to get all valid connections for a node
fn get_node_connections(
    graph: &OverlapTwinGraph,
    node_idx: NodeIndex,
) -> Vec<(NodeIndex, EdgeIndex, bool)> {
    let mut connections = Vec::new();

    if let Some(node) = graph.nodes.get(&node_idx) {
        // Process edges maintaining in/out distinction
        for &edge_idx in &node.out_edges {
            if let Some(Some(edge)) = graph.edges.get(edge_idx) {
                let (next_node, forward) = if edge.node1 == node_idx {
                    (edge.node2, edge.forward1)
                } else {
                    (edge.node1, edge.forward2)
                };
                connections.push((next_node, edge_idx, forward));
            }
        }

        for &edge_idx in &node.in_edges {
            if let Some(Some(edge)) = graph.edges.get(edge_idx) {
                let (next_node, forward) = if edge.node1 == node_idx {
                    (edge.node2, !edge.forward1) // Flip orientation for in-edges
                } else {
                    (edge.node1, !edge.forward2) // Flip orientation for in-edges
                };
                connections.push((next_node, edge_idx, forward));
            }
        }
    }

    connections
}

// Helper to find the next node in a path, considering bidirectionality
fn get_next_unvisited(
    graph: &OverlapTwinGraph,
    node: NodeIndex,
    visited_edges: &std::collections::HashSet<EdgeIndex>,
) -> Option<(NodeIndex, EdgeIndex, bool)> {
    let connections = get_node_connections(graph, node);

    // Find first unvisited connection
    for (next_node, edge_idx, forward) in connections {
        if !visited_edges.contains(&edge_idx) {
            return Some((next_node, edge_idx, forward));
        }
    }
    None
}

fn find_non_branching_paths(graph: &OverlapTwinGraph) -> Vec<(Vec<NodeIndex>, Vec<EdgeIndex>)> {
    let mut paths = Vec::new();
    let mut visited_nodes = FxHashSet::default();

    // Start from each potential unitig endpoint
    for (&start_node, _) in &graph.nodes {
        // Skip if we've already processed this node
        if visited_nodes.contains(&start_node) {
            continue;
        }

        // Start a new path if this is an endpoint or branch point
        let mut path_both = [(vec![], vec![]), (vec![], vec![])];
        visited_nodes.insert(start_node);

        let edges_to_search = traverse_graph_for_unitigs(graph, start_node, None);

        for (i, edge_idx) in edges_to_search.into_iter().enumerate() {
            let mut current_edges = Vec::new();
            let mut curr_edge = edge_idx;
            let mut current_path = vec![];
            let mut prev_node = Some(start_node);
            loop {
                let edge = graph.edges.get(curr_edge).unwrap().as_ref().unwrap();
                let next_node = graph.other_node(prev_node.unwrap(), edge);
                if visited_nodes.contains(&next_node) {
                    path_both[i] = (current_path, current_edges);
                    break;
                }
                let orientation = edge.get_orientation(prev_node.unwrap(), next_node);
                if valid_unitig_connection(graph, prev_node.unwrap(), next_node, orientation) {
                    current_edges.push(curr_edge);
                    current_path.push(next_node);
                    visited_nodes.insert(next_node);
                    let current_node = next_node;
                    let next = traverse_graph_for_unitigs(graph, current_node, prev_node);
                    if next.len() != 1 {
                        panic!("Unitig path is not linear");
                    }
                    curr_edge = next[0];
                    prev_node = Some(current_node);
                } else {
                    path_both[i] = (current_path, current_edges);
                    break;
                }
            }
        }

        let mut unitig_node_path = vec![];
        let mut unitig_edge_path = vec![];
        for node in path_both[0].0.iter().rev() {
            unitig_node_path.push(*node);
        }
        for edge in path_both[0].1.iter().rev() {
            unitig_edge_path.push(*edge);
        }
        unitig_node_path.push(start_node);
        for node in path_both[1].0.iter() {
            unitig_node_path.push(*node);
        }
        for edge in path_both[1].1.iter() {
            unitig_edge_path.push(*edge);
        }
        paths.push((unitig_node_path, unitig_edge_path));
    }

    paths
}

//Checks the connextion if it's valid. Returns None if it's not, e.g. the next node is a branch
//point that can not be integrated
fn valid_unitig_connection(
    graph: &OverlapTwinGraph,
    curr: NodeIndex,
    next_node: NodeIndex,
    orientation: (bool, bool),
) -> bool {
    let deg1;
    if orientation.0 {
        deg1 = graph.nodes.get(&curr).unwrap().out_edges.len();
    } else {
        deg1 = graph.nodes.get(&curr).unwrap().in_edges.len();
    }
    let deg2;
    if orientation.1 {
        deg2 = graph.nodes.get(&next_node).unwrap().in_edges.len();
    } else {
        deg2 = graph.nodes.get(&next_node).unwrap().out_edges.len();
    }

    if deg1 == 1 && deg2 == 1 {
        return true;
    } else {
        return false;
    }
}
