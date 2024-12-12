use crate::types::*;
use fxhash::FxHashMap;
use fxhash::FxHashSet;

// Common node trait - both ReadData and UnitigNode share these properties
pub trait GraphNode {
    fn in_edges(&self) -> &[EdgeIndex];
    fn out_edges(&self) -> &[EdgeIndex];
    fn both_edges(&self) -> impl Iterator<Item = &EdgeIndex> {
        self.in_edges().iter().chain(self.out_edges().iter())
    }

    fn edges_direction(&self, direction: &Direction) -> &[EdgeIndex] {
        match direction {
            Direction::Incoming => self.in_edges(),
            Direction::Outgoing => self.out_edges(),
        }
    }

    fn edges_direction_reverse(&self, direction: &Direction) -> &[EdgeIndex] {
        match direction {
            Direction::Incoming => self.out_edges(),
            Direction::Outgoing => self.in_edges(),
        }
    }

    fn in_edges_mut(&mut self) -> &mut Vec<EdgeIndex>;
    fn out_edges_mut(&mut self) -> &mut Vec<EdgeIndex>;
}

// Common edge trait - both ReadOverlapEdgeTwin and UnitigEdge share these
pub trait GraphEdge {
    fn node1(&self) -> NodeIndex;
    fn node2(&self) -> NodeIndex;
    fn orientation1(&self) -> bool;
    fn orientation2(&self) -> bool;
    fn other_node(&self, node: NodeIndex) -> NodeIndex {
        if node == self.node1() {
            return self.node2();
        } else if node == self.node2() {
            return self.node1();
        } else {
            panic!("Node not in edge");
        }
    }
    fn get_orientation(&self, node1: NodeIndex, node2: NodeIndex) -> (bool, bool){
        if self.node1() == node1 && self.node2() == node2{
            (self.orientation1(), self.orientation2())
        }
        else if self.node1() == node2 && self.node2() == node1{
            (!self.orientation2(), !self.orientation1())
        }
        else{
            panic!("Nodes do not match edge")
        }
    }
    //o>--->o is outgoing for n1, incoming for n2
    //o<---->o is incoming for n1, incoming for n2. etc
    fn node_edge_direction(&self, node: &NodeIndex) -> Direction {
        if self.node1() == *node {
            if self.orientation1() {
                Direction::Outgoing
            } else {
                Direction::Incoming
            }
        } else if self.node2() == *node {
            if self.orientation2() {
                Direction::Incoming
            } else {
                Direction::Outgoing
            }
        } else {
            panic!("Edge does not connect to node");
        }
    }
}

// Base implementation of a bidirected graph that both can use
#[derive(Debug, Clone)]
pub struct BidirectedGraph<N, E> {
    pub nodes: FxHashMap<NodeIndex, N>,
    pub edges: Vec<Option<E>>,
}

impl<N: GraphNode + std::fmt::Debug, E: GraphEdge + std::fmt::Debug> BidirectedGraph<N, E> {

    pub fn remove_edges(&mut self, edges: FxHashSet<EdgeIndex>) {
        for node in self.nodes.values_mut() {
            node.in_edges_mut().retain(|x| !edges.contains(x));
            node.out_edges_mut().retain(|x| !edges.contains(x));
        }
        for edge_id in edges {
            self.edges[edge_id] = None;
        }
    }

    //Checks the connextion if it's valid. Returns None if it's not, e.g. the next node is a branch
    //point that can not be integrated
    fn valid_unitig_connection(
        &self,
        curr: NodeIndex,
        next_node: NodeIndex,
        orientation: (bool, bool),
    ) -> bool {
        let deg1;
        if orientation.0 {
            deg1 = self.nodes.get(&curr).unwrap().out_edges().len();
        } else {
            deg1 = self.nodes.get(&curr).unwrap().in_edges().len();
        }
        let deg2;
        if orientation.1 {
            deg2 = self.nodes.get(&next_node).unwrap().in_edges().len();
        } else {
            deg2 = self.nodes.get(&next_node).unwrap().out_edges().len();
        }

        if deg1 == 1 && deg2 == 1 {
            return true;
        } else {
            return false;
        }
    }

    fn traverse_inout_nonbranch(&self, node_idx: NodeIndex, previous_node: Option<NodeIndex>) -> Vec<EdgeIndex> {
        let node = self.nodes.get(&node_idx).unwrap();
        let in_deg = node.in_edges().len();
        let out_deg = node.out_edges().len();
        if in_deg == 1 && out_deg == 1 {
            if previous_node.is_none() {
                return vec![node.in_edges()[0], node.out_edges()[0]];
            } else {
                for &edge_idx in node.both_edges() {
                    let edge = self.edges.get(edge_idx).unwrap().as_ref().unwrap();
                    if edge.node1() != previous_node.unwrap() && edge.node2() != previous_node.unwrap() {
                        return vec![edge_idx];
                    }
                }
                //Circular o --> x
                //          ^  // 
                for &edge in node.both_edges() {
                    let edge = self.edges.get(edge).unwrap().as_ref().unwrap();
                    log::trace!("{:?}",&edge);
                }
                log::trace!("{}, {}, {:?}, {:?}, {:?}", in_deg, out_deg, node.in_edges(), node.out_edges(), previous_node);
                return vec![];
            }
        } else if in_deg == 1 && out_deg != 1 {
            return vec![node.in_edges()[0]];
        } else if in_deg != 1 && out_deg == 1 {
            return vec![node.out_edges()[0]];
        } else {
            return vec![];
        }
    }

    pub fn find_non_branching_paths(&self) -> Vec<(Vec<NodeIndex>, Vec<EdgeIndex>)> {
        let mut paths = Vec::new();
        let mut visited_nodes = FxHashSet::default();

        // Start from each potential unitig endpoint
        for (&start_node, _) in &self.nodes {
            // Skip if we've already processed this node
            if visited_nodes.contains(&start_node) {
                continue;
            }

            // Start a new path if this is an endpoint or branch point
            let mut path_both = [(vec![], vec![]), (vec![], vec![])];
            visited_nodes.insert(start_node);

            let edges_to_search = self.traverse_inout_nonbranch(start_node, None);

            for (i, edge_idx) in edges_to_search.into_iter().enumerate() {
                let mut current_edges = Vec::new();
                let mut curr_edge = edge_idx;
                let mut current_path = vec![];
                let mut prev_node = Some(start_node);
                loop {
                    let edge = self.edges.get(curr_edge).unwrap().as_ref().unwrap();
                    let next_node = edge.other_node(prev_node.unwrap());
                    if visited_nodes.contains(&next_node) {
                        path_both[i] = (current_path, current_edges);
                        break;
                    }
                    let orientation = edge.get_orientation(prev_node.unwrap(), next_node);
                    if self.valid_unitig_connection(prev_node.unwrap(), next_node, orientation) {
                        current_edges.push(curr_edge);
                        current_path.push(next_node);
                        visited_nodes.insert(next_node);
                        let current_node = next_node;
                        let next = self.traverse_inout_nonbranch(current_node, prev_node);
                        if next.len() != 1 {
                            log::trace!("{} Unitig path is not linear", start_node);
                            path_both[i] = (current_path, current_edges);
                            break;
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

        // Implementation that works for both graph types
    }

    pub fn remove_nodes(&mut self, nodes_to_remove: &[NodeIndex], keep_as_output: bool) {
        let mut touched_nodes = FxHashSet::default();
        let remove_set = FxHashSet::from_iter(nodes_to_remove.iter());
        for &node_idx in nodes_to_remove {
            if let Some(node) = self.nodes.get(&node_idx) {
                let edge_lists = self.nodes.get(&node_idx).unwrap().both_edges();
                for &edge_idx in edge_lists {
                    if let Some(Some(edge)) = self.edges.get(edge_idx) {
                        let other_idx = if edge.node1() == node_idx {
                            edge.node2()
                        } else {
                            edge.node1()
                        };
                        touched_nodes.insert(other_idx);
                    }
                }
                // Remove edges
                for &edge_idx in node.both_edges(){
                    self.edges[edge_idx] = None;
                }
            }
        }
        for node_idx in touched_nodes{
            let node = self.nodes.get_mut(&node_idx).unwrap();
            node.in_edges_mut().retain(|&edge| {
                if let Some(Some(edge)) = self.edges.get(edge) {
                    !remove_set.contains(&edge.node1())
                } else {
                    false
                }
            });
            node.out_edges_mut().retain(|&edge| {
                if let Some(Some(edge)) = self.edges.get(edge) {
                    !remove_set.contains(&edge.node2())
                } else {
                    false
                }
            });
        }

        // Remove nodes
        for node_idx in nodes_to_remove {
            self.nodes.remove(node_idx);
        }
    }

}

pub fn orientation_list<E: GraphEdge>(node_indices: &[NodeIndex], edges: &[E]) -> Vec<bool> {
    let mut orientations = Vec::new();
    if node_indices.len() == 0 {
        panic!("No reads in a graph node; something went wrong.");
    }
    else if node_indices.len() == 1{
        return vec![true];
    }
    else {
        for i in 0..node_indices.len() - 1{
            let node = &node_indices[i];
            let next_node = &node_indices[i + 1];
            let edge = &edges[i];
            let orientation = edge.get_orientation(*node, *next_node);
            if i == 0{
                orientations.push(orientation.0);
            }
            orientations.push(orientation.1);
        }
    }
    return orientations
}

