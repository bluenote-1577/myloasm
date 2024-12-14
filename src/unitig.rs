use crate::cli::Cli;
use bio_seq::prelude::*;
use crate::graph::*;
use crate::mapping;
use crate::twin_graph::*;
use crate::types::*;
use fishers_exact::fishers_exact;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use rust_lapper::Lapper;
use std::io::BufWriter;
use std::io::Write;
use std::collections::VecDeque;
use std::panic;

#[derive(Debug, Clone, PartialEq, Default)]
pub struct UnitigNode {
    pub read_indices_ori: Vec<(NodeIndex, bool)>,
    pub internal_overlaps: Vec<ReadOverlapEdgeTwin>,
    pub read_names: Vec<String>,
    in_edges: Vec<EdgeIndex>,
    out_edges: Vec<EdgeIndex>,
    pub node_id: NodeIndex,
    pub node_hash_id: NodeIndex,
    pub mean_read_depth: Option<f64>, //Depth as measured by mean depth of reads within unitig
    base_info: BaseInfo,
    mapping_info: MappingInfo,
}

pub trait NodeSequence {
    fn base_seq(&self) -> &Seq<Dna>;
    fn read_positions_internal(&self) -> &Vec<(usize, usize)>;
    fn left_cut(&self) -> usize;
    fn right_cut(&self) -> usize;
    fn length(&self) -> usize;
    fn base_info_present(&self) -> bool;
    fn set_info(&mut self, info: BaseInfo);
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BaseInfo {
    pub base_seq: Seq<Dna>,
    pub read_positions_internal: Vec<(usize, usize)>,
    pub length: usize,
    pub left_cut: usize,
    pub right_cut: usize,
    pub present: bool,
}

impl PartialEq for MappingInfo {
    fn eq(&self, other: &Self) -> bool {
        self.median_depth == other.median_depth
            && self.mean_depth == other.mean_depth
            && self.mapping_boundaries.intervals == other.mapping_boundaries.intervals
    }
}

impl Default for MappingInfo {
    fn default() -> Self {
        MappingInfo {
            median_depth: 0.0,
            mean_depth: 0.0,
            mapping_boundaries: Lapper::new(vec![]),
            present: false,
            length: 0,
        }
    }
}

impl NodeMapping for UnitigNode {
    fn median_mapping_depth(&self) -> f64 {
        self.mapping_info.median_depth
    }
    fn mean_mapping_depth(&self) -> f64 {
        self.mapping_info.mean_depth
    }
    fn mapping_boundaries(&self) -> &Lapper<u32, SmallTwinOl> {
        &self.mapping_info.mapping_boundaries
    }
    fn set_mapping_info(&mut self, mapping_info: MappingInfo) {
        self.mapping_info = mapping_info;
    }
    fn mapping_info_present(&self) -> bool {
        self.mapping_info.present
    }
    fn reference_length(&self) -> usize {
        self.mapping_info.length
    }
    fn mapped_indices(&self) -> Vec<usize> {
        self.mapping_info.mapping_boundaries.iter().map(|x| x.val.query_id as usize).collect::<Vec<usize>>()
    }
}

impl NodeSequence for UnitigNode {
    fn base_seq(&self) -> &Seq<Dna> {
        &self.base_info.base_seq
    }
    fn read_positions_internal(&self) -> &Vec<(usize, usize)> {
        &self.base_info.read_positions_internal
    }
    fn left_cut(&self) -> usize {
        self.base_info.left_cut
    }
    fn right_cut(&self) -> usize {
        self.base_info.right_cut
    }
    fn length(&self) -> usize {
        self.base_info.length
    }
    fn base_info_present(&self) -> bool {
        self.base_info.present
    }
    fn set_info(&mut self, info: BaseInfo) {
        self.base_info = info;
    }
}

impl GraphNode for UnitigNode {
    fn in_edges(&self) -> &[EdgeIndex] {
        &self.in_edges
    }
    fn out_edges(&self) -> &[EdgeIndex] {
        &self.out_edges
    }
    fn in_edges_mut(&mut self) -> &mut Vec<EdgeIndex> {
        &mut self.in_edges
    }
    fn out_edges_mut(&mut self) -> &mut Vec<EdgeIndex> {
        &mut self.out_edges
    }
}

impl UnitigNode {}

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
    fn node1(&self) -> NodeIndex {
        self.from_unitig
    }
    fn node2(&self) -> NodeIndex {
        self.to_unitig
    }
    fn orientation1(&self) -> bool {
        self.f1
    }
    fn orientation2(&self) -> bool {
        self.f2
    }
}

impl<'a> GraphEdge for &'a UnitigEdge {
    fn node1(&self) -> NodeIndex {
        self.from_unitig
    }
    fn node2(&self) -> NodeIndex {
        self.to_unitig
    }
    fn orientation1(&self) -> bool {
        self.f1
    }
    fn orientation2(&self) -> bool {
        self.f2
    }
}

pub type UnitigGraph = BidirectedGraph<UnitigNode, UnitigEdge>;

impl UnitigGraph {
    fn re_unitig(&mut self) {
        let mut new_unitig_graph = UnitigGraph {
            nodes: FxHashMap::default(),
            edges: Vec::new(),
        };

        let mut old_nodes_to_new_nodes = FxHashMap::default();
        let mut old_terminal_edges: FxHashSet<NodeIndex> = FxHashSet::default();
        let mut old_nodes_new_ori: FxHashMap<NodeIndex, bool> = FxHashMap::default();
        let unitig_nbp = self.find_non_branching_paths();
        for (nodepath, edgepath) in unitig_nbp {
            // let node_id_path = nodepath
            //     .iter()
            //     .map(|&node_ind| self.nodes[&node_ind].read_indices_ori[0].0)
            //     .collect::<Vec<_>>();
            //eprintln!("Unitig path: {:?}, edgepath {:?}", node_id_path, edgepath);
            //eprintln!("First unit in len {:?} out len {:?} | last unit in len {:?} out len {:?}", 
            // self.nodes[&nodepath[0]].in_edges(), 
            // self.nodes[&nodepath[0]].out_edges(),
            // self.nodes[&nodepath[nodepath.len() - 1]].in_edges(),
            // self.nodes[&nodepath[nodepath.len() - 1]].out_edges());

            let mut new_read_indices_ori = vec![];
            let mut new_internal_overlaps = vec![];
            let mut new_read_names = vec![];
            let overlaps = edgepath
                .iter()
                .map(|&edge_idx| self.edges[edge_idx].as_ref().unwrap())
                .collect::<Vec<_>>();
            let mut new_unitig = UnitigNode::default();
            let unitig_path_ori = orientation_list(&nodepath, &overlaps);
            for (i, (&ori, &node_ind)) in unitig_path_ori.iter().zip(nodepath.iter()).enumerate() {
                old_nodes_to_new_nodes.insert(node_ind, new_unitig_graph.nodes.len());
                old_nodes_new_ori.insert(node_ind, ori);
                let unitig = self.nodes.get_mut(&node_ind).unwrap();
                let mut rdio = std::mem::take(&mut unitig.read_indices_ori);
                let internal_overlaps = std::mem::take(&mut unitig.internal_overlaps);
                let read_names = std::mem::take(&mut unitig.read_names);
                if ori {
                    new_read_indices_ori.extend(rdio);
                    new_internal_overlaps.extend(internal_overlaps);
                    new_read_names.extend(read_names);
                } else {
                    rdio = rdio
                        .into_iter()
                        .rev()
                        .map(|(ind, ori)| (ind, !ori))
                        .collect();
                    new_read_indices_ori.extend(rdio);
                    new_internal_overlaps.extend(internal_overlaps.into_iter().rev());
                    new_read_names.extend(read_names.into_iter().rev());
                }
                if i < nodepath.len() - 1 {
                    new_internal_overlaps
                        .push(self.edges[edgepath[i]].as_ref().unwrap().overlap.clone());
                }
                if let Some(mean_read_depth) = unitig.mean_read_depth {
                    let mut mean_read_depth_unitig = new_unitig.mean_read_depth.unwrap_or(0.);
                    mean_read_depth_unitig += mean_read_depth / nodepath.len() as f64;
                    new_unitig.mean_read_depth = Some(mean_read_depth_unitig);
                }
            }
            new_unitig.read_indices_ori = new_read_indices_ori;
            new_unitig.internal_overlaps = new_internal_overlaps;
            new_unitig.node_id = new_unitig.read_indices_ori[0].0;
            new_unitig.node_hash_id = new_unitig_graph.nodes.len();
            new_unitig.read_names = new_read_names;
            new_unitig_graph
                .nodes
                .insert(new_unitig_graph.nodes.len(), new_unitig);

            let last_node = &self.nodes[&nodepath[nodepath.len() - 1]];
            if last_node.in_edges().len() > 1 {
                old_terminal_edges.extend(last_node.in_edges());
            }
            if last_node.out_edges().len() > 1 {
                old_terminal_edges.extend(last_node.out_edges());
            }
            let first_node = &self.nodes[&nodepath[0]];
            if first_node.in_edges().len() > 1 {
                old_terminal_edges.extend(first_node.in_edges());
            }
            if first_node.out_edges().len() > 1 {
                old_terminal_edges.extend(first_node.out_edges());
            }

            //Must be single circular unitig
            if nodepath.len() == 1
                && first_node.in_edges().len() == 1
                && first_node.out_edges().len() == 1
                && first_node.in_edges()[0] == first_node.out_edges()[0]
            {
                old_terminal_edges.extend(first_node.in_edges());
            }

            //Check if circular unitig (more than 1 node); this condition does not imply
            // x -> o -> o -> o -> x
            if first_node.in_edges().len() == 1 &&
            first_node.out_edges().len() == 1 &&
            last_node.in_edges().len() == 1 &&
            last_node.out_edges().len() == 1 &&
            nodepath.len() > 1{
                //CHeck if there is an edge from the last node to the first node that isn't the same as the last edgepath edge
                for &edge_id in last_node.both_edges(){
                    if edge_id != edgepath[edgepath.len() - 1]{
                        let edge = self.edges[edge_id].as_ref().unwrap();
                        if edge.from_unitig == last_node.node_hash_id && edge.to_unitig == first_node.node_hash_id{
                            old_terminal_edges.insert(edge_id);
                        }
                        else if edge.from_unitig == first_node.node_hash_id && edge.to_unitig == last_node.node_hash_id{
                            old_terminal_edges.insert(edge_id);
                        }
                    }

                }
            }
        }

        //Re-edge the graph using the old unitig edges but mapped to the new nodes
        for old_edge_id in old_terminal_edges {
            let old_edge = self.edges[old_edge_id].as_ref().unwrap();
            let new_node1_ind = old_nodes_to_new_nodes[&old_edge.node1()];
            let new_node2_ind = old_nodes_to_new_nodes[&old_edge.node2()];
            let old_orientations = self.edges[old_edge_id]
                .as_ref()
                .unwrap()
                .get_orientation(old_edge.node1(), old_edge.node2());
            let new_orientations = (
                old_nodes_new_ori[&old_edge.node1()],
                old_nodes_new_ori[&old_edge.node2()],
            );

            let new_unitig_edge = UnitigEdge {
                overlap: self.edges[old_edge_id].as_ref().unwrap().overlap.clone(),
                from_read_idx: self.edges[old_edge_id].as_ref().unwrap().from_read_idx,
                to_read_idx: self.edges[old_edge_id].as_ref().unwrap().to_read_idx,
                from_unitig: new_node1_ind,
                to_unitig: new_node2_ind,
                f1: new_orientations.0 == old_orientations.0,
                f2: new_orientations.1 == old_orientations.1,
            };

            let new_n1 = new_unitig_graph.nodes.get_mut(&new_node1_ind).unwrap();

            if new_unitig_edge.f1 {
                new_n1.out_edges.push(new_unitig_graph.edges.len());
            } else {
                new_n1.in_edges.push(new_unitig_graph.edges.len());
            }

            let new_n2 = new_unitig_graph.nodes.get_mut(&new_node2_ind).unwrap();
            if new_unitig_edge.f2 {
                new_n2.in_edges.push(new_unitig_graph.edges.len());
            } else {
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

            let oris = orientation_list(&nodevec, &overlaps);

            let read_indices_ori: Vec<(usize, bool)> =
                (0..nodevec.len()).map(|i| (nodevec[i], oris[i])).collect();
            let node_id = read_indices_ori[0].0;
            let read_names = read_indices_ori
                .iter()
                .map(|(idx, _)| reads[*idx].id.clone())
                .collect();

            let depth = read_indices_ori.iter().map(|(idx, _)| reads[*idx].depth.unwrap()).sum::<f64>()
                / read_indices_ori.len() as f64;
            let unitig = UnitigNode {
                read_indices_ori,
                read_names,
                internal_overlaps: overlaps,
                in_edges: Vec::new(),
                out_edges: Vec::new(),
                node_id: node_id,
                node_hash_id: unitig_graph.nodes.len(),
                mean_read_depth: Some(depth),
                base_info: BaseInfo::default(),
                mapping_info: MappingInfo::default(),
            };
            unitig_graph.nodes.insert(unitig_graph.nodes.len(), unitig);
        }

        let mut terminal_edge_to_unitigs = FxHashMap::default();
        let mut unitigs_to_terminal_edges = FxHashMap::default();
        let mut edge_assignments = vec![];
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

        let mut seen_edges = FxHashSet::default();
        //Check for overlaps between terminal reads of unitigs
        for (unitig_id, term) in unitigs_to_terminal_edges {
            for (edge_id, read_id, read_ori) in term {
                if seen_edges.insert(edge_id) == false {
                    continue;
                }
                let other_unitigs = terminal_edge_to_unitigs.get(&edge_id).unwrap();
                for (unitig_id2, read_id2, read_ori2) in other_unitigs {
                    //TODO I think this should be disabled. Terminal edges should be able to look back onto itself. 
                    //if unitig_id == *unitig_id2 {
                    //    continue;
                    //}
                    let overlap_edge = overlap_graph.edges[*edge_id].as_ref().unwrap();
                    
                    //Ensure that the edge corresponds to the unitigs
                    if overlap_edge.node1() != read_id || overlap_edge.node2() != *read_id2 {
                        if overlap_edge.node1() != *read_id2 || overlap_edge.node2() != read_id {
                            continue;
                        }
                    }

                    //let sorted_uni_ids = if unitig_id < *unitig_id2 {
                    //    (*unitig_id, **unitig_id2)
                    //} else {
                    //    (**unitig_id2, *unitig_id)
                    //};

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

    pub fn to_fasta<T>(&self, filename: T, args: &Cli) 
    where 
        T: AsRef<std::path::Path>,
    {
        let mut bufwriter = BufWriter::new(std::fs::File::create(filename).unwrap());
        for (_id, unitig) in self.nodes.iter() {
            if unitig.read_indices_ori.len() < args.min_reads_contig
                && unitig.in_edges().len() + unitig.out_edges().len() == 0
                {
                    log::trace!(
                        "Unitig {} is disconnected with < 3 reads",
                        unitig.read_indices_ori[0].0
                    );
                    continue;
                }
            let mut base_seq = unitig.base_seq();
            let empty_seq = dna!("A").to_owned();
            if base_seq.len() == 0 {
                base_seq = &empty_seq;
            }
            let mut name = ">u".to_string() + &unitig.read_indices_ori[0].0.to_string();
            let depth = if unitig.mapping_info_present() {
                unitig.mean_mapping_depth()
            } else {
                unitig.mean_read_depth.unwrap_or(0.01)
            };
            name += &format!(
                " len:{} depth:{:.1}",
                unitig.length(),
                depth
            );
            if base_seq.len() == 1{
                name += " EMPTY";
            }
            writeln!(bufwriter, "{}", name).unwrap();
            writeln!(bufwriter, "{}", base_seq).unwrap();
        }
    }

    // Convert to GFA format
    pub fn to_gfa<T>(&mut self, filename: T, output_readgroups: bool, output_sequences: bool, reads: &[TwinRead], args: &Cli)
    where
        T: AsRef<std::path::Path>,
    {
        let edge_file = filename.as_ref().with_extension("edges");
        let mut bufwriter = BufWriter::new(std::fs::File::create(filename).unwrap());
        let mut edgewriter = BufWriter::new(std::fs::File::create(edge_file).unwrap());
        let mut gfa = String::new();
        // Header
        gfa.push_str("H\tVN:Z:1.0\n");

        // Segments
        for (_id, unitig) in self.nodes.iter_mut() {
            if unitig.read_indices_ori.len() < args.min_reads_contig
                && unitig.in_edges().len() + unitig.out_edges().len() == 0
            {
                log::trace!(
                    "Unitig {} is disconnected with < 3 reads",
                    unitig.read_indices_ori[0].0
                );
                continue;
            }

            if !unitig.base_info_present() {
                panic!("Unitig does not have base info");
            }

            let mut base_seq = unitig.base_seq();
            let empty_seq = Seq::new();
            if !output_sequences {
                base_seq = &empty_seq;
            } 

            let depth = if unitig.mapping_info_present() {
                unitig.mean_mapping_depth()
            } else {
                unitig.mean_read_depth.unwrap_or(0.01)
            };

            gfa.push_str(&format!(
                "S\tu{}\t{}\tLN:i:{}\tDP:f:{:.1}\n",
                unitig.read_indices_ori[0].0,
                //String::from_utf8(unitig.raw_consensus()).unwrap(),
                base_seq,
                unitig.length(),
                depth,
            ));

            if output_readgroups {
                let mut curr_pos = 0;
                let mut count = 0;
                for (read_idx, read_ori) in unitig.read_indices_ori.iter() {
                    let ori_string = if *read_ori { "+" } else { "-" };
                    let read = &reads[*read_idx];
                    let range = unitig.read_positions_internal()[count];
                    let length = range.1 - range.0;
                    gfa.push_str(&format!(
                        "a\tu{},{}-{}\t{}\t{}\t{}\t{}\t{}\n",
                        unitig.read_indices_ori[0].0,
                        range.0,
                        range.1,
                        read_idx,
                        curr_pos,
                        read.id,
                        ori_string,
                        length
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
            writeln!(
                edgewriter,
                "{} {} {} {} OL:{} SNP_SHARE:{} SNP_DIFF:{} READ1: {} {} READ2:{} {}",
                id1,
                from_orient,
                id2,
                to_orient,
                edge.overlap.overlap_len_bases,
                edge.overlap.shared_snpmers,
                edge.overlap.diff_snpmers,
                edge.overlap.node1,
                if edge.overlap.forward1 {"+"} else {"-"},
                edge.overlap.node2,
                if edge.overlap.forward2 {"+"} else {"-"},
            )
            .unwrap();

            gfa.push_str(&format!(
                "L\tu{}\t{}\tu{}\t{}\t{}M\n",
                id1,
                from_orient,
                id2,
                to_orient,
                //edge.overlap.overlap_len_bases
                0
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
                if e.from_unitig == e.to_unitig {
                    continue;
                }
                else if e.from_read_idx == unitig.read_indices_ori[0].0 {
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
                if e.from_unitig == e.to_unitig{
                    continue
                }
                else if e.from_read_idx == unitig.read_indices_ori[0].0 {
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

    fn cut_overlap_boundaries(&self) -> FxHashMap<NodeIndex, (usize, usize)> {
        let mut visited_nodes: FxHashSet<NodeIndex> = FxHashSet::default();
        let mut changes = vec![];
        for (id, _) in self.nodes.iter() {
            if visited_nodes.contains(id) {
                continue;
            }
            let mut explore_vec: Vec<(usize, Option<&UnitigEdge>)> = vec![(*id, None)];
            let mut visited_nodes_dir = FxHashSet::default();
            let mut visited_nodes_cc = FxHashSet::default();
            loop {
                if explore_vec.is_empty() {
                    break;
                }

                let (node, inc_edge) = explore_vec.pop().unwrap();
                if let Some(edge) = inc_edge {
                    let direction_into_node = edge.node_edge_direction(&node);
                    if visited_nodes_dir.contains(&(node, direction_into_node.clone())) {
                        continue;
                    }
                    visited_nodes_dir.insert((node, direction_into_node));
                } else {
                    visited_nodes_dir.insert((node, Direction::Incoming));
                    visited_nodes_dir.insert((node, Direction::Outgoing));
                }

                let (left_cut, right_cut) = self.cut_node_and_inc_edge(&node, inc_edge);
                changes.push((node, left_cut, right_cut));
                visited_nodes_cc.insert(node);
                let unitig = self.nodes.get(&node).unwrap();
                for edge in unitig.in_edges.iter().chain(unitig.out_edges.iter()) {
                    let edge = self.edges[*edge].as_ref().unwrap();
                    let other_node = if edge.from_unitig == node {
                        edge.to_unitig
                    } else {
                        edge.from_unitig
                    };
                    explore_vec.push((other_node, Some(edge)));
                }
            }

            visited_nodes.extend(&visited_nodes_cc);
        }
        let mut return_cut_map = FxHashMap::default();
        for (node, left_cut, right_cut) in changes {
            let lr_cut = return_cut_map.entry(node).or_insert((0, 0));
            lr_cut.0 += left_cut;
            lr_cut.1 += right_cut;
        }
        return return_cut_map;
    }

    pub fn remove_tips(&mut self, length: usize, num_reads: usize, keep: bool) {
        let node_to_sizeread_map = self.get_all_connected_components();
        let mut unitigs_to_remove = Vec::new();
        let mut debug_ids = vec![];
        // First pass: find all dead ends
        let dead_ends: Vec<NodeIndex> = self
            .nodes
            .keys()
            .filter(|&node_idx| {
                let unitig = &self.nodes[node_idx];
                unitig.in_edges.is_empty() || unitig.out_edges.is_empty()
            })
            .copied()
            .collect();

        for dead_end_ind in dead_ends {
            if let Some(unitig) = self.nodes.get(&dead_end_ind) {
                let (bp_size_cc, reads_in_cc) = node_to_sizeread_map[&dead_end_ind];
                if unitig.length() <= length.min(bp_size_cc / 10)
                    || unitig.read_indices_ori.len() <= num_reads.min(reads_in_cc / 10)
                {
                    unitigs_to_remove.push(dead_end_ind);
                    debug_ids.push(unitig.read_indices_ori[0].0);
                }
            }
        }

        log::info!("Removing {} tips", unitigs_to_remove.len());
        log::trace!("Unitigs to remove: {:?}", debug_ids);
        self.remove_nodes(&unitigs_to_remove, keep);

        self.re_unitig();
    }

    pub fn pop_bubbles(&mut self, max_length: usize) {
        let node_ids = self.nodes.keys().copied().collect::<Vec<_>>();
        let mut visited: FxHashSet<NodeIndex> = FxHashSet::default();
        let mut num_bubbles = 0;
        for n_id in node_ids {
            if visited.contains(&n_id) {
                continue;
            }
            if self.nodes[&n_id].in_edges().len() > 1 {
                let opt = self.get_bubble_remove_nodes(Direction::Incoming, n_id, max_length);
                if let Some((remove_nodes, _visited_nodes)) = opt {
                    num_bubbles += 1;
                    visited.extend(&remove_nodes);
                    self.remove_nodes(&remove_nodes, false);
                    //visited.extend(visited_nodes);
                }
            }

            //Need to check again because node may be removed.
            if visited.contains(&n_id) {
                continue;
            }
            if self.nodes[&n_id].out_edges().len() > 1 {
                let opt = self.get_bubble_remove_nodes(Direction::Outgoing, n_id, max_length);
                if let Some((remove_nodes, _visited_nodes)) = opt {
                    num_bubbles += 1;
                    visited.extend(&remove_nodes);
                    self.remove_nodes(&remove_nodes, false);
                    //visited.extend(visited_nodes);
                }
            }
        }
        log::info!("Removed {} bubbles at max length {}", num_bubbles, max_length);
        self.re_unitig();
    }

    fn get_bubble_remove_nodes(
        &self,
        right_direction: Direction,
        n_id: NodeIndex,
        max_length: usize,
    ) -> Option<(Vec<NodeIndex>, FxHashSet<NodeIndex>)> {
        let mut seen_vertices = FxHashSet::default();
        seen_vertices.insert(n_id);
        let mut stack = vec![(n_id, right_direction)];
        let mut depth_length: FxHashMap<NodeIndex, f64> = FxHashMap::default();
        let mut distances: FxHashMap<NodeIndex, usize> = FxHashMap::default();
        let mut num_outstanding = 0;
        let mut post_visit_edgecount = FxHashMap::default();
        let mut traceback = FxHashMap::default();
        traceback.insert(n_id, None);

        log::trace!(
            "Starting node bubble {} ({})",
            self.nodes[&n_id].read_indices_ori[0].0,
            n_id
        );

        while !stack.is_empty() {
            let v = stack.pop().unwrap();
            let prev_score = if depth_length.contains_key(&v.0) {
                depth_length[&v.0]
            } else {
                0.
            };
            let prev_dist = if distances.contains_key(&v.0) {
                distances[&v.0]
            } else {
                0
            };
            for edge_ind in self.nodes[&v.0].edges_direction(&v.1) {
                let edge = self.edges[*edge_ind].as_ref().unwrap();
                let other_node_ind = edge.other_node(v.0);
                let other_node = &self.nodes[&other_node_ind];
                let num_reads = other_node.read_indices_ori.len() as f64;
                seen_vertices.insert(other_node_ind);
                let other_node_leftdir = edge.node_edge_direction(&other_node_ind);
                // This is a problem for circular contigs... 
                // if other_node_ind == n_id {
                //     return None;
                // }
                if other_node_ind == n_id{ 
                    // The node comes back to its own direction
                    if other_node_leftdir == right_direction {
                        return None;
                    }
                } 
                if prev_dist > max_length {
                    return None;
                }
                if !distances.contains_key(&other_node_ind) {
                    let other_node_inc_dir = other_node.edges_direction(&other_node_leftdir);
                    post_visit_edgecount.insert(other_node_ind, other_node_inc_dir.len() as i64);
                    num_outstanding += 1;
                    distances.insert(other_node_ind, prev_dist + other_node.length());
                    traceback.insert(other_node_ind, Some(v.0));
                    depth_length.insert(
                        other_node_ind,
                        prev_score + other_node.mean_read_depth.unwrap() * num_reads,
                    );
                } else {
                    if prev_dist + other_node.length() < distances[&other_node_ind] {
                        distances.insert(other_node_ind, prev_dist + other_node.length());
                    }
                    if prev_score + other_node.mean_read_depth.unwrap() * num_reads
                        > depth_length[&other_node_ind]
                    {
                        depth_length.insert(
                            other_node_ind,
                            prev_score + other_node.mean_read_depth.unwrap() * num_reads,
                        );
                        traceback.insert(other_node_ind, Some(v.0));
                    }
                }
                post_visit_edgecount.insert(
                    other_node_ind,
                    post_visit_edgecount[&other_node_ind] as i64 - 1,
                );
                if post_visit_edgecount[&other_node_ind] == 0 {
                    num_outstanding -= 1;
                    if other_node
                        .edges_direction_reverse(&other_node_leftdir)
                        .len()
                        != 0
                        || (num_outstanding == 0 && stack.len() == 0)
                    {
                        stack.push((other_node_ind, other_node_leftdir.reverse()));
                    }
                }

                log::trace!(
                    "Other node: {:?} ({})  num_outstanding: {:?}, post_visit: {:?}",
                    other_node_ind,
                    other_node.read_indices_ori[0].0,
                    num_outstanding,
                    &post_visit_edgecount
                );
            }

            if stack.len() == 1 && num_outstanding == 0 {
                let right_dir = &stack[0].1;
                let end_node = &self.nodes[&stack[0].0];
                if end_node.edges_direction_reverse(&right_dir).len() > 1 {
                    let mut pushed_vertices = FxHashSet::default();
                    let start_node_id = self.nodes[&n_id].node_id;
                    let end_node_id = end_node.node_id;
                    log::debug!("Bubble found between {} and {}", start_node_id, end_node_id);
                    log::trace!("{:?}", &seen_vertices);
                    let mut path = vec![];
                    let mut debug_path = vec![];
                    //Traceback from the end node to the start node
                    let mut curr_node = stack[0].0;
                    while let Some(prev_node) = traceback[&curr_node] {
                        //Traceback could be circular
                        if pushed_vertices.contains(&curr_node) {
                            break;
                        }
                        path.push(curr_node);
                        debug_path.push(self.nodes[&curr_node].node_id);
                        pushed_vertices.insert(curr_node);
                        curr_node = prev_node;
                    }
                    path.push(curr_node);
                    debug_path.push(self.nodes[&curr_node].node_id);
                    let remove_vertices = seen_vertices
                        .difference(&path.iter().cloned().collect::<FxHashSet<_>>())
                        .cloned()
                        .collect();
                    log::debug!("Best path {:?}", &debug_path);
                    if !debug_path.contains(&start_node_id) || !debug_path.contains(&end_node_id) {
                        log::debug!("ERROR: Path between {} and {} does not contain start or end node", start_node_id, end_node_id);
                        return None;
                    }

                    let first_last_matching_c1 = debug_path.first().unwrap() == &start_node_id && debug_path.last().unwrap() == &end_node_id;
                    let first_last_matching_c2 = debug_path.first().unwrap() == &end_node_id && debug_path.last().unwrap() == &start_node_id;

                    if !first_last_matching_c1 && !first_last_matching_c2 {
                        log::debug!("ERROR: Path between {} and {} does not start and end at the correct nodes", start_node_id, end_node_id);
                        return None
                    }
                    return Some((remove_vertices, seen_vertices));
                }
            }
        }

        return None;
    }

    fn cut_node_and_inc_edge(&self, node: &NodeIndex, inc: Option<&UnitigEdge>) -> (usize, usize) {
        let cut_dir_edges;
        let direction;
        let unitig = self.nodes.get(node).unwrap();
        if unitig.in_edges.len() + unitig.out_edges.len() == 0 {
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
            //Circular contig
            let contig = &self.nodes.get(node).unwrap();
            if contig.in_edges().len() == 1 && contig.out_edges().len() == 1 {
                let in_edge = self.edges[contig.in_edges()[0]].as_ref().unwrap();
                let out_edge = self.edges[contig.out_edges()[0]].as_ref().unwrap();
                if in_edge.to_unitig == out_edge.from_unitig {
                    let overlap_length = in_edge.overlap.overlap_len_bases;
                    return (0, overlap_length);
                }
            }

            //Otherwise
            return (0, 0);
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
                let last_read = unitig.read_indices_ori[unitig.read_indices_ori.len() - 1].0;
                if edge.from_read_idx == first_read {
                    return (edge.overlap.overlap_len_bases, 0);
                } else if edge.to_read_idx == first_read {
                    return (edge.overlap.overlap_len_bases, 0);
                } else if edge.from_read_idx == last_read {
                    return (0, edge.overlap.overlap_len_bases);
                } else if edge.to_read_idx == last_read {
                    return (0, edge.overlap.overlap_len_bases);
                } else {
                    dbg!(edge, &unitig.read_indices_ori);
                    dbg!(&self.nodes[&edge.from_unitig]);
                    dbg!(&self.nodes[&edge.to_unitig]);
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
                component_length += unitig.length();
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

    pub fn cut_chimeric_regions(&mut self) {
        let mut breakpoints = FxHashMap::default();
        let mut remove_node_ids = vec![];
        for node in self.nodes.values() {
            let bp = mapping::cov_mapping_breakpoints(node);
            if !bp.is_empty() {
                breakpoints.insert(node.node_hash_id, bp);
            }
        }
        let mut redirected_edges = FxHashMap::default();
        log::debug!(
            "{} possible chimeric regions in unitig graph",
            breakpoints.len()
        );

        for (unitig_hash_id, breakpoints) in breakpoints {
            remove_node_ids.push(unitig_hash_id);
            redirected_edges.insert((unitig_hash_id, Direction::Incoming), None);
            redirected_edges.insert((unitig_hash_id, Direction::Outgoing), None);
            log::trace!(
                "Unitig: {:?} breakpoints: {:?}",
                unitig_hash_id,
                breakpoints
            );
            if self.nodes[&unitig_hash_id].read_indices_ori.len() == 1 {
                if self.nodes[&unitig_hash_id].read_names.len() == 0 {
                    panic!("{:?}", self.nodes[&unitig_hash_id]);
                }
                log::trace!("Cut read id: {}", self.nodes[&unitig_hash_id].read_names[0]);
                continue;
            }
            let redir = self.cut_chimeric_unitig(unitig_hash_id, breakpoints);
            redirected_edges.extend(redir);
        }

        self.redirect_and_remove(redirected_edges);
        self.remove_nodes(&remove_node_ids, false);
        self.re_unitig();
    }

    // Calculate the consensus sequence of all unitigs
    pub fn get_sequence_info(&mut self, reads: &[TwinRead], config: &GetSequenceInfoConfig) {
        let blunted = config.blunted;
        let dna_seq_info = config.dna_seq_info;
        let best_overlap_chunk = config.best_overlap_chunk;
        let cut_map;
        if blunted {
            cut_map = self.cut_overlap_boundaries();
        } else {
            cut_map = self
                .nodes
                .keys()
                .map(|k| (*k, (0, 0)))
                .collect::<FxHashMap<_, _>>();
        }
        let mut baseinfos = vec![];
        for (key, node) in self.nodes.iter() {
            let left_cut = cut_map[key].0;
            let right_cut = cut_map[key].1;
            let base_info;
            if best_overlap_chunk{
                base_info = get_base_info_mapchunks(left_cut, right_cut, node, reads, dna_seq_info)
            }
            else{
                base_info = get_base_info_overlaps(left_cut, right_cut, node, reads, dna_seq_info)
            }
            baseinfos.push((*key, base_info));
        }

        for (key, base_info) in baseinfos {
            self.nodes.get_mut(&key).unwrap().base_info = base_info;
        }
    }

    fn cut_chimeric_unitig(
        &mut self,
        unitig_hash_id: NodeIndex,
        breakpoints: Vec<Breakpoints>,
    ) -> FxHashMap<(NodeIndex, Direction), Option<NodeIndex>> {
        let mut cut_up_nodes = vec![];
        let num_nodes = self.nodes.len();
        let node = self.nodes.get(&unitig_hash_id).unwrap();
        log::trace!("breakpoints: {:?}", breakpoints);

        // Cut a multiple-read unitig with a chimera
        let mut cut_indices = vec![];
        let mut breakpoint_positions = breakpoints.iter().map(|bp| (bp.pos1, bp.pos2)).collect::<Vec<_>>();
        breakpoint_positions.sort();
        let mut cumulative = 0;
        let mut bpoint_ind = 0;
        for (l, (read_s, read_e)) in node.read_positions_internal().iter().enumerate() {
            cumulative += read_e - read_s;
            while cumulative > breakpoint_positions[bpoint_ind].0 {
                if !cut_indices.contains(&l) {
                    cut_indices.push(l);
                    log::trace!("Cut read id: {}", node.read_names[l]);
                }
                bpoint_ind += 1;
                if bpoint_ind == breakpoints.len() {
                    break;
                }
            }
            if bpoint_ind == breakpoints.len() {
                break;
            }
        }

        //cut_indices are ith read. internal overlaps represent read to read+1.
        let mut vec_new_read_indices_ori = vec![];
        let mut vec_new_internal_overlaps = vec![];
        let mut vec_new_indices = vec![];

        let mut readios = vec![];
        let mut ioverlaps = vec![];
        let mut indices = vec![];

        let mut last_break = 0;
        for i in 0..node.read_indices_ori.len() {
            if cut_indices.contains(&i) {
                last_break = i + 1;
                vec_new_read_indices_ori.push(readios);
                vec_new_internal_overlaps.push(ioverlaps);
                vec_new_indices.push(indices);
                readios = vec![];
                ioverlaps = vec![];
                indices = vec![];
                continue;
            }
            readios.push(node.read_indices_ori[i].clone());
            indices.push(i);
            if i > last_break {
                ioverlaps.push(node.internal_overlaps[i - 1].clone());
            }
        }

        vec_new_read_indices_ori.push(readios);
        vec_new_internal_overlaps.push(ioverlaps);
        vec_new_indices.push(indices);

        debug_assert!(vec_new_internal_overlaps.len() == vec_new_read_indices_ori.len());
        debug_assert!(vec_new_indices.len() == vec_new_read_indices_ori.len());
        let num_new_nodes = vec_new_read_indices_ori.len();

        let mut redirected_nodes = FxHashMap::default();
        let mut add_counter = 0;
        for (i, (rio, nio)) in vec_new_read_indices_ori
            .into_iter()
            .zip(vec_new_internal_overlaps.into_iter())
            .enumerate()
        {
            //Cut the read at the end of a unitig; remove entirely
            let new_node_hash_id = num_nodes + add_counter;
            let redirect = if rio.is_empty() {
                None
            } else {
                Some(new_node_hash_id)
            };
            if i == 0 {
                redirected_nodes.insert((unitig_hash_id, Direction::Incoming), redirect);
            }
            if i == num_new_nodes - 1 {
                redirected_nodes.insert((unitig_hash_id, Direction::Outgoing), redirect);
            }
            if rio.is_empty() {
                continue;
            }
            let new_node = UnitigNode {
                node_id: rio[0].0,
                read_indices_ori: rio,
                internal_overlaps: nio,
                read_names: vec_new_indices[i]
                    .iter()
                    .map(|&idx| self.nodes[&unitig_hash_id].read_names[idx].clone())
                    .collect(),
                in_edges: vec![],
                out_edges: vec![],
                node_hash_id: new_node_hash_id,
                mean_read_depth: self.nodes[&unitig_hash_id].mean_read_depth,
                base_info: BaseInfo::default(),
                mapping_info: MappingInfo::default(),
            };
            log::trace!(
                "New node: {}, {}, {}",
                new_node.node_id,
                new_node.node_hash_id,
                new_node.read_indices_ori.len()
            );
            cut_up_nodes.push(new_node);
            add_counter += 1;
        }
        for node in cut_up_nodes {
            if self.nodes.contains_key(&node.node_hash_id) {
                panic!("Node already exists");
            }
            self.nodes.insert(node.node_hash_id, node);
        }
        return redirected_nodes;
    }

    fn redirect_and_remove(
        &mut self,
        redirected_edges: FxHashMap<(NodeIndex, Direction), Option<NodeIndex>>,
    ) {
        let mut invalid_edges = FxHashSet::default();
        for (i, edge) in self.edges.iter_mut().enumerate() {
            if let Some(edge) = edge {
                let n1dir = edge.node_edge_direction(&edge.node1());
                let n2dir = edge.node_edge_direction(&edge.node2());

                if let Some(new_node_opt) = redirected_edges.get(&(edge.node1(), n1dir)) {
                    if let Some(new_node) = new_node_opt {
                        edge.from_unitig = *new_node;
                        if n1dir == Direction::Incoming {
                            self.nodes.get_mut(&new_node).unwrap().in_edges.push(i);
                        } else {
                            self.nodes.get_mut(&new_node).unwrap().out_edges.push(i);
                        }
                    } else {
                        invalid_edges.insert(i);
                    }
                }
                if let Some(new_node_opt) = redirected_edges.get(&(edge.node2(), n2dir)) {
                    if let Some(new_node) = new_node_opt {
                        edge.to_unitig = *new_node;
                        if n2dir == Direction::Incoming {
                            self.nodes.get_mut(&new_node).unwrap().in_edges.push(i);
                        } else {
                            self.nodes.get_mut(&new_node).unwrap().out_edges.push(i);
                        }
                    } else {
                        invalid_edges.insert(i);
                    }
                }
            }
        }

        for node in self.nodes.values_mut() {
            node.in_edges.retain(|x| !invalid_edges.contains(x));
            node.out_edges.retain(|x| !invalid_edges.contains(x));
        }
        self.edges.iter_mut().enumerate().for_each(|(i, edge)| {
            if invalid_edges.contains(&i) {
                *edge = None;
            }
        });

        for node in self.nodes.values() {
            for edge in node.in_edges.iter().chain(node.out_edges.iter()) {
                if self.edges[*edge].is_none() {
                    dbg!(edge);
                    dbg!(&node);
                    panic!("Edge is none");
                }
            }
        }

        for (key, _) in redirected_edges.keys() {
            self.nodes.remove(&key);
        }
    }

    
    fn _add_edges(&mut self, edges: Vec<UnitigEdge>) {
        for edge in edges {
            let f1 = edge.f1;
            let f2 = edge.f2;
            let new_edge_id = self.edges.len();
            let from_node = self.nodes.get_mut(&edge.from_unitig).unwrap();
            if f1 {
                from_node.out_edges.push(new_edge_id);
            } else {
                from_node.in_edges.push(new_edge_id);
            }
            let to_node = self.nodes.get_mut(&edge.to_unitig).unwrap();
            if f2 {
                to_node.in_edges.push(new_edge_id);
            } else {
                to_node.out_edges.push(new_edge_id);
            }
            let new_edge = Some(edge);
            self.edges.push(new_edge);
        }
    }

    pub fn cut_z_edges(&mut self, _args: &Cli){
        let mut edges_to_remove = FxHashSet::default();
        for unitig in self.nodes.values(){
            if unitig.in_edges().len() > 1{
                let z_edges = self.get_z_edges(unitig, Direction::Incoming);
                //TODO will handle case with > 1 z edge later 
                if z_edges.len() == 1{
                    for edge_id in z_edges{
                        edges_to_remove.insert(edge_id);
                    }
                }
            }
            if unitig.out_edges().len() > 1{
                let z_edges = self.get_z_edges(unitig, Direction::Outgoing);
                //TODO will handle case with > 1 z edge later 
                if z_edges.len() == 1{
                    for edge_id in z_edges{
                        edges_to_remove.insert(edge_id);
                    }
                }
            }
        }

        for edge in edges_to_remove.iter(){
            let edge = self.edges[*edge].as_ref().unwrap();
            log::trace!("Removing edge: {} -- {} ", self.nodes[&edge.from_unitig].node_id, self.nodes[&edge.to_unitig].node_id);
        }
        log::info!("Cut {} z-edges", edges_to_remove.len());
        self.remove_edges(edges_to_remove);
        self.re_unitig();
    }

    fn get_z_edges(&self, unitig: &UnitigNode, direction: Direction) -> Vec<EdgeIndex>{
        let edges;
        if direction == Direction::Outgoing{
            edges = unitig.out_edges()
        } else {
            edges = unitig.in_edges()
        }

        let mut junction_edges: Vec<EdgeIndex> = vec![];
        //let mut first_hit_nodes = FxHashSet::default();
        let mut used_edges = FxHashSet::default();
        let mut edges_to_explore: VecDeque<(EdgeIndex, &UnitigNode)> = VecDeque::new();
        for edge in edges{
            edges_to_explore.push_back((*edge, unitig));
            used_edges.insert(*edge);
        }
        junction_edges.extend(edges);
        while edges_to_explore.len() > 0{
            let (edge_id, new_unitig) = edges_to_explore.pop_front().unwrap();
            let edge = self.edges[edge_id].as_ref().unwrap();
            let other_node = edge.other_node(new_unitig.node_hash_id);
            let relative_direction_other = edge.node_edge_direction(&other_node);
            let other_edges = self.nodes[&other_node].edges_direction(&relative_direction_other);
            for &candidate_edge_id in other_edges{
                if !used_edges.contains(&candidate_edge_id){
                    let other_unitig = &self.nodes[&other_node];
                    edges_to_explore.push_back((candidate_edge_id, other_unitig));
                    junction_edges.push(candidate_edge_id);
                    used_edges.insert(candidate_edge_id);
                }
            }
        }

        let mut z_edges = vec![];

        for edge_id in junction_edges{
            let edge = &self.edges[edge_id].as_ref().unwrap();
            let n1 = &self.nodes[&edge.from_unitig];
            let n2 = &self.nodes[&edge.to_unitig];
            let n1_direction_edges = n1.edges_direction(&edge.node_edge_direction(&edge.from_unitig));
            let n2_direction_edges = n2.edges_direction(&edge.node_edge_direction(&edge.to_unitig));
            let mut n1_comp = FxHashSet::default();
            let mut n2_comp = FxHashSet::default();
            if n1_direction_edges.len() == 1 || n2_direction_edges.len() == 1{
                continue;
            }
            let edge_length = edge.overlap.overlap_len_bases;

            let mut exist_larger_edge = false; 
            for &edge_id2 in n1_direction_edges{
                if edge_id2 == edge_id{
                    continue;
                }
                let edge2 = self.edges[edge_id2].as_ref().unwrap();
                if edge2.overlap.overlap_len_bases > edge_length{
                    exist_larger_edge = true;
                }
                let other_node = edge2.other_node(n1.node_hash_id);
                let direction = edge2.node_edge_direction(&other_node);
                n1_comp.insert((other_node, direction));
            }

            for &edge_id2 in n2_direction_edges{
                if edge_id2 == edge_id{
                    continue;
                }
                let edge2 = self.edges[edge_id2].as_ref().unwrap();

                if edge2.overlap.overlap_len_bases > edge_length{
                    exist_larger_edge = true;
                }

                let other_node = edge2.other_node(n2.node_hash_id);
                let direction = edge2.node_edge_direction(&other_node);
                n2_comp.insert((other_node, direction));
            }

            // Need a edge that supports a z-edge to
            // be larger than the z-edge
            if !exist_larger_edge{
                continue;
            }

            let mut n1_connect_n2 = false;
            for (node, direction) in n1_comp.iter(){
                let unitig = &self.nodes[&node];
                //check if any edges hit anything in n2
                
                let edges = unitig.edges_direction(&direction);
                for &edge_id3 in edges{
                    let edge3 = self.edges[edge_id3].as_ref().unwrap();
                    let node3 = edge3.other_node(*node);
                    if n2_comp.contains(&(node3, edge3.node_edge_direction(&node3))){
                        //dbg!(n1.node_id, n2.node_id, &n1_comp, &n2_comp);
                        n1_connect_n2 = true;
                        break;
                    }
                }
            }

            if !n1_connect_n2{
                z_edges.push(edge_id);
            }
        }
        return z_edges
    }

    fn safe_given_forward_back(&self, unitig: &UnitigNode, edge: &UnitigEdge) -> bool{

        let max_length_forward_search = 100_000;
        let length_back_safe = 100_000;
        let mut safe = false;

        let mut seen_nodes = FxHashSet::default();
        let mut to_search = VecDeque::new();

        let other_node = edge.other_node(unitig.node_hash_id);
        let other_node_back_dir = edge.node_edge_direction(&other_node);

        //Node and back-direction
        to_search.push_back((other_node, other_node_back_dir));
        seen_nodes.insert(unitig.node_hash_id);
        let mut travelled = 0;

        let mut forbidden_nodes = FxHashSet::default();
        forbidden_nodes.insert(unitig.node_hash_id);

        //Search forward until a node has sufficient "backing"
        let mut num_iters = 0;
        while travelled < max_length_forward_search{
            if to_search.len() == 0{
                break;
            }
            let (node, direction) = to_search.pop_front().unwrap();
            //println!("DBG! {}", node);
            if seen_nodes.contains(&node){
                continue;
            }
            seen_nodes.insert(node);
            let unitig = &self.nodes[&node];

            //Check if the node has sufficient "backing"
            let back_found = self.search_dir_until_safe(unitig, direction, length_back_safe, &forbidden_nodes);
            if back_found{
                safe = true;
                break;
            }

            forbidden_nodes.insert(node);
            let forward_edges = unitig.edges_direction(&direction.reverse());
            for f_edge in forward_edges{
                let f_edge = self.edges[*f_edge].as_ref().unwrap();
                let f_node = f_edge.other_node(node);
                to_search.push_front((f_node, f_edge.node_edge_direction(&f_node)));
            }

            travelled += unitig.length();
            num_iters += 1;
        }

        return safe;
    }

    fn _forward_and_back(&self, unitig: &UnitigNode, direction: Direction) -> Vec<EdgeIndex>{
        let edges;
        if direction == Direction::Outgoing{
            edges = unitig.out_edges()
        } else {
            edges = unitig.in_edges()
        }

        let mut z_edges = vec![];
        let mut first_hit_nodes = FxHashSet::default();
        let mut used_edges = FxHashSet::default();
        for edge_id in edges{
            let edge = self.edges[*edge_id].as_ref().unwrap();
            let other_node = edge.other_node(unitig.node_hash_id);
            let relative_direction_other = edge.node_edge_direction(&other_node);
            first_hit_nodes.insert((other_node, relative_direction_other));
            used_edges.insert(*edge_id);
        }

        //Assume degree > 1 so  o -- x
        //                      |
        //                      x
        for edge_id in edges{
            let edge = self.edges[*edge_id].as_ref().unwrap();
            let other_node = edge.other_node(unitig.node_hash_id);
            let other_unitig = &self.nodes[&other_node];
            let relative_direction_other = edge.node_edge_direction(&other_node);
            let other_edges = other_unitig.edges_direction(&relative_direction_other);

            //Look for nodes that extend past the second node but not adjacent to first
            //      o -- o
            //      |
            // x -- o
            for edge_id2 in other_edges{
                if used_edges.contains(edge_id2){
                    continue;
                }
                let third_node = self.edges[*edge_id2].as_ref().unwrap().other_node(other_node);
                let third_unitig = &self.nodes[&third_node];
                let relative_direction_third = self.edges[*edge_id2].as_ref().unwrap().node_edge_direction(&third_node);
                let third_edges = third_unitig.edges_direction(&relative_direction_third);
                used_edges.insert(*edge_id2);

                // If the third node has >= 1 edge that does not connect to the first node 
                // or any of the second nodes,
                // then the second edge is a z-edge 
                let mut second_is_z = false;
                for edge_id3 in third_edges{
                    if used_edges.contains(edge_id3){
                        continue;
                    }
                    second_is_z = true;
                    let fourth_node = self.edges[*edge_id3].as_ref().unwrap().other_node(third_node);
                    let relative_direction_fourth = self.edges[*edge_id3].as_ref().unwrap().node_edge_direction(&fourth_node);
                    if fourth_node == unitig.node_hash_id || first_hit_nodes.contains(&(fourth_node, relative_direction_fourth)){
                        second_is_z = false;
                    }
                }
                if second_is_z{
                    z_edges.push(*edge_id2);
                }
            }
        }
        return z_edges;
    }

    pub fn resolve_bridged_repeats<T>(&mut self, _args: &Cli, ol_thresh: f64, out_file: T) 
    where T: AsRef<std::path::Path>,
    {
        let mut unitig_edge_file = BufWriter::new(std::fs::File::create(out_file).unwrap());
        let mut removed_edges = FxHashSet::default();
        for (_n_id, unitig) in self.nodes.iter() {
            for (_, edges) in [&unitig.in_edges(), &unitig.out_edges()]
                .into_iter()
                .enumerate()
            {
                let mut putative_removal = vec![];
                // try to cut
                if edges.len() > 1 {
                    let mut overlaps = vec![];
                    // let mut endpoint_covs = vec![];

                     for edge_id in edges.iter() {
                         let edge = &self.edges[*edge_id].as_ref().unwrap();
                         let ol = edge.overlap.overlap_len_bases;
                         overlaps.push(ol);
                    }

                    let max_ol = *overlaps.iter().max().unwrap();
                    for i in 0..edges.len() {
                        let edge = &self.edges[edges[i]].as_ref().unwrap();
                        let uni1 = &self.nodes[&edge.from_unitig];
                        let uni2 = &self.nodes[&edge.to_unitig];

                        let direction1 = edge.node_edge_direction(&edge.from_unitig);
                        let direction2 = edge.node_edge_direction(&edge.to_unitig);
                        let num_edges_1 = uni1.edges_direction(&direction1).len();
                        let num_edges_2 = uni2.edges_direction(&direction2).len();

                        //println!("DBG, {}-{}: {} {}", uni1.node_hash_id, uni2.node_hash_id, num_edges_1, num_edges_2);
                        let safe = self.safe_given_forward_back(unitig, edge);

                        let ol_score = overlaps[i] as f64 / max_ol as f64; // higher better
                        let mut smallest_less_pvalue = 1.;
                        for j in 0..edges.len() {
                            if i == j {
                                continue;
                            }
                            let other_edge = &self.edges[edges[j]].as_ref().unwrap();
                            let contingency_table = [
                                edge.overlap.shared_snpmers as u32,
                                other_edge.overlap.shared_snpmers as u32,
                                edge.overlap.diff_snpmers as u32,
                                other_edge.overlap.diff_snpmers as u32,
                            ];
                            let p_value = fishers_exact(&contingency_table).unwrap().less_pvalue;
                            if p_value < smallest_less_pvalue {
                                smallest_less_pvalue = p_value;
                            }
                        }

                        // Overlap is, coverage is twice, or the fisher pval is significantly smaller
                        // when compared against the most strain-specific best edges
                        let cut;
                        //let ol_thresh = if !safe { ol_thresh / 2. } else { ol_thresh};
                        if (ol_score < ol_thresh && safe)
                            || (smallest_less_pvalue < 0.05 && safe)
                        {
                            cut = true;
                            putative_removal.push(edges[i]);
                        } else {
                            cut = false;
                        }
                        writeln!(
                            unitig_edge_file,
                            "{}-{}, {} {} snp_share:{}, snp_diff:{}, ol_length:{}, ol_score:{}, specific_score:{}, safe:{}, removed:{}",
                            uni1.read_indices_ori[0].0,
                            uni2.read_indices_ori[0].0,
                            if edge.f1 {"+"} else {"-"},
                            if edge.f2 {"+"} else {"-"},
                            edge.overlap.shared_snpmers,
                            edge.overlap.diff_snpmers,
                            edge.overlap.overlap_len_bases,
                            ol_score,
                            smallest_less_pvalue,
                            safe,
                            cut
                        ).unwrap();
                    }
                }
                if putative_removal.len() < edges.len() {
                    for edge in putative_removal {
                        removed_edges.insert(edge);
                    }
                }
            }
        }
        log::info!("Cutting {} edges that have low relative confidence ({})", removed_edges.len(), ol_thresh);
        self.remove_edges(removed_edges);
        self.re_unitig();
    }

    pub fn print_statistics(self){
        //Print n50, name of largest contig, largest contig size, number of contigs

        let mut contig_sizes = self.nodes.iter().map(|(_, node)| node.length()).collect::<Vec<_>>();
        contig_sizes.sort();
        let n50 = contig_sizes.iter().sum::<usize>() / 2;
        let contig_sizes = contig_sizes.iter().rev();
        let mut curr_sum = 0;
        let mut n50_size = 0;
        for size in contig_sizes {
            curr_sum += size;
            if curr_sum >= n50 {
                n50_size = *size;
                break;
            }
        }

        let largest_contig = self.nodes.iter().max_by_key(|(_, node)| node.length());
        if let Some(largest_contig) = largest_contig {

            let largest_contig_size = largest_contig.1.length();
            let num_contigs = self.nodes.len();
            log::info!("-------------- Assembly statistics --------------");
            log::info!("N50: {}", n50_size);
            log::info!("Largest contig has size: {}", largest_contig_size);
            log::info!("Number of contigs: {}", num_contigs);
            log::info!("Total bases within assembly is {}", n50 * 2);
            log::info!("-------------------------------------------------");

        } else {
            log::info!("WARNING: No contigs found.");
            return;
        }
        
    }

    fn search_dir_until_safe(&self, unitig: &UnitigNode, direction: Direction, length: usize, forbidden_nodes: &FxHashSet<NodeIndex>) -> bool{
        let mut nodes_and_distances :FxHashMap<NodeIndex,usize> = FxHashMap::default();
        let mut to_search = VecDeque::new();
        let mut achieved_length = false;

        to_search.push_back((unitig.node_hash_id, direction));
        nodes_and_distances.insert(unitig.node_hash_id, 0);
        let mut iteration = 0;

        //Search in the direction until a path has been found with length > length
        while !to_search.is_empty(){
            let (node, dir) = to_search.pop_front().unwrap();
            //println!("DBG BACK {} DIST {}", node, nodes_and_distances[&node]);
            let edges = self.nodes[&node].edges_direction(&dir);
            let curr_min_dist = nodes_and_distances[&node];
            for edge_id in edges{
                let edge = self.edges[*edge_id].as_ref().unwrap();
                let other_node = edge.other_node(node);

                if forbidden_nodes.contains(&other_node){
                    continue;
                }

                let other_unitig = &self.nodes[&other_node];

                if nodes_and_distances.contains_key(&other_node){
                    *nodes_and_distances.get_mut(&other_node).unwrap() = nodes_and_distances[&other_node].min(curr_min_dist + other_unitig.length());
                }
                else{
                    nodes_and_distances.insert(other_node, curr_min_dist + other_unitig.length());
                    let other_dir = edge.node_edge_direction(&other_node).reverse();
                    to_search.push_back((other_node, other_dir));
                }

                if nodes_and_distances[&other_node] > length{
                    achieved_length = true;
                    break;
                }

                iteration += 1;
                if iteration > 100{
                    achieved_length = true;
                    break;
                }
                
            }
        }
        return achieved_length;
    }
}

fn overhang_and_overlap_list(node: &UnitigNode) -> (Vec<usize>, Vec<usize>) {
    let mut overhangs = Vec::new();
    let mut overlaps = Vec::new();
    for i in 0..node.read_indices_ori.len() - 1{
        let node1 = node.read_indices_ori[i].0;
        let node2 = node.read_indices_ori[i + 1].0;
        let internal_edge = &node.internal_overlaps[i];
        if node1 == internal_edge.node1 && node2 == internal_edge.node2{
            overhangs.push(internal_edge.hang1);
            overhangs.push(internal_edge.hang2);
            overlaps.push(internal_edge.overlap1_len);
            overlaps.push(internal_edge.overlap2_len);
        } else if node1 == internal_edge.node2 && node2 == internal_edge.node1{
            overhangs.push(internal_edge.hang2);
            overhangs.push(internal_edge.hang1);
            overlaps.push(internal_edge.overlap2_len);
            overlaps.push(internal_edge.overlap1_len);
        } else {
            panic!("Internal overlap does not match read indices");
        }
    }
    return (overhangs, overlaps);
}

fn get_base_info_overlaps(
    left_cut: usize,
    right_cut: usize,
    node: &UnitigNode,
    reads: &[TwinRead],
    dna_seq_info: bool,
) -> BaseInfo {
    let mut length = 0;
    let mut base_seq = Seq::new();
    let mut ranges = vec![];
    let mut carryover = left_cut;
    let (overhangs, internal_ol_len) = overhang_and_overlap_list(node);
    let mut hang_ind = 0;
    for (i, indori) in node.read_indices_ori.iter().enumerate() {
        let ind = indori.0;
        let ori = indori.1;
        let mut range;
        if node.read_indices_ori.len() == 1 {
            let end;
            if right_cut > reads[ind].base_length {
                end = 0;
            } else {
                end = reads[ind].base_length - right_cut;
            }
            range = (carryover, end);
        } else if i == 0 {
            let mut end = reads[ind].base_length as i64 + 1 - internal_ol_len[0] as i64 - overhangs[0] as i64 - overhangs[1] as i64;
            if end < 0 {
                end = 0;
            }
            range = (carryover, end as usize);
            hang_ind += 2;
        } else if i == node.read_indices_ori.len() - 1 {
            range = (carryover, reads[ind].base_length - right_cut.min(reads[ind].base_length));
            hang_ind += 2;
        } else {
            let mut end = reads[ind].base_length as i64 + 1 - internal_ol_len[hang_ind] as i64 - overhangs[hang_ind] as i64 - overhangs[hang_ind + 1] as i64;
            if end < 0 {
                end = 0;
            }
            range = (carryover, end as usize);
            hang_ind += 2;
        }
        if range.0 >= range.1 {
            carryover -= range.0 - range.1;
            range = (0, 0);
        } else {
            carryover = 0;
        }
        length += range.1 - range.0;
        ranges.push(range);
        if dna_seq_info {
            if ori {
                base_seq.append(&reads[ind].dna_seq[range.0..range.1]);
            } else {
                base_seq.append(&reads[ind].dna_seq.revcomp()[range.0..range.1]);
            }
        }
    }

    let base_info = BaseInfo {
        base_seq: base_seq,
        read_positions_internal: ranges,
        length: length,
        left_cut: left_cut,
        right_cut: right_cut,
        present: true,
    };
    return base_info;
}


fn get_base_info_mapchunks(
    _left_cut: usize,
    _right_cut: usize,
    _node: &UnitigNode,
    _reads: &[TwinRead],
    _dna_seq_info: bool,
) -> BaseInfo {
    return BaseInfo::default();
}
