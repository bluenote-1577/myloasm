use crate::cli::Cli;
use crate::constants::QUANTILE_UNITIG_WEIGHT;
use crate::graph::*;
use crate::twin_graph::*;
use crate::types::*;
use crate::unitig_utils::*;
use bio_seq::prelude::*;
use fishers_exact::fishers_exact;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use rust_lapper::Lapper;
use std::collections::VecDeque;
use std::io::BufWriter;
use std::io::Write;
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
    pub min_read_depth: Option<f64>, //Depth as measured by median of minimum depth of reads within unitig
    pub median_read_depth: Option<f64>, //Depth as measured by median of median depth of reads within unitig
    pub unique_length: Option<usize>, // Length of the unitig that is not covered by any other unitig's overlap; can be 0
    pub read_min_depths: Vec<(f64, usize)>,
    pub read_median_depths: Vec<(f64, usize)>,

    // If the unitig/contig is considered as an alternate; not implemented TODO
    pub alternate: bool,
    base_info: BaseInfo,
    mapping_info: MappingInfo,
}

pub trait NodeSequence {
    fn base_seq(&self) -> &Seq<Dna>;
    fn read_positions_internal(&self) -> &Vec<(usize, usize)>;
    fn left_cut(&self) -> usize;
    fn right_cut(&self) -> usize;
    fn cut_length(&self) -> usize;
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
        self.mapping_boundaries.intervals == other.mapping_boundaries.intervals
            && self.length == other.length
    }
}

impl Default for MappingInfo {
    fn default() -> Self {
        MappingInfo {
            median_depth: 0.0,
            minimum_depth: 0.0,
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
    fn min_mapping_depth(&self) -> f64 {
        self.mapping_info.minimum_depth
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
        self.mapping_info
            .mapping_boundaries
            .iter()
            .map(|x| x.val.query_id as usize)
            .collect::<Vec<usize>>()
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
    fn cut_length(&self) -> usize {
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
    pub fn new() -> Self {
        UnitigGraph {
            nodes: FxHashMap::default(),
            edges: Vec::new(),
        }
    }

    pub fn clear_edges(&mut self) {
        self.edges.clear();
        for (_, node) in self.nodes.iter_mut() {
            node.in_edges.clear();
            node.out_edges.clear();
        }
    }

    pub fn re_unitig(&mut self) {
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
            let mut median_read_depths = vec![];
            let mut min_read_depths = vec![];
            let overlaps = edgepath
                .iter()
                .map(|&edge_idx| self.edges[edge_idx].as_ref().unwrap())
                .collect::<Vec<_>>();
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

                min_read_depths.extend(unitig.read_min_depths.iter().map(|x| *x));
                median_read_depths.extend(unitig.read_median_depths.iter().map(|x| *x));
            }
            let median_min_depth = median_weight(&mut min_read_depths, QUANTILE_UNITIG_WEIGHT);
            let median_median_depth = median_weight(&mut median_read_depths, QUANTILE_UNITIG_WEIGHT);
            let new_unitig = UnitigNode {
                internal_overlaps: new_internal_overlaps,
                read_names: new_read_names,
                in_edges: Vec::new(),
                out_edges: Vec::new(),
                node_id: new_read_indices_ori[0].0,
                read_indices_ori: new_read_indices_ori,
                node_hash_id: new_unitig_graph.nodes.len(),
                min_read_depth: median_min_depth,
                median_read_depth: median_median_depth,
                base_info: BaseInfo::default(),
                mapping_info: MappingInfo::default(),
                read_min_depths: min_read_depths,
                read_median_depths: median_read_depths,
                unique_length: None,
                alternate: false,
            };
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
            if first_node.in_edges().len() == 1
                && first_node.out_edges().len() == 1
                && last_node.in_edges().len() == 1
                && last_node.out_edges().len() == 1
                && nodepath.len() > 1
            {
                //CHeck if there is an edge from the last node to the first node that isn't the same as the last edgepath edge
                for &edge_id in last_node.both_edges() {
                    if edge_id != edgepath[edgepath.len() - 1] {
                        let edge = self.edges[edge_id].as_ref().unwrap();
                        if edge.from_unitig == last_node.node_hash_id
                            && edge.to_unitig == first_node.node_hash_id
                        {
                            old_terminal_edges.insert(edge_id);
                        } else if edge.from_unitig == first_node.node_hash_id
                            && edge.to_unitig == last_node.node_hash_id
                        {
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
    pub fn from_overlaps(reads: &[TwinRead], overlaps: Vec<OverlapConfig>, args: &Cli) -> Self {
        let overlap_graph = read_graph_from_overlaps_twin(overlaps, args);
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

            let mut median_min_depth_vec = read_indices_ori
                .iter()
                .map(|(idx, _)| {
                    let read = &reads[*idx];
                    let depth = read.min_depth.unwrap();
                    let length = read.base_length;
                    (depth as f64, length)
                })
                .collect::<Vec<(f64, usize)>>();
            let median_min_depth = median_weight(&mut median_min_depth_vec, QUANTILE_UNITIG_WEIGHT);

            let mut median_median_depth_vec = read_indices_ori
                .iter()
                .map(|(idx, _)| {
                    let read = &reads[*idx];
                    let depth = read.median_depth.unwrap();
                    let length = read.base_length;
                    (depth as f64, length)
                })
                .collect::<Vec<(f64, usize)>>();
            let median_median_depth = median_weight(&mut median_median_depth_vec, QUANTILE_UNITIG_WEIGHT);

            let unitig = UnitigNode {
                read_indices_ori,
                read_names,
                internal_overlaps: overlaps,
                in_edges: Vec::new(),
                out_edges: Vec::new(),
                node_id: node_id,
                node_hash_id: unitig_graph.nodes.len(),
                min_read_depth: median_min_depth,
                median_read_depth: median_median_depth,
                base_info: BaseInfo::default(),
                mapping_info: MappingInfo::default(),
                read_median_depths: median_median_depth_vec,
                read_min_depths: median_min_depth_vec,
                unique_length: None,
                alternate: false,
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

        if cfg!(debug_assertions) {
            unitig_graph.test_consistent_left_right_edges();
        }

        let unitig_conf = GetSequenceInfoConfig {
            blunted: false,
            dna_seq_info: false,
            best_overlap_chunk: false,
        };
        unitig_graph.get_sequence_info(&reads, &unitig_conf);
        let out_file = std::path::Path::new(&args.output_dir).join("unitig_graph.gfa");
        unitig_graph.to_gfa(out_file, true, true, &reads, &args);

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
                unitig.median_mapping_depth()
            } else {
                unitig.median_read_depth.unwrap_or(-1.)
            };
            name += &format!(" len:{} depth:{:.1}", unitig.cut_length(), depth);
            if base_seq.len() == 1 {
                name += " EMPTY";
            }
            writeln!(bufwriter, "{}", name).unwrap();
            writeln!(bufwriter, "{}", base_seq).unwrap();
        }
    }

    // Convert to GFA format
    pub fn to_gfa<T>(
        &mut self,
        filename: T,
        output_readgroups: bool,
        output_sequences: bool,
        reads: &[TwinRead],
        args: &Cli,
    ) where
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
            if unitig.read_indices_ori.len() <= args.min_reads_contig
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

            let min_mapping_depth;
            let median_map_depth;
            if unitig.mapping_info_present() {
                median_map_depth = unitig.median_mapping_depth();
                min_mapping_depth = unitig.min_mapping_depth();
            } else {
                median_map_depth = -1.;
                min_mapping_depth = -1.;
            };

            let median_read_depth = unitig.median_read_depth.unwrap_or(-1.);
            let min_read_depth = unitig.min_read_depth.unwrap_or(-1.);

            gfa.push_str(&format!(
                "S\tu{}\t{}\tLN:i:{}\tDP:f:{:.1}\tMED_PRIM_DP:f:{:.1}\tMAP_DP:f:{:.1}\tMED_MAP_DP:f:{:.1}\n",
                unitig.read_indices_ori[0].0,
                //String::from_utf8(unitig.raw_consensus()).unwrap(),
                base_seq,
                unitig.cut_length(),
                min_read_depth,
                median_read_depth,
                min_mapping_depth,
                median_map_depth,
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
                        "a\tu{},{}-{}\t{}\t{}\t{}\t{}\t{}\tDP:{}\n",
                        unitig.read_indices_ori[0].0,
                        range.0,
                        range.1,
                        read_idx,
                        curr_pos,
                        read.id,
                        ori_string,
                        length,
                        read.min_depth.unwrap()
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
                "u{} {} u{} {} OL:{} SNP_SHARE:{} SNP_DIFF:{} READ1: {} {} READ2:{} {}",
                id1,
                from_orient,
                id2,
                to_orient,
                edge.overlap.overlap_len_bases,
                edge.overlap.shared_snpmers,
                edge.overlap.diff_snpmers,
                edge.overlap.node1,
                if edge.overlap.forward1 { "+" } else { "-" },
                edge.overlap.node2,
                if edge.overlap.forward2 { "+" } else { "-" },
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
                } else if e.from_read_idx == unitig.read_indices_ori[0].0 {
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
                if e.from_unitig == e.to_unitig {
                    continue;
                } else if e.from_read_idx == unitig.read_indices_ori[0].0 {
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

    fn remove_tips_internal(&mut self, length: usize, num_reads: usize, keep: bool) {
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
                if unitig.unique_length.unwrap() <= length.min(bp_size_cc / 10)
                    || unitig.read_indices_ori.len() <= num_reads.min(reads_in_cc / 10)
                {
                    //println!("Unitig {} is a dead end; unique_length {}, bp_size_cc {}, reads_in_cc {}", unitig.read_indices_ori[0].0, unitig.unique_length.unwrap(), bp_size_cc, reads_in_cc);
                    unitigs_to_remove.push(dead_end_ind);
                    debug_ids.push(unitig.read_indices_ori[0].0);
                }
            }
        }

        log::debug!("Removing {} tips", unitigs_to_remove.len());
        log::trace!("Unitigs to remove: {:?}", debug_ids);
        self.remove_nodes(&unitigs_to_remove, keep);
    }

    pub fn remove_tips(&mut self, length: usize, num_reads: usize, keep: bool) {
        self.remove_tips_internal(length, num_reads, keep);
        self.re_unitig();
    }

    pub fn pop_bubbles(&mut self, max_length: usize, max_number_nodes: Option<usize>, keep: bool) {
        let max_number_nodes = max_number_nodes.unwrap_or(usize::MAX);
        let node_ids = self.nodes.keys().copied().collect::<Vec<_>>();
        let mut visited: FxHashSet<NodeIndex> = FxHashSet::default();
        let mut num_bubbles = 0;
        for n_id in node_ids {
            if visited.contains(&n_id) {
                continue;
            }
            for direction in [Direction::Incoming, Direction::Outgoing].iter() {
                if visited.contains(&n_id) {
                    continue;
                }
                if self.nodes[&n_id].edges_direction(direction).len() > 1 {
                    if let Some(bubble_result) = self.double_bubble_remove_nodes(*direction, n_id, max_length, max_number_nodes){
                        visited.extend(&bubble_result.remove_nodes);
                        //TODO
                        self.remove_edges(bubble_result.remove_edges);
                        self.remove_nodes(&bubble_result.remove_nodes, keep);
                        num_bubbles += 1;
                    }
                }
            }
        }
        log::debug!(
            "Removed {} bubbles at max length {}",
            num_bubbles,
            max_length
        );
        self.re_unitig();
    }

    fn double_bubble_remove_nodes(
        &self, 
        direction: Direction, 
        n_id: NodeIndex,
        max_length: usize,
        max_number_nodes: usize
    ) -> Option<BubblePopResult> {

        let opt = self.get_bubble_remove_nodes(direction, n_id, max_length, max_number_nodes);
        if let Some(bubble_result) = opt {
            if let Some(bubble_result_back) = self.get_bubble_remove_nodes(
                bubble_result.end_direction,
                bubble_result.sink_hash_id,
                max_length,
                max_number_nodes
            ) {
                if bubble_result.remove_edges == bubble_result_back.remove_edges {
                    return Some(bubble_result);
                }
            }
        }

        return None
    }

    fn get_bubble_remove_nodes(
        &self,
        right_direction: Direction,
        n_id: NodeIndex,
        max_length: usize,
        max_number_nodes: usize
    ) -> Option<BubblePopResult> {
        let mut seen_vertices = FxHashSet::default();
        let mut seen_edges = FxHashSet::default();
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

            // Includes end and beginning node. 
            if distances.len() > max_number_nodes{
                return None;
            }
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
                seen_edges.insert(edge_ind);
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
                if other_node_ind == n_id {
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
                    distances.insert(
                        other_node_ind,
                        prev_dist + other_node.unique_length.unwrap(),
                    );
                    traceback.insert(other_node_ind, Some((v.0, edge_ind)));
                    depth_length.insert(
                        other_node_ind,
                        prev_score + other_node.min_read_depth.unwrap() * num_reads,
                    );
                } else {
                    if prev_dist + other_node.unique_length.unwrap() < distances[&other_node_ind] {
                        distances.insert(
                            other_node_ind,
                            prev_dist + other_node.unique_length.unwrap(),
                        );
                    }
                    if prev_score + other_node.min_read_depth.unwrap() * num_reads
                        > depth_length[&other_node_ind]
                    {
                        depth_length.insert(
                            other_node_ind,
                            prev_score + other_node.min_read_depth.unwrap() * num_reads,
                        );
                        traceback.insert(other_node_ind, Some((v.0, edge_ind)));
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
                    let mut good_edges = FxHashSet::default();
                    let mut pushed_vertices = FxHashSet::default();
                    let start_node_id = self.nodes[&n_id].node_id;
                    let start_node_hash_id = self.nodes[&n_id].node_hash_id;
                    let end_node_id = end_node.node_id;
                    let end_node_hash_id = end_node.node_hash_id;
                    log::debug!("Bubble found between {} and {}", start_node_id, end_node_id);
                    log::trace!("{:?}", &seen_vertices);
                    let mut path = vec![];
                    let mut debug_path = vec![];

                    //Traceback from the end node to the start node
                    let mut curr_node = stack[0].0;
                    while let Some((prev_node, good_edge)) = traceback[&curr_node] {
                        //Traceback could be circular
                        if pushed_vertices.contains(&curr_node) {
                            break;
                        }
                        path.push(curr_node);
                        good_edges.insert(good_edge);
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
                    let remove_edges = seen_edges
                        .difference(&good_edges)
                        .map(|x| **x)
                        .collect::<FxHashSet<_>>();

                    log::debug!("Best path {:?}", &debug_path);
                    if !debug_path.contains(&start_node_id) || !debug_path.contains(&end_node_id) {
                        log::debug!(
                            "ERROR: Path between {} and {} does not contain start or end node",
                            start_node_id,
                            end_node_id
                        );
                        return None;
                    }

                    let first_last_matching_c1 = debug_path.first().unwrap() == &start_node_id
                        && debug_path.last().unwrap() == &end_node_id;
                    let first_last_matching_c2 = debug_path.first().unwrap() == &end_node_id
                        && debug_path.last().unwrap() == &start_node_id;

                    if !first_last_matching_c1 && !first_last_matching_c2 {
                        log::trace!("Path between {} and {} does not start and end at the correct nodes", start_node_id, end_node_id);
                        return None;
                    }
                    return Some(
                        BubblePopResult{
                            original_direction: right_direction,
                            end_direction: stack[0].1.reverse(),
                            source_hash_id: start_node_hash_id,
                            sink_hash_id: end_node_hash_id,
                            remove_nodes: remove_vertices,
                            remove_edges: remove_edges,
                        }
                    );
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
                component_length += unitig.cut_length();
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
            if best_overlap_chunk {
                base_info = get_base_info_mapchunks(left_cut, right_cut, node, reads, dna_seq_info)
            } else {
                base_info = get_base_info_overlaps(left_cut, right_cut, node, reads, dna_seq_info)
            }
            baseinfos.push((*key, base_info));
        }

        for (key, base_info) in baseinfos {
            self.nodes.get_mut(&key).unwrap().base_info = base_info;
        }

        if !blunted {
            self.get_unique_lengths();
        }
    }

    fn _redirect_and_remove(
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

    pub fn cut_z_edges(&mut self, _args: &Cli) {
        let mut edges_to_remove = FxHashSet::default();
        for unitig in self.nodes.values() {
            if unitig.in_edges().len() > 1 {
                let z_edges = self.get_z_edges(unitig, Direction::Incoming);
                //TODO will handle case with > 1 z edge later
                if z_edges.len() == 1 {
                    for edge_id in z_edges {
                        edges_to_remove.insert(edge_id);
                    }
                }
            }
            if unitig.out_edges().len() > 1 {
                let z_edges = self.get_z_edges(unitig, Direction::Outgoing);
                //TODO will handle case with > 1 z edge later
                if z_edges.len() == 1 {
                    for edge_id in z_edges {
                        edges_to_remove.insert(edge_id);
                    }
                }
            }
        }

        for edge in edges_to_remove.iter() {
            let edge = self.edges[*edge].as_ref().unwrap();
            log::trace!(
                "Removing edge: {} -- {} ",
                self.nodes[&edge.from_unitig].node_id,
                self.nodes[&edge.to_unitig].node_id
            );
        }
        log::debug!("Cut {} z-edges", edges_to_remove.len());
        self.remove_edges(edges_to_remove);
        self.re_unitig();
    }

    fn get_z_edges(&self, unitig: &UnitigNode, direction: Direction) -> Vec<EdgeIndex> {
        let edges;
        if direction == Direction::Outgoing {
            edges = unitig.out_edges()
        } else {
            edges = unitig.in_edges()
        }

        let mut junction_edges: Vec<EdgeIndex> = vec![];
        //let mut first_hit_nodes = FxHashSet::default();
        let mut used_edges = FxHashSet::default();
        let mut edges_to_explore: VecDeque<(EdgeIndex, &UnitigNode)> = VecDeque::new();
        for edge in edges {
            edges_to_explore.push_back((*edge, unitig));
            used_edges.insert(*edge);
        }
        junction_edges.extend(edges);
        while edges_to_explore.len() > 0 {
            let (edge_id, new_unitig) = edges_to_explore.pop_front().unwrap();
            let edge = self.edges[edge_id].as_ref().unwrap();
            let other_node = edge.other_node(new_unitig.node_hash_id);
            let relative_direction_other = edge.node_edge_direction(&other_node);
            let other_edges = self.nodes[&other_node].edges_direction(&relative_direction_other);
            for &candidate_edge_id in other_edges {
                if !used_edges.contains(&candidate_edge_id) {
                    let other_unitig = &self.nodes[&other_node];
                    edges_to_explore.push_back((candidate_edge_id, other_unitig));
                    junction_edges.push(candidate_edge_id);
                    used_edges.insert(candidate_edge_id);
                }
            }
        }

        let mut z_edges = vec![];

        for edge_id in junction_edges {
            let edge = &self.edges[edge_id].as_ref().unwrap();
            let n1 = &self.nodes[&edge.from_unitig];
            let n2 = &self.nodes[&edge.to_unitig];
            let n1_direction_edges =
                n1.edges_direction(&edge.node_edge_direction(&edge.from_unitig));
            let n2_direction_edges = n2.edges_direction(&edge.node_edge_direction(&edge.to_unitig));
            let mut n1_comp = FxHashSet::default();
            let mut n2_comp = FxHashSet::default();
            if n1_direction_edges.len() == 1 || n2_direction_edges.len() == 1 {
                continue;
            }
            let edge_length = edge.overlap.overlap_len_bases;

            let mut exist_larger_edge = false;
            for &edge_id2 in n1_direction_edges {
                if edge_id2 == edge_id {
                    continue;
                }
                let edge2 = self.edges[edge_id2].as_ref().unwrap();
                if edge2.overlap.overlap_len_bases > edge_length {
                    exist_larger_edge = true;
                }
                let other_node = edge2.other_node(n1.node_hash_id);
                let direction = edge2.node_edge_direction(&other_node);
                n1_comp.insert((other_node, direction));
            }

            for &edge_id2 in n2_direction_edges {
                if edge_id2 == edge_id {
                    continue;
                }
                let edge2 = self.edges[edge_id2].as_ref().unwrap();

                if edge2.overlap.overlap_len_bases > edge_length {
                    exist_larger_edge = true;
                }

                let other_node = edge2.other_node(n2.node_hash_id);
                let direction = edge2.node_edge_direction(&other_node);
                n2_comp.insert((other_node, direction));
            }

            // Need a edge that supports a z-edge to
            // be larger than the z-edge
            if !exist_larger_edge {
                continue;
            }

            let mut n1_connect_n2 = false;
            for (node, direction) in n1_comp.iter() {
                let unitig = &self.nodes[&node];
                //check if any edges hit anything in n2

                let edges = unitig.edges_direction(&direction);
                for &edge_id3 in edges {
                    let edge3 = self.edges[edge_id3].as_ref().unwrap();
                    let node3 = edge3.other_node(*node);
                    if n2_comp.contains(&(node3, edge3.node_edge_direction(&node3))) {
                        //dbg!(n1.node_id, n2.node_id, &n1_comp, &n2_comp);
                        n1_connect_n2 = true;
                        break;
                    }
                }
            }

            if !n1_connect_n2 {
                z_edges.push(edge_id);
            }
        }
        return z_edges;
    }

    fn safe_given_forward_back(
        &self,
        unitig: &UnitigNode,
        edge: &UnitigEdge,
        max_forward: usize,
        max_reads_forward: usize,
        safe_length_back: usize,
        safety_cov_edge_ratio: Option<f64>,
        removed_edges: &FxHashSet<EdgeIndex>,
    ) -> bool {

        //Special case for circular contigs; we can cut the circular edge 
        // if the contig looks like 
        //                ->
        //     <-->(self)o 
        //                 ->                  
        if edge.from_unitig == edge.to_unitig {
            if self.nodes[&edge.from_unitig].in_edges().len() > 1
                && self.nodes[&edge.from_unitig].out_edges().len() > 1
            {
                return true;
            }
        }

        let max_length_forward_search = max_forward;
        let length_back_safe = safe_length_back;
        let mut safe = false;

        let mut seen_nodes = FxHashSet::default();
        let mut to_search = VecDeque::new();

        let other_node = edge.other_node(unitig.node_hash_id);
        let other_node_back_dir = edge.node_edge_direction(&other_node);

        //Node and back-direction
        to_search.push_back((other_node, other_node_back_dir));
        seen_nodes.insert(unitig.node_hash_id);
        let mut travelled = 0;
        let mut travelled_reads = 0;

        let mut forbidden_nodes = FxHashSet::default();
        forbidden_nodes.insert(unitig.node_hash_id);

        //Search forward until a node has sufficient "backing"
        while travelled < max_length_forward_search && travelled_reads < max_reads_forward {
            if to_search.len() == 0 {
                break;
            }
            let (node, direction) = to_search.pop_front().unwrap();
            //println!("DBG! {}", node);
            if seen_nodes.contains(&node) {
                continue;
            }
            seen_nodes.insert(node);
            let unitig = &self.nodes[&node];

            //Check if the node has sufficient "backing"
            let back_found =
                self.search_dir_until_safe(unitig, direction, length_back_safe, safety_cov_edge_ratio, &forbidden_nodes, removed_edges);
            if back_found {
                safe = true;
                break;
            }

            forbidden_nodes.insert(node);
            let forward_edges = self.get_safe_edges_from_cov_threshold(unitig.edges_direction(&direction.reverse()), safety_cov_edge_ratio);
            for f_edge in forward_edges {
                if removed_edges.contains(&f_edge) {
                    continue;
                }
                let f_edge = self.edges[f_edge].as_ref().unwrap();
                let f_node = f_edge.other_node(node);
                to_search.push_front((f_node, f_edge.node_edge_direction(&f_node)));
            }

            travelled_reads += unitig.read_indices_ori.len();
            travelled += unitig.unique_length.unwrap();
        }

        return safe;
    }

    fn _forward_and_back(&self, unitig: &UnitigNode, direction: Direction) -> Vec<EdgeIndex> {
        let edges;
        if direction == Direction::Outgoing {
            edges = unitig.out_edges()
        } else {
            edges = unitig.in_edges()
        }

        let mut z_edges = vec![];
        let mut first_hit_nodes = FxHashSet::default();
        let mut used_edges = FxHashSet::default();
        for edge_id in edges {
            let edge = self.edges[*edge_id].as_ref().unwrap();
            let other_node = edge.other_node(unitig.node_hash_id);
            let relative_direction_other = edge.node_edge_direction(&other_node);
            first_hit_nodes.insert((other_node, relative_direction_other));
            used_edges.insert(*edge_id);
        }

        //Assume degree > 1 so  o -- x
        //                      |
        //                      x
        for edge_id in edges {
            let edge = self.edges[*edge_id].as_ref().unwrap();
            let other_node = edge.other_node(unitig.node_hash_id);
            let other_unitig = &self.nodes[&other_node];
            let relative_direction_other = edge.node_edge_direction(&other_node);
            let other_edges = other_unitig.edges_direction(&relative_direction_other);

            //Look for nodes that extend past the second node but not adjacent to first
            //      o -- o
            //      |
            // x -- o
            for edge_id2 in other_edges {
                if used_edges.contains(edge_id2) {
                    continue;
                }
                let third_node = self.edges[*edge_id2]
                    .as_ref()
                    .unwrap()
                    .other_node(other_node);
                let third_unitig = &self.nodes[&third_node];
                let relative_direction_third = self.edges[*edge_id2]
                    .as_ref()
                    .unwrap()
                    .node_edge_direction(&third_node);
                let third_edges = third_unitig.edges_direction(&relative_direction_third);
                used_edges.insert(*edge_id2);

                // If the third node has >= 1 edge that does not connect to the first node
                // or any of the second nodes,
                // then the second edge is a z-edge
                let mut second_is_z = false;
                for edge_id3 in third_edges {
                    if used_edges.contains(edge_id3) {
                        continue;
                    }
                    second_is_z = true;
                    let fourth_node = self.edges[*edge_id3]
                        .as_ref()
                        .unwrap()
                        .other_node(third_node);
                    let relative_direction_fourth = self.edges[*edge_id3]
                        .as_ref()
                        .unwrap()
                        .node_edge_direction(&fourth_node);
                    if fourth_node == unitig.node_hash_id
                        || first_hit_nodes.contains(&(fourth_node, relative_direction_fourth))
                    {
                        second_is_z = false;
                    }
                }
                if second_is_z {
                    z_edges.push(*edge_id2);
                }
            }
        }
        return z_edges;
    }


    pub fn resolve_bridged_repeats<T>(
        &mut self,
        _args: &Cli,
        ol_thresh: f64,
        unitig_cov_ratio_cut: Option<f64>,
        safety_cov_edge_ratio: Option<f64>,
        out_file: T,
        max_forward: usize,
        max_reads_forward: usize,
        safe_length_back: usize,
    ) where
        T: AsRef<std::path::Path>,
    {
        let mut unitig_edge_file = BufWriter::new(std::fs::File::create(out_file).unwrap());
        let mut removed_edges = FxHashSet::default();
        let mut all_edges_sorted = vec![];

        for (i,edge) in self.edges.iter().enumerate(){
            if let Some(edge) = edge{
                let uni1 = &self.nodes[&edge.from_unitig];
                let uni2 = &self.nodes[&edge.to_unitig];
                let overlap_len = edge.overlap.overlap_len_bases;

                let cov_ratio = 1./pseudocount_cov(uni1.min_read_depth.unwrap(), uni2.min_read_depth.unwrap());
                let score = cov_ratio * overlap_len as f64;
                all_edges_sorted.push((i, score));
                // smaller score is worse. small overlap len and large diff
            }
        }

        all_edges_sorted.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let all_edges_sorted = all_edges_sorted.iter().map(|x| x.0).collect::<Vec<_>>();

        for edge_id in all_edges_sorted {
            self.safely_cut_edge(
                edge_id,
                &mut removed_edges,
                ol_thresh,
                unitig_cov_ratio_cut,
                safety_cov_edge_ratio,
                max_forward,
                max_reads_forward,
                safe_length_back,
                &mut unitig_edge_file,
            );
        }
        log::debug!(
            "Cutting {} edges that have low relative confidence ({})",
            removed_edges.len(),
            ol_thresh
        );
        self.remove_edges(removed_edges);
        self.re_unitig();
    }

    pub fn print_statistics(self, args: &Cli) {
        //Print n50, name of largest contig, largest contig size, number of contigs

        let mut contig_sizes = self
            .nodes
            .iter()
            .filter(|(_, node)| {
                    if node.read_indices_ori.len() <= args.min_reads_contig
                        && node.in_edges().len() + node.out_edges().len() == 0{
                            return false;
                    }
                    else{
                        return true;
                    }
                }
            )
            .map(|(_, node)| node.cut_length())
            .collect::<Vec<_>>();
        contig_sizes.sort();
        let n50 = contig_sizes.iter().sum::<usize>() / 2;
        let contig_sizes_iterev = contig_sizes.iter().rev();
        let mut curr_sum = 0;
        let mut n50_size = 0;
        for size in contig_sizes_iterev{
            curr_sum += size;
            if curr_sum >= n50 {
                n50_size = *size;
                break;
            }
        }

        let largest_contig = self.nodes.iter().max_by_key(|(_, node)| node.cut_length());
        if let Some(largest_contig) = largest_contig {
            let largest_contig_size = largest_contig.1.cut_length();
            let num_contigs = contig_sizes.len();
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

    // Search a unitig in Direction until a path has been found with length > length
    // subject to the constraint that the path does not contain any of the forbidden nodes or already-removed edges
    fn search_dir_until_safe(
        &self,
        unitig: &UnitigNode,
        direction: Direction,
        length: usize,
        safety_cov_edge_ratio: Option<f64>, 
        forbidden_nodes: &FxHashSet<NodeIndex>,
        removed_edges: &FxHashSet<EdgeIndex>
    ) -> bool {
        let mut nodes_and_distances: FxHashMap<NodeIndex, usize> = FxHashMap::default();
        let mut to_search = VecDeque::new();
        let mut achieved_length = false;

        //Push_front for dfs
        to_search.push_front((unitig.node_hash_id, direction));
        nodes_and_distances.insert(unitig.node_hash_id, 0);
        let mut iteration = 0;

        //Search in the direction until a path has been found with length > length
        while !to_search.is_empty() {
            let (node, dir) = to_search.pop_front().unwrap();
            //println!("DBG BACK {} DIST {}", node, nodes_and_distances[&node]);
            let edges = self.get_safe_edges_from_cov_threshold(self.nodes[&node].edges_direction(&dir), safety_cov_edge_ratio);
            let curr_min_dist = nodes_and_distances[&node];
            for edge_id in edges {
                let edge = self.edges[edge_id].as_ref().unwrap();
                let other_node = edge.other_node(node);

                if forbidden_nodes.contains(&other_node) || removed_edges.contains(&edge_id) {
                    continue;
                }

                let other_unitig = &self.nodes[&other_node];

                if nodes_and_distances.contains_key(&other_node) {
                    *nodes_and_distances.get_mut(&other_node).unwrap() = nodes_and_distances
                        [&other_node]
                        .min(curr_min_dist + other_unitig.unique_length.unwrap());
                } else {
                    nodes_and_distances.insert(
                        other_node,
                        curr_min_dist + other_unitig.unique_length.unwrap(),
                    );
                    let other_dir = edge.node_edge_direction(&other_node).reverse();
                    to_search.push_back((other_node, other_dir));
                }

                if nodes_and_distances[&other_node] > length {
                    achieved_length = true;
                    break;
                }

                iteration += 1;
                if iteration > 100 {
                    achieved_length = true;
                    break;
                }
            }
        }
        return achieved_length;
    }

    //Assume that the graph has not been blunted.
    fn get_unique_lengths(&mut self) {
        for node in self.nodes.values_mut() {
            if node.right_cut() > 0 || node.left_cut() > 0 {
                panic!("Node has been blunted. Cannot calculate unique length.");
            }
            let non_self_edges_in = node.in_edges().iter().filter(|&x| {
                let edge = self.edges[*x].as_ref().unwrap();
                edge.from_unitig != edge.to_unitig
            });
            let non_self_edges_out = node.out_edges().iter().filter(|&x| {
                let edge = self.edges[*x].as_ref().unwrap();
                edge.from_unitig != edge.to_unitig
            });

            let in_max_overlap = non_self_edges_in
                .map(|x| self.edges[*x].as_ref().unwrap().overlap.overlap_len_bases)
                .max()
                .unwrap_or(0);
            let out_max_overlap = non_self_edges_out
                .map(|x| self.edges[*x].as_ref().unwrap().overlap.overlap_len_bases)
                .max()
                .unwrap_or(0);
            node.unique_length = Some(
                (node.cut_length() as i64 - in_max_overlap.max(out_max_overlap) as i64).max(0)
                    as usize,
            );
        }
    }

    pub fn cut_coverage(&mut self, cov : f64){
        let mut nodes_to_remove = FxHashSet::default();
        for node in self.nodes.values(){
            if node.min_read_depth.unwrap() < cov{
                nodes_to_remove.insert(node.node_hash_id);
            }
        }
        self.remove_nodes(&nodes_to_remove.into_iter().collect::<Vec<_>>(), false);
        self.re_unitig();
    }

    pub fn progressive_cut_lowest(&mut self) -> usize{
        let mut lowest_cov = f64::MAX;
        for node in self.nodes.values(){
            if node.both_edges().collect::<Vec<_>>().len() == 0{
                continue;
            }
            //circular
            if node.in_edges == node.out_edges{
                continue;
            }
            if node.min_read_depth.unwrap() < lowest_cov{
                lowest_cov = node.min_read_depth.unwrap();
            }
        }
        let mut nodes_to_remove = FxHashSet::default();
        for node in self.nodes.values(){
            if node.min_read_depth.unwrap() == lowest_cov{
                if node.in_edges == node.out_edges{
                    continue;
                }
                nodes_to_remove.insert(node.node_hash_id);
            }
        }
        let num_removed_nodes = nodes_to_remove.len();

        self.remove_nodes(&nodes_to_remove.into_iter().collect::<Vec<_>>(), true);
        self.re_unitig();
        return num_removed_nodes;
    }

    fn safely_cut_edge<T>(
        &self,
        edge_id: EdgeIndex,
        removed_edges: &mut FxHashSet<EdgeIndex>,
        ol_thresh: f64,
        unitig_cov_ratio_cut: Option<f64>,
        safety_cov_edge_ratio: Option<f64>,
        max_forward: usize,
        max_reads_forward: usize,
        safe_length_back: usize,
        unitig_edge_file: &mut BufWriter<T>,
    ) where
        T: Write,
    {
        let edge = &self.edges[edge_id].as_ref().unwrap();
        let unitig_terminals = [edge.from_unitig, edge.to_unitig];

        for unitig_id in unitig_terminals.iter() {
            if removed_edges.contains(&edge_id) {
                continue;
            }

            let unitig = &self.nodes[unitig_id];
            let direction = edge.node_edge_direction(unitig_id);
            let edges = unitig.edges_direction(&direction);

            let non_cut_edges = edges
                .iter()
                .filter(|&x| !removed_edges.contains(x))
                .map(|x| *x)
                .collect::<Vec<_>>();

            if non_cut_edges.len() <= 1 {
                continue;
            }

            let mut overlaps = vec![];
            for &edge_id in non_cut_edges.iter() {
                let edge = &self.edges[edge_id].as_ref().unwrap();
                let ol = edge.overlap.overlap_len_bases;
                overlaps.push(ol);
            }

            let max_ol = *overlaps.iter().max().unwrap();

            let uni1 = &self.nodes[&edge.from_unitig];
            let uni2 = &self.nodes[&edge.to_unitig];

            let safe = self.safe_given_forward_back(
                unitig,
                edge,
                max_forward,
                max_reads_forward,
                safe_length_back,
                safety_cov_edge_ratio,
                removed_edges,
            );

            if !safe {
                continue;
            }

            let ol_score = edge.overlap.overlap_len_bases as f64 / max_ol as f64; // higher better

            let cov_pseudo_val = pseudocount_cov(
                uni1.min_read_depth.unwrap(),
                uni2.min_read_depth.unwrap(),
            );

            let cut;
            if (ol_score < ol_thresh)
                || (unitig_cov_ratio_cut.is_some() && cov_pseudo_val > unitig_cov_ratio_cut.unwrap())
            {
                cut = true;
                removed_edges.insert(edge_id);
            } else {
                cut = false;
            }
            writeln!(
                unitig_edge_file,
                "u{}-u{}, cut:{} safe:{} snp_share:{}, snp_diff:{}, ol_length:{}, ol_score:{}, specific_score:{}",
                uni1.node_id,
                uni2.node_id,
                cut,
                safe,
                edge.overlap.shared_snpmers,
                edge.overlap.diff_snpmers,
                edge.overlap.overlap_len_bases,
                ol_score,
                cov_pseudo_val
            )
            .unwrap();
        }
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    // Helper struct to make test graph construction easier
    struct MockUnitigBuilder {
        nodes: FxHashMap<NodeIndex, UnitigNode>,
        edges: Vec<Option<UnitigEdge>>,
        next_node_id: NodeIndex,
        corresponding_reads: Vec<TwinRead>,
    }

    impl MockUnitigBuilder {
        fn new() -> Self {
            Self {
                nodes: FxHashMap::default(),
                edges: Vec::new(),
                next_node_id: 0,
                corresponding_reads: Vec::new(),
            }
        }

        // Add a node with specified number of reads and coverage
        fn add_node(&mut self, num_reads: usize, min_depth: f64) -> NodeIndex {

            let mut internal_overlaps = vec![];
            let read_indices: Vec<(NodeIndex, bool)> = (0..num_reads)
                .map(|i| (i as NodeIndex + self.corresponding_reads.len(), true))
                .collect();

            for i in 0..num_reads - 1{
                let generic_internal_overlap = ReadOverlapEdgeTwin {
                    node1: self.corresponding_reads.len() as NodeIndex + i,
                    node2: self.corresponding_reads.len() as NodeIndex + 1 + i,
                    hang1: 0,
                    hang2: 0,
                    overlap1_len: 1000,
                    overlap2_len: 1000,
                    forward1: true,
                    forward2: true,
                    overlap_len_bases: 1000,
                    shared_minimizers: 100,
                    diff_snpmers: 0,
                    shared_snpmers: 10,
                };
                internal_overlaps.push(generic_internal_overlap);
            }

            for _ in 0..num_reads{
                let generic_read = TwinRead{
                    minimizers: vec![],
                    snpmers: vec![],
                    dna_seq: Seq::new(),
                    base_length: 2000,
                    k: 21,
                    est_id: None,
                    id: "na".to_string(),
                    min_depth: Some(min_depth),
                    median_depth: Some(2.*min_depth),
                    split_chimera: false,
                };
                
                self.corresponding_reads.push(generic_read);
            }

            let node_id = self.next_node_id;
            self.next_node_id += 1;

            // Create mock read indices
            let node = UnitigNode {
                read_indices_ori: read_indices,
                //internal overlaps only used for base-level information... not needed for topology tests
                internal_overlaps,
                read_names: vec!["na".to_string(); num_reads],
                in_edges: vec![],
                out_edges: vec![],
                node_id: node_id,
                node_hash_id: node_id,
                min_read_depth: Some(min_depth),
                median_read_depth: Some(2.*min_depth),
                unique_length: None,
                base_info: BaseInfo::default(),
                mapping_info: MappingInfo::default(),
                read_min_depths: vec![(min_depth, 1000);num_reads],
                read_median_depths: vec![(min_depth, 1000); num_reads],
                alternate: false,
            };

            self.nodes.insert(node_id, node);
            node_id
        }

        // Add an edge between nodes
        fn add_edge(&mut self, from: NodeIndex, to: NodeIndex, overlap_len: usize, f1: bool, f2: bool) {
            let from_read_idx = self.nodes[&from].read_indices_ori.last().unwrap().0;
            let to_read_idx = self.nodes[&to].read_indices_ori.first().unwrap().0;
            let edge = UnitigEdge {
                overlap: ReadOverlapEdgeTwin {
                    node1: from_read_idx,
                    node2: to_read_idx,
                    hang1: 0,
                    hang2: 0,
                    overlap1_len: overlap_len,
                    overlap2_len: overlap_len,
                    forward1: f1,
                    forward2: f2,
                    overlap_len_bases: overlap_len,
                    shared_minimizers: 100,
                    diff_snpmers: 0,
                    shared_snpmers: 10,
                },
                from_read_idx,
                to_read_idx,
                from_unitig: from,
                to_unitig: to,
                f1,
                f2,
            };

            let edge_idx = self.edges.len();
            
            // Update node edge lists
            if f1 {
                self.nodes.get_mut(&from).unwrap().out_edges.push(edge_idx);
            } else {
                self.nodes.get_mut(&from).unwrap().in_edges.push(edge_idx);
            }
            
            if f2 {
                self.nodes.get_mut(&to).unwrap().in_edges.push(edge_idx);
            } else {
                self.nodes.get_mut(&to).unwrap().out_edges.push(edge_idx);
            }

            self.edges.push(Some(edge));
        }

        fn build(self) -> (UnitigGraph, Vec<TwinRead>) {
           let mut g = UnitigGraph {
                nodes: self.nodes,
                edges: self.edges,
            };
            let reads = self.corresponding_reads;
            g.get_sequence_info(&reads, &GetSequenceInfoConfig::default());
            return (g, reads);
        }
    }

    // Test tip removal
    #[test]
    fn test_remove_tips() {
        // Create a graph with a tip:
        // n1 -> n2 
        //  ↘
        //   n4 (tip)
        let mut builder = MockUnitigBuilder::new();
        
        let n1 = builder.add_node(10, 10.0);
        let n2 = builder.add_node(10, 10.0);
        let n4 = builder.add_node(1, 5.0); // tip node

        builder.add_edge(n1, n2, 100, true, true);
        builder.add_edge(n1, n4, 100, true, true);

        let (mut graph, _reads) = builder.build();
        assert!(graph.nodes.len() == 3);
        
        // Remove tips
        graph.remove_tips_internal(500, 2, false);
        assert!(graph.nodes.len() == 2);
        assert!(!graph.nodes.contains_key(&n4));
    }

    // Test bubble detection
    #[test]
    fn test_bubble_detection_simple() {
        // Create a simple bubble:
        // n0 -> n1 -> n3
        //  ↘        ↗  
        //   -> n2 ->
        let mut builder = MockUnitigBuilder::new();
        
        let n0 = builder.add_node(3, 10.0);
        let n1 = builder.add_node(3, 10.0);
        let n2 = builder.add_node(3, 5.0);
        let n3 = builder.add_node(3, 10.0);

        builder.add_edge(n0, n1, 100, true, true);
        builder.add_edge(n0, n2, 100, true, true);
        builder.add_edge(n1, n3, 100, true, true);
        builder.add_edge(n2, n3, 100, true, true);

        let (graph, _reads) = builder.build();
        
        // Test bubble detection from n0
        let result = graph.double_bubble_remove_nodes(Direction::Outgoing, n0, 5000, usize::MAX);
        assert!(result.is_some());
        
        let nodes_to_remove = result.unwrap().remove_nodes;
        // Should remove the lower coverage path
        assert!(graph.nodes[&nodes_to_remove[0]].min_read_depth.unwrap() == 5.0);
        assert_eq!(nodes_to_remove.len(), 1);

        //Bidirected shenanigans
        // Create a simple bubble:
            // n0 -> n1 >-< n3
            //  ↘        ↗  
            //   -> n2 ->

        let mut builder = MockUnitigBuilder::new();
        
        let n0 = builder.add_node(3, 10.0);
        let n1 = builder.add_node(3, 10.0);
        let n2 = builder.add_node(3, 5.0);
        let n3 = builder.add_node(3, 10.0);

        builder.add_edge(n0, n1, 100, true, true);
        builder.add_edge(n0, n2, 100, true, true);
        builder.add_edge(n1, n3, 100, true, false);
        builder.add_edge(n2, n3, 100, true, false);

        let (graph, _reads) = builder.build();
        
        // Test bubble detection from n0
        let result = graph.double_bubble_remove_nodes(Direction::Outgoing, n0, 5000, usize::MAX);
        assert!(result.is_some());
        
        let nodes_to_remove = result.unwrap().remove_nodes;
        // Should remove the lower coverage path
        assert!(graph.nodes[&nodes_to_remove[0]].min_read_depth.unwrap() == 5.0);
        assert_eq!(nodes_to_remove.len(), 1);

    }

    #[test]
    fn test_bubble_detection_simple_cicular() {
        // Create a simple bubble:
        //--- n0 -> n1 -> n0 ----
        //    ↘        ↗  
        //      -> n2 ->
        let mut builder = MockUnitigBuilder::new();
        
        let n0 = builder.add_node(300, 10.0);
        let n1 = builder.add_node(3, 10.0);
        let n2 = builder.add_node(3, 5.0);

        builder.add_edge(n0, n1, 100, true, true);
        builder.add_edge(n0, n2, 100, true, true);
        builder.add_edge(n1, n0, 100, true, true);
        builder.add_edge(n2, n0, 100, true, true);

        let (graph, _reads) = builder.build();
        
        // Test bubble detection from n0
        let result = graph.double_bubble_remove_nodes(Direction::Outgoing, n0, 5000, usize::MAX);
        assert!(result.is_some());
        
        let nodes_to_remove = result.unwrap().remove_nodes;
        // Should remove the lower coverage path
        assert!(graph.nodes[&nodes_to_remove[0]].min_read_depth.unwrap() == 5.0);
        assert_eq!(nodes_to_remove.len(), 1);
    }

    #[test]
    fn test_bubble_detection_simple_extender() {
        // Create a simple bubble:
        //n4 -> n0 -> n1 -> n3 -> n5
            //  ↘        ↗  
            //   -> n2 ->
        let mut builder = MockUnitigBuilder::new();
        
        let n0 = builder.add_node(3, 10.0);
        let n1 = builder.add_node(3, 10.0);
        let n2 = builder.add_node(3, 5.0);
        let n3 = builder.add_node(3, 10.0);
        let n4 = builder.add_node(3, 10.0);
        let n5 = builder.add_node(3, 10.0);

        builder.add_edge(n0, n1, 100, true, true);
        builder.add_edge(n0, n2, 100, true, true);
        builder.add_edge(n1, n3, 100, true, true);
        builder.add_edge(n2, n3, 100, true, true);
        builder.add_edge(n4, n0, 100, true, true);
        builder.add_edge(n3, n5, 100, true, true);

        let (graph, _reads) = builder.build();
        
        // Test bubble detection from n0
        let result = graph.double_bubble_remove_nodes(Direction::Outgoing, n0, 5000, usize::MAX);
        assert!(result.is_some());
        
        let nodes_to_remove = result.unwrap().remove_nodes;
        // Should remove the lower coverage path
        assert!(graph.nodes[&nodes_to_remove[0]].min_read_depth.unwrap() == 5.0);
        assert_eq!(nodes_to_remove.len(), 1);
    }

    #[test]
    fn test_bubble_detection_open_ended() {
        // Create a simple bubble:
        // n0 -> n1 -> n3 
        //  ↘     ↘ ↗  
        //   -> n2 -> n4 
        let mut builder = MockUnitigBuilder::new();
        
        let n0 = builder.add_node(3, 10.0);
        let n1 = builder.add_node(3, 10.0);
        let n2 = builder.add_node(3, 5.0);
        let n3 = builder.add_node(3, 10.0);
        let n4 = builder.add_node(3, 5.0);

        builder.add_edge(n0, n1, 100, true, true);
        builder.add_edge(n0, n2, 100, true, true);
        builder.add_edge(n1, n3, 100, true, true);
        builder.add_edge(n2, n3, 100, true, true);
        builder.add_edge(n1, n4, 100, true, true);
        builder.add_edge(n2, n4, 100, true, true);

        let (graph, _reads) = builder.build();
        
        //Should not be a bubble
        let result1 = graph.double_bubble_remove_nodes(Direction::Outgoing, n0, 5000, usize::MAX);
        let result2 = graph.double_bubble_remove_nodes(Direction::Incoming, n3, 5000, usize::MAX);
        let result3 = graph.double_bubble_remove_nodes(Direction::Incoming, n4, 5000, usize::MAX);
        dbg!(&result1, &result2, &result3);

        assert!(result1.is_none());
        assert!(result2.is_none());
        assert!(result3.is_none());
    }

    #[test]
    fn test_bubble_detection_redirect() {
        // Create a simple bubble:
        // n0 -> n1 -> n3 <- -> n4
        //  ↘      ↗      | 
        //   -> n2  ------|
        let mut builder = MockUnitigBuilder::new();
        
        let n0 = builder.add_node(3, 10.0);
        let n1 = builder.add_node(3, 10.0);
        let n2 = builder.add_node(3, 5.0);
        let n3 = builder.add_node(3, 10.0);
        let n4 = builder.add_node(3, 5.0);

        builder.add_edge(n0, n1, 100, true, true);
        builder.add_edge(n0, n2, 100, true, true);
        builder.add_edge(n1, n3, 100, true, true);
        builder.add_edge(n2, n3, 100, true, true);
        builder.add_edge(n2, n3, 100, true, false);
        builder.add_edge(n3, n4, 100, true, true);

        let (graph, _reads) = builder.build();
        
        //Should not be a bubble
        let result1 = graph.double_bubble_remove_nodes(Direction::Outgoing, n0, 5000, usize::MAX);
        let result2 = graph.double_bubble_remove_nodes(Direction::Incoming, n3, 5000, usize::MAX);
        let result3 = graph.double_bubble_remove_nodes(Direction::Incoming, n4, 5000, usize::MAX);
        dbg!(&result1, &result2, &result3);

        assert!(result1.is_none());
        assert!(result2.is_none());
        assert!(result3.is_none());
    }

    #[test]
    fn test_bubble_detection_self_loop() {
        // Create a simple bubble:
        // n0 -> n1 -> n3
        //  ↘        ↗  
        //   -> n2 ->
        //     ↘^ 
        let mut builder = MockUnitigBuilder::new();
        
        let n0 = builder.add_node(3, 10.0);
        let n1 = builder.add_node(3, 10.0);
        let n2 = builder.add_node(3, 5.0);
        let n3 = builder.add_node(3, 10.0);

        builder.add_edge(n0, n1, 100, true, true);
        builder.add_edge(n0, n2, 100, true, true);
        builder.add_edge(n1, n3, 100, true, true);
        builder.add_edge(n2, n3, 100, true, true);
        builder.add_edge(n2, n2, 100, true, true);

        let (graph, _reads) = builder.build();
        
        // Test bubble detection from n0
        let result = graph.double_bubble_remove_nodes(Direction::Outgoing, n0, 5000, usize::MAX);
        assert!(result.is_none());
    }

    #[test]
    fn test_bubble_detection_complex() {
        // Create a simple bubble:
        // n0 -> n1 -> n3 -> 
        //  ↘        ↗      n5
        //   -> n2 -> n4 ↗ 
        let mut builder = MockUnitigBuilder::new();
        
        let n0 = builder.add_node(3, 10.0);
        let n1 = builder.add_node(3, 10.0);
        let n2 = builder.add_node(3, 5.0);
        let n3 = builder.add_node(3, 10.0);
        let n4 = builder.add_node(3, 20.0);
        let n5 = builder.add_node(3, 10.0);

        builder.add_edge(n0, n1, 100, true, true);
        builder.add_edge(n0, n2, 100, true, true);
        builder.add_edge(n1, n3, 100, true, true);
        builder.add_edge(n2, n3, 100, true, true);
        builder.add_edge(n2, n4, 100, true, true);
        builder.add_edge(n4, n5, 100, true, true);
        builder.add_edge(n3, n5, 100, true, true);

        let (graph, _reads) = builder.build();
        
        // Test bubble detection from n0
        let result = graph.double_bubble_remove_nodes(Direction::Incoming, n5, 50000, usize::MAX);
        assert!(result.is_some());
        let nodes_to_remove = result.unwrap().remove_nodes;
        assert!(nodes_to_remove.len() == 2);
        assert!(nodes_to_remove.contains(&n3));
        assert!(nodes_to_remove.contains(&n1));
    }

    // Test safe path detection
    #[test]
    fn test_safe_path_detection() {
        // Create a graph with a safe path:
        // ----> n1 -> n2 
        //       ↘     
        // n4 ->  n5 ----- 
        let mut builder = MockUnitigBuilder::new();
        
        let n1 = builder.add_node(100, 10.0);
        let n2 = builder.add_node(100, 10.0);
        let n4 = builder.add_node(100, 10.0);
        let n5 = builder.add_node(100, 10.0);

        builder.add_edge(n1, n2, 100, true, true);
        builder.add_edge(n1, n5, 100, true, true);
        builder.add_edge(n4, n5, 100, true, true);
        
        let (graph, _reads) = builder.build();

        // Test if path through n1->n5 is safe
        let edge = graph.edges[1].as_ref().unwrap();
        let result = graph.safe_given_forward_back(
            &graph.nodes[&n1],
            edge,
            2000, // max_forward
            5,    // max_reads_forward
            1000, // safe_length_back
            None,
            &FxHashSet::default()
        );
        
        assert!(result);

        let edge = graph.edges[0].as_ref().unwrap();
        let result = graph.safe_given_forward_back(
            &graph.nodes[&n1],
            edge,
            2000, // max_forward
            5,    // max_reads_forward
            1000, // safe_length_back
            None,
            &FxHashSet::default()
        );
        
        assert!(!result);

    }

    #[test]
    fn safely_cut_edge_test_basic_x(){
        let mut builder = MockUnitigBuilder::new();

        // Make sure not to cut the safe edge after multiple iterations
        //   1k
        // 1 ---> 2
        //    X    (3->2 : 10k), (3->2 : 1k)
        // 3 ---> 4
        //   10k
        
        let n1 = builder.add_node(100, 10.0);
        let n2 = builder.add_node(100, 10.0);
        let n3 = builder.add_node(100, 10.0);
        let n4 = builder.add_node(100, 10.0);

        builder.add_edge(n1, n2, 1000, true, true);
        builder.add_edge(n1, n4, 10000, true, true);
        builder.add_edge(n3, n4, 10000, true, true);
        builder.add_edge(n3, n2, 1000, true, true);
        
        let (graph, _reads) = builder.build();

        let mut removed_edges = FxHashSet::default();
        //let mut unitig_edge_file = BufWriter::new(std::fs::File::create("unitig_edge_file").unwrap());
        //create empty mock file
        let mut unitig_edge_file = BufWriter::new(std::io::sink());

        //Cut 0 (ol 1000) but don't cut 1 (ol 1000)
        let edge_order = vec![0, 3, 2, 1];

        for edge_id in edge_order{
            graph.safely_cut_edge(
                edge_id,
                &mut removed_edges,
                0.5,
                None,
                None,
                2000,
                5,
                1000,
                &mut unitig_edge_file,
            );
        }

        assert!(removed_edges.contains(&0));
        assert!(!removed_edges.contains(&3));
    }

    #[test]
    fn safely_cut_edge_test_basic_x_cov(){
        let mut builder = MockUnitigBuilder::new();

        // Make sure not to cut the safe edge after multiple iterations
        //          1k
        // 1 (100x) ---> 2 (10x)
        //           X 
        // 3 (10x) ---> 4 (100x)
        //          10k
        
        let n1 = builder.add_node(100, 100.0);
        let n2 = builder.add_node(100, 10.0);
        let n3 = builder.add_node(100, 10.0);
        let n4 = builder.add_node(100, 100.0);

        builder.add_edge(n1, n2, 1000, true, true);
        builder.add_edge(n1, n4, 10000, true, true);
        builder.add_edge(n3, n4, 10000, true, true);
        builder.add_edge(n3, n2, 1000, true, true);
        
        let (graph, _reads) = builder.build();

        let mut removed_edges = FxHashSet::default();
        //let mut unitig_edge_file = BufWriter::new(std::fs::File::create("unitig_edge_file").unwrap());
        //create empty mock file
        let mut unitig_edge_file = BufWriter::new(std::io::sink());

        //Cut 0 (ol 1000) but don't cut 1 (ol 1000)
        let edge_order = vec![0, 3, 2, 1];

        for edge_id in edge_order{
            graph.safely_cut_edge(
                edge_id,
                &mut removed_edges,
                0.5,
                Some(3.),
                None,
                2000,
                5,
                1000,
                &mut unitig_edge_file,
            );
        }

        assert!(removed_edges.contains(&0));
        assert!(!removed_edges.contains(&3));

        removed_edges.clear();

        // Edge 1 is cut due to cov, edge 0 is cut due to ol ratio and cov and is compatible
        let edge_order = vec![1,2,3,0];

        for edge_id in edge_order{
            graph.safely_cut_edge(
                edge_id,
                &mut removed_edges,
                0.5,
                Some(3.),
                None,
                2000,
                5,
                1000,
                &mut unitig_edge_file,
            );
        }

        dbg!(&removed_edges);

        assert!(removed_edges.contains(&0));
        assert!(removed_edges.contains(&2));
        assert!(!removed_edges.contains(&1));
        assert!(!removed_edges.contains(&3));

    }
}

