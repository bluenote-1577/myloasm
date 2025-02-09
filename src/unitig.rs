use crate::cli::Cli;
use crate::constants::ID_THRESHOLD_ITERS;
use crate::constants::PSEUDOCOUNT;
use crate::constants::QUANTILE_UNITIG_WEIGHT;
use crate::graph::*;
use crate::twin_graph::*;
use crate::types::*;
use crate::unitig_utils::*;
use bio_seq::prelude::*;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use rand::prelude::*;
use rayon::prelude::*;
use rust_lapper::Lapper;
use std::collections::VecDeque;
use std::io::BufWriter;
use std::io::Write;
use std::panic;
use std::sync::Mutex;

#[derive(Debug, Clone, PartialEq, Default)]
pub struct UnitigNode {
    pub read_indices_ori: Vec<(NodeIndex, bool)>,
    pub internal_overlaps: Vec<ReadOverlapEdgeTwin>,
    pub read_names: Vec<String>,
    in_edges: Vec<EdgeIndex>,
    out_edges: Vec<EdgeIndex>,
    pub node_id: NodeIndex,
    pub node_hash_id: NodeIndex,
    pub min_read_depth_multi: Option<MultiCov>, //Depth as measured by median of minimum depth of reads within unitig
    pub median_read_depth: Option<f64>, //Depth as measured by median of median depth of reads within unitig
    pub unique_length: Option<usize>, // Length of the unitig that is not covered by any other unitig's overlap; can be 0
    pub read_min_depths_multi: Vec<(MultiCov, usize)>,
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
    fn edge_id_est(&self, c: usize) -> f64 {
        self.overlap.edge_id_est(c)
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
    fn edge_id_est(&self, c: usize) -> f64 {
        self.overlap.edge_id_est(c)
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

                min_read_depths.extend(unitig.read_min_depths_multi.iter().map(|x| *x));
                median_read_depths.extend(unitig.read_median_depths.iter().map(|x| *x));
            }
            let median_min_depth =
                median_weight_multi(&mut min_read_depths, QUANTILE_UNITIG_WEIGHT);
            let median_median_depth =
                median_weight(&mut median_read_depths, QUANTILE_UNITIG_WEIGHT);
            let new_unitig = UnitigNode {
                internal_overlaps: new_internal_overlaps,
                read_names: new_read_names,
                in_edges: Vec::new(),
                out_edges: Vec::new(),
                node_id: new_read_indices_ori[0].0,
                read_indices_ori: new_read_indices_ori,
                node_hash_id: new_unitig_graph.nodes.len(),
                min_read_depth_multi: median_min_depth,
                median_read_depth: median_median_depth,
                base_info: BaseInfo::default(),
                mapping_info: MappingInfo::default(),
                read_min_depths_multi: min_read_depths,
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
    pub fn from_overlaps(reads: &[TwinRead], overlaps: Vec<OverlapConfig>, args: &Cli) -> (Self, OverlapAdjMap) {

        let mut adj_map = FxHashMap::default();
        for overlap in overlaps.iter() {
            if overlap.overlap1_len > 1000 {
                adj_map.entry(overlap.read_i).or_insert_with(Vec::new).push(overlap.read_j);
                adj_map.entry(overlap.read_j).or_insert_with(Vec::new).push(overlap.read_i);
            }
        }
        let overlap_adj_map = OverlapAdjMap {
            adj_map,
        };

        let overlap_graph = read_graph_from_overlaps_twin(overlaps, reads, args);
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
                    let depth_multi = read.min_depth_multi.as_ref().unwrap().clone();
                    let length = read.base_length;
                    (depth_multi, length)
                })
                .collect::<Vec<(MultiCov, usize)>>();
            let median_min_depth =
                median_weight_multi(&mut median_min_depth_vec, QUANTILE_UNITIG_WEIGHT);

            let mut median_median_depth_vec = read_indices_ori
                .iter()
                .map(|(idx, _)| {
                    let read = &reads[*idx];
                    let depth = read.median_depth.unwrap();
                    let length = read.base_length;
                    (depth as f64, length)
                })
                .collect::<Vec<(f64, usize)>>();
            let median_median_depth =
                median_weight(&mut median_median_depth_vec, QUANTILE_UNITIG_WEIGHT);

            let unitig = UnitigNode {
                read_indices_ori,
                read_names,
                internal_overlaps: overlaps,
                in_edges: Vec::new(),
                out_edges: Vec::new(),
                node_id: node_id,
                node_hash_id: unitig_graph.nodes.len(),
                min_read_depth_multi: median_min_depth,
                median_read_depth: median_median_depth,
                base_info: BaseInfo::default(),
                mapping_info: MappingInfo::default(),
                read_median_depths: median_median_depth_vec,
                read_min_depths_multi: median_min_depth_vec,
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

        (unitig_graph, overlap_adj_map)
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
            name += &format!(" len:{}_depth:{:?}", unitig.cut_length(), unitig.min_read_depth_multi.unwrap());
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
            let min_read_depth_multi = unitig
                .min_read_depth_multi
                .unwrap_or([-1.; ID_THRESHOLD_ITERS]);

            gfa.push_str(&format!(
                "S\tu{}\t{}\tLN:i:{}\tDP:f:{:.1}\tDP2:f:{:.1}\tDP3:f:{:.1}\tMED_PRIM_DP:f:{:.1}\tMAP_DP:f:{:.1}\tMED_MAP_DP:f:{:.1}\n",
                unitig.read_indices_ori[0].0,
                //String::from_utf8(unitig.raw_consensus()).unwrap(),
                base_seq,
                unitig.cut_length(),
                min_read_depth_multi[0],
                min_read_depth_multi[1],
                min_read_depth_multi[2],
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
                        "a\tu{},{}-{}\t{}\t{}\t{}\t{}\t{}\tDP1:{},DP2:{},DP3:{}\n",
                        unitig.read_indices_ori[0].0,
                        range.0,
                        range.1,
                        read_idx,
                        curr_pos,
                        read.id,
                        ori_string,
                        length,
                        read.min_depth_multi.as_ref().unwrap()[0],
                        read.min_depth_multi.as_ref().unwrap()[1],
                        read.min_depth_multi.as_ref().unwrap()[2],
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

        log::trace!("Removing {} tips", unitigs_to_remove.len());
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
                    if let Some(bubble_result) = self.double_bubble_remove_nodes(
                        *direction,
                        n_id,
                        max_length,
                        max_number_nodes,
                    ) {
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
            "BUBBLE: Removed {} bubbles at max length {}",
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
        max_number_nodes: usize,
    ) -> Option<BubblePopResult> {
        let opt = self.get_bubble_remove_nodes(direction, n_id, max_length, max_number_nodes);
        if let Some(bubble_result) = opt {
            if let Some(bubble_result_back) = self.get_bubble_remove_nodes(
                bubble_result.end_direction,
                bubble_result.sink_hash_id,
                max_length,
                max_number_nodes,
            ) {
                if bubble_result.remove_edges == bubble_result_back.remove_edges {
                    return Some(bubble_result);
                }
            }
        }

        return None;
    }

    fn get_bubble_remove_nodes(
        &self,
        right_direction: Direction,
        n_id: NodeIndex,
        max_length: usize,
        max_number_nodes: usize,
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
            if distances.len() > max_number_nodes {
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
                        prev_score
                            + other_node.min_read_depth_multi.unwrap().iter().sum::<f64>()
                                * num_reads,
                    );
                } else {
                    if prev_dist + other_node.unique_length.unwrap() < distances[&other_node_ind] {
                        distances.insert(
                            other_node_ind,
                            prev_dist + other_node.unique_length.unwrap(),
                        );
                    }
                    if prev_score
                        + other_node.min_read_depth_multi.unwrap().iter().sum::<f64>() * num_reads
                        > depth_length[&other_node_ind]
                    {
                        depth_length.insert(
                            other_node_ind,
                            prev_score
                                + other_node.min_read_depth_multi.unwrap().iter().sum::<f64>()
                                    * num_reads,
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
                    log::trace!("Bubble found between {} and {}", start_node_id, end_node_id);
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

                    log::trace!("Best path {:?}", &debug_path);
                    if !debug_path.contains(&start_node_id) || !debug_path.contains(&end_node_id) {
                        log::trace!(
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
                        log::trace!(
                            "Path between {} and {} does not start and end at the correct nodes",
                            start_node_id,
                            end_node_id
                        );
                        return None;
                    }
                    return Some(BubblePopResult {
                        original_direction: right_direction,
                        end_direction: stack[0].1.reverse(),
                        source_hash_id: start_node_hash_id,
                        sink_hash_id: end_node_hash_id,
                        remove_nodes: remove_vertices,
                        remove_edges: remove_edges,
                    });
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
        log::debug!("Z-EDGE: Cut {} z-edges", edges_to_remove.len());
        self.remove_edges(edges_to_remove);
        self.re_unitig();
    }

    fn get_z_edges(&self, unitig: &UnitigNode, direction: Direction) -> Vec<EdgeIndex> {
        let edges;

        let mut non_z_edge_lengths = vec![];
        let mut z_edge_len = None;
        let mut z_edge_id = None;
        let mut unique_edge_ids = FxHashSet::default();

        if direction == Direction::Outgoing {
            edges = unitig.out_edges()
        } else {
            edges = unitig.in_edges()
        }
        unique_edge_ids.extend(edges.iter());

        if edges.len() != 2 {
            return vec![];
        }

        let mut potential_z_edge = false;
        let mut side_edge_1 = false;
        for edge_id in edges.iter() {
            let edge = self.edges[*edge_id].as_ref().unwrap();
            let other_node_id = edge.other_node(unitig.node_hash_id);
            let other_unitig = &self.nodes[&other_node_id];
            let dir = edge.node_edge_direction(&other_node_id);
            let other_edges = other_unitig.edges_direction(&dir);
            if other_edges.len() == 2 {
                z_edge_len = Some(edge.overlap.overlap_len_bases);
                z_edge_id = Some(*edge_id);
                let side_edge_2_id = other_edges
                    .iter()
                    .filter(|x| **x != *edge_id)
                    .collect::<Vec<_>>();
                if side_edge_2_id.len() == 1 {
                    let side_edge_2 = self.edges[*side_edge_2_id[0]].as_ref().unwrap();
                    unique_edge_ids.insert(*side_edge_2_id[0]);
                    non_z_edge_lengths.push(side_edge_2.overlap.overlap_len_bases);
                    potential_z_edge = true;
                } else {
                    return vec![];
                }
            } else if other_edges.len() == 1 {
                side_edge_1 = true;
                non_z_edge_lengths.push(edge.overlap.overlap_len_bases);
            } else {
                return vec![];
            }
        }

        if potential_z_edge && side_edge_1 && unique_edge_ids.len() == 3 {
            if non_z_edge_lengths[0] > z_edge_len.unwrap()
                && non_z_edge_lengths[1] > z_edge_len.unwrap()
            {
                return vec![z_edge_id.unwrap()];
            }
        }

        return vec![];
    }

    fn _get_z_edges_old(&self, unitig: &UnitigNode, direction: Direction) -> Vec<EdgeIndex> {
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
        starting_unitig: &UnitigNode,
        starting_edge: &UnitigEdge,
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
        if starting_edge.from_unitig == starting_edge.to_unitig {
            if self.nodes[&starting_edge.from_unitig].in_edges().len() > 1
                && self.nodes[&starting_edge.from_unitig].out_edges().len() > 1
            {
                return true;
            }
        }

        let max_length_forward_search = max_forward;
        let length_back_safe = safe_length_back;
        let mut safe = false;

        let mut seen_nodes = FxHashSet::default();
        let mut to_search = VecDeque::new();

        let other_node = starting_edge.other_node(starting_unitig.node_hash_id);
        let other_node_back_dir = starting_edge.node_edge_direction(&other_node);

        //Node and back-direction
        to_search.push_back((other_node, other_node_back_dir));
        seen_nodes.insert(starting_unitig.node_hash_id);
        let mut travelled = 0;
        let mut travelled_reads = 0;

        let mut forbidden_nodes = FxHashSet::default();
        forbidden_nodes.insert(starting_unitig.node_hash_id);

        //Search forward until a node has sufficient "backing"
        while travelled < max_length_forward_search && travelled_reads < max_reads_forward {
            //TODO break...? small tips
            if to_search.len() == 0 {
                safe = true;
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
            let back_found = self.search_dir_until_safe(
                unitig,
                direction,
                length_back_safe,
                safety_cov_edge_ratio,
                &forbidden_nodes,
                removed_edges,
                Some((starting_unitig, starting_edge)),
            );
            if back_found {
                safe = true;
                break;
            }

            forbidden_nodes.insert(node);
            let forward_edges = self.get_safe_edges_from_cov_threshold(
                unitig.edges_direction(&direction.reverse()),
                safety_cov_edge_ratio,
            );
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
        args: &Cli,
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

        for (i, edge) in self.edges.iter().enumerate() {
            if let Some(edge) = edge {
                let uni1 = &self.nodes[&edge.from_unitig];
                let uni2 = &self.nodes[&edge.to_unitig];
                let overlap_len = edge.overlap.overlap_len_bases;

                let cov_ratio = 1.
                    / pseudocount_cov_multi(
                        uni1.min_read_depth_multi.unwrap(),
                        uni2.min_read_depth_multi.unwrap(),
                    );
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
                None,
                true,
                max_forward,
                max_reads_forward,
                safe_length_back,
                &mut unitig_edge_file,
                args.c
            );
        }
        log::debug!(
            "BRIDGED: Cutting {} edges that have low relative confidence ({})",
            removed_edges.len(),
            ol_thresh
        );
        self.remove_edges(removed_edges);
        self.re_unitig();
    }

    pub fn print_statistics(&self, args: &Cli) {
        //Print n50, name of largest contig, largest contig size, number of contigs

        let mut contig_sizes = self
            .nodes
            .iter()
            .filter(|(_, node)| {
                if node.read_indices_ori.len() <= args.min_reads_contig
                    && node.in_edges().len() + node.out_edges().len() == 0
                {
                    return false;
                } else {
                    return true;
                }
            })
            .map(|(_, node)| node.cut_length())
            .collect::<Vec<_>>();
        contig_sizes.sort();
        let n50 = contig_sizes.iter().sum::<usize>() / 2;
        let contig_sizes_iterev = contig_sizes.iter().rev();
        let mut curr_sum = 0;
        let mut n50_size = 0;
        for size in contig_sizes_iterev {
            curr_sum += size;
            if curr_sum >= n50 {
                n50_size = *size;
                break;
            }
        }

        let num_circ_contigs_geq_1m = self
            .nodes
            .iter()
            .filter(|(_, node)| {
                node.is_circular()
                    && node.cut_length() >= 1_000_000
            })
            .count();

        let num_circ_contigs_geq_100k = self
            .nodes
            .iter()
            .filter(|(_, node)| {
                node.is_circular()
                    && node.cut_length() >= 100_000
            })
            .count();

        let num_1m_contigs = self
            .nodes
            .iter()
            .filter(|(_, node)| {
                node.cut_length() >= 1_000_000
            })
            .count();

        let largest_contig = self.nodes.iter().max_by_key(|(_, node)| node.cut_length());
        if let Some(largest_contig) = largest_contig {
            let largest_contig_size = largest_contig.1.cut_length();
            let num_contigs = contig_sizes.len();
            log::info!("-------------- Assembly statistics --------------");
            log::info!("N50: {}", n50_size);
            log::info!("Largest contig has size: {}", largest_contig_size);
            log::info!("Number of contigs: {}", num_contigs);
            log::info!("Number of circular contigs >= 1M: {}", num_circ_contigs_geq_1m);
            log::info!("Number of contigs (circular or non-circular)>= 1M: {}", num_1m_contigs);
            log::info!("Number of circular contigs >= 100k: {}", num_circ_contigs_geq_100k);
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
        removed_edges: &FxHashSet<EdgeIndex>,
        starting_node_edge: Option<(&UnitigNode, &UnitigEdge)>,
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
            let edges = self.get_safe_edges_from_cov_threshold(
                self.nodes[&node].edges_direction(&dir),
                safety_cov_edge_ratio,
            );
            let curr_min_dist = nodes_and_distances[&node];
            for edge_id in edges {
                let edge = self.edges[edge_id].as_ref().unwrap();
                let other_node = edge.other_node(node);

                // Circular condition -- if the starting node is the same as the searched node
                // but the directions are different then we check if the length of the circular contig is
                // sufficient.
                if let Some((start_node, start_edge)) = starting_node_edge {
                    if start_node.node_hash_id == other_node {
                        let direction_start =
                            start_edge.node_edge_direction(&start_node.node_hash_id);
                        let direction_search = edge.node_edge_direction(&other_node);
                        if start_node.node_hash_id == other_node
                            && direction_start != direction_search
                        {
                            assert!(start_node.base_info_present());
                            if start_node.base_info.length > length {
                                return true;
                            }
                        }
                    }
                }

                if forbidden_nodes.contains(&other_node) || removed_edges.contains(&edge_id) {
                    continue;
                }

                let other_unitig = &self.nodes[&other_node];

                if nodes_and_distances.contains_key(&other_node) {
                    *nodes_and_distances.get_mut(&other_node).unwrap() = nodes_and_distances
                        [&other_node]
                        .min(curr_min_dist + other_unitig.base_info.length);
                } else {
                    nodes_and_distances
                        .insert(other_node, curr_min_dist + other_unitig.base_info.length);
                    let other_dir = edge.node_edge_direction(&other_node).reverse();
                    to_search.push_back((other_node, other_dir));
                }

                if nodes_and_distances[&other_node] > length {
                    achieved_length = true;
                    break;
                }

                iteration += 1;
                if iteration > 25 {
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

    pub fn cut_coverage(&mut self, cov: f64) {
        let mut nodes_to_remove = FxHashSet::default();
        for node in self.nodes.values() {
            if node.min_read_depth_multi.unwrap().iter().sum::<f64>() / (ID_THRESHOLD_ITERS as f64)
                < cov
            {
                nodes_to_remove.insert(node.node_hash_id);
            }
        }
        self.remove_nodes(&nodes_to_remove.into_iter().collect::<Vec<_>>(), false);
        self.re_unitig();
    }

    fn safely_cut_edge<T>(
        &self,
        edge_id: EdgeIndex,
        removed_edges: &mut FxHashSet<EdgeIndex>,
        ol_thresh: f64,
        unitig_cov_ratio_cut: Option<f64>,
        safety_cov_edge_ratio: Option<f64>,
        strain_repeat_map: Option<&FxHashMap<NodeIndex, FxHashSet<NodeIndex>>>,
        snpmer_id_est_safety: bool,
        max_forward: usize,
        max_reads_forward: usize,
        safe_length_back: usize,
        unitig_edge_file: &mut T,
        c: usize,
    ) where
        T: Write,
    {
        let edge = &self.edges[edge_id].as_ref().unwrap();
        let unitig_terminals = [edge.from_unitig, edge.to_unitig];

        writeln!(
            unitig_edge_file,
            "u{}-u{}, ol:{}",
            self.nodes[&edge.from_unitig].node_id,
            self.nodes[&edge.to_unitig].node_id,
            edge.overlap.overlap_len_bases
        )
        .unwrap();

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
                let fsv = edge.edge_id_est(c);
                overlaps.push((ol, fsv));
            }

            let (max_ol, max_ols_fsv) = *overlaps.iter().max_by(|a, b| a.0.cmp(&b.0)).unwrap();
            let ol_score = edge.overlap.overlap_len_bases as f64 / max_ol as f64; // higher better

            // We don't cut if the fsv of the cut edge is higher than the longest edge
            let fsv_edge = edge.edge_id_est(c);
            if fsv_edge > max_ols_fsv && snpmer_id_est_safety{
                continue;
            }

            let uni1 = &self.nodes[&edge.from_unitig];
            let uni2 = &self.nodes[&edge.to_unitig];

            let cov_pseudo_val = pseudocount_cov_multi(
                uni1.min_read_depth_multi.unwrap(),
                uni2.min_read_depth_multi.unwrap(),
            );

            let mut potential_cut = false;
            if (ol_score < ol_thresh)
                || (unitig_cov_ratio_cut.is_some()
                    && cov_pseudo_val > unitig_cov_ratio_cut.unwrap())
            {
                potential_cut = true;
            }

            if !potential_cut {
                continue;
            }

            // Check all edges besides the current edge. 
            // If the current edge's destination is a strain repeat
            // of any other edge, this edge cut is safe. 
            let mut strain_repeat_safe = false;
            if let Some(strain_repeat_map) = strain_repeat_map {
                let mut possible_repeats: FxHashSet<NodeIndex> = FxHashSet::default();
                for &edge_id_iter in non_cut_edges.iter(){
                    if edge_id_iter == edge_id{
                        continue;
                    } 
                    let other_node_id = self.edges[edge_id_iter].as_ref().unwrap().other_node(unitig.node_hash_id);
                    if let Some(repeats) = strain_repeat_map.get(&other_node_id){
                        possible_repeats.extend(&mut repeats.iter());
                    }
                }
                let current_other_node_id = edge.other_node(unitig.node_hash_id);
                if possible_repeats.contains(&current_other_node_id){
                    strain_repeat_safe = true;
                }
            }


            // Make sure another path, besides this edge, is safe.
            let mut forward_search_forbidden_nodes = FxHashSet::default();
            forward_search_forbidden_nodes.insert(unitig.node_hash_id);
            forward_search_forbidden_nodes.insert(edge.other_node(unitig.node_hash_id));
            let mut safe_if_cut = self.search_dir_until_safe(
                unitig,
                direction,
                safe_length_back,
                safety_cov_edge_ratio,
                &forward_search_forbidden_nodes,
                removed_edges,
                Some((unitig, edge)),
            );

            //Check circularity condition, as safe_if_cut is true if cutting -> circular contig
            if !safe_if_cut {
                for edge_id_check in non_cut_edges.iter() {
                    if edge_id_check == &edge_id {
                        continue;
                    }
                    let edge = self.edges[*edge_id_check].as_ref().unwrap();
                    if edge.from_unitig == edge.to_unitig {
                        safe_if_cut = true;
                        break;
                    }
                }
            }

            drop(forward_search_forbidden_nodes);

            if !safe_if_cut && !strain_repeat_safe {
                writeln!(
                    unitig_edge_file,
                    "u{}-u{}, not safe if cut",
                    self.nodes[&edge.from_unitig].node_id, self.nodes[&edge.to_unitig].node_id,
                )
                .unwrap();
                continue;
            }

            let safe = self.safe_given_forward_back(
                unitig,
                edge,
                max_forward,
                max_reads_forward,
                safe_length_back,
                safety_cov_edge_ratio,
                removed_edges,
            );

            if (safe || strain_repeat_safe)
                && ((ol_score < ol_thresh)
                    || (unitig_cov_ratio_cut.is_some()
                        && cov_pseudo_val > unitig_cov_ratio_cut.unwrap()))
            {
                removed_edges.insert(edge_id);
            } 
            writeln!(
                unitig_edge_file,
                "u{}-u{}, cut:{} safe:{} snp_share:{}, snp_diff:{}, ol_length:{}, ol_score:{}, specific_score:{}, sr_safe:{}",
                uni1.node_id,
                uni2.node_id,
                true,
                safe,
                edge.overlap.shared_snpmers,
                edge.overlap.diff_snpmers,
                edge.overlap.overlap_len_bases,
                ol_score,
                cov_pseudo_val,
                strain_repeat_safe,
            )
            .unwrap();
        }
    }

    fn traverse_walk_probabilistic(
        &self,
        direction: Direction,
        unitig_id: Option<NodeIndex>,
        temperature: f64,
        _max_steps: usize,
        coverages: &mut Vec<(MultiCov, usize)>,
        rng: &mut StdRng,
    ) -> Option<EdgeIndex> {
        let unitig_id = match unitig_id {
            Some(unitig_id) => unitig_id,
            None => return None,
        };

        let unitig = &self.nodes[&unitig_id];
        let edges = unitig.edges_direction(&direction);
        if edges.is_empty() {
            return None;
        }
        let mut edge_probs = vec![];
        let mut total_prob = 0.0;
        let max_ol = edges
            .iter()
            .map(|x| self.edges[*x].as_ref().unwrap().overlap.overlap_len_bases)
            .max()
            .unwrap();
        for edge_id in edges {
            let edge = self.edges[*edge_id].as_ref().unwrap();
            let other_node = edge.other_node(unitig.node_hash_id);
            let other_unitig = &self.nodes[&other_node];

            let other_covs = &other_unitig.read_min_depths_multi;
            let cov_different = quantile_dist(coverages, other_covs);
            let ol = edge.overlap.overlap_len_bases;
            let ol_ratio = ol as f64 / max_ol as f64;

            let cov1;
            let cov2;

            if edge.f1 {
                cov1 = self.nodes[&edge.from_unitig]
                    .read_min_depths_multi
                    .last()
                    .unwrap()
                    .0[ID_THRESHOLD_ITERS - 1];
            } else {
                cov1 = self.nodes[&edge.from_unitig]
                    .read_min_depths_multi
                    .first()
                    .unwrap()
                    .0[ID_THRESHOLD_ITERS - 1];
            }

            if edge.f2 {
                cov2 = self.nodes[&edge.to_unitig]
                    .read_min_depths_multi
                    .first()
                    .unwrap()
                    .0[ID_THRESHOLD_ITERS - 1];
            } else {
                cov2 = self.nodes[&edge.to_unitig]
                    .read_min_depths_multi
                    .last()
                    .unwrap()
                    .0[ID_THRESHOLD_ITERS - 1];
            }

            let cov_ratio_weight = 1.0 - PSEUDOCOUNT / ((cov1.min(cov2)).max(1.0) + PSEUDOCOUNT);

            if let Some(cov_different) = cov_different {
                let edge_prob = 2.713_f64.powf(
                    ((cov_ratio_weight * (ol_ratio - 1.) - cov_different) / temperature).min(100.0),
                );
                total_prob += edge_prob;
                edge_probs.push((edge_id, edge_prob));
            }
        }

        //Use RNG_SEED
        let random_val: f64 = rng.gen(); // generates a float between 0 and 1
        let mut current_prob = 0.0;
        for (edge_id, edge_prob) in edge_probs {
            current_prob += edge_prob;
            if current_prob >= random_val * total_prob {
                let other_unitig_id = self.edges[*edge_id]
                    .as_ref()
                    .unwrap()
                    .other_node(unitig.node_hash_id);
                let other_unitig = &self.nodes[&other_unitig_id];
                coverages.extend(other_unitig.read_min_depths_multi.clone());
                return Some(*edge_id);
            }
        }

        panic!("Shouldn't get here");
    }

    pub fn get_strain_repeats(&self, overlap_adj_map: &OverlapAdjMap, _args: &Cli) -> FxHashMap<NodeIndex, FxHashSet<NodeIndex>> {
        let min_contig_length = 100_000;
        let mut strain_repeats = FxHashMap::default();
        let mut read_to_contig_map = FxHashMap::default();
        for (contig_id, contig) in self.nodes.iter() {
            if contig.cut_length() < min_contig_length {
                continue;
            }
            for (read_id,_) in contig.read_indices_ori.iter() {
                read_to_contig_map.insert(*read_id, *contig_id);
            }
        }
        for (contig_id, contig) in self.nodes.iter() {
            //Only large enough contigs
            if contig.cut_length() < min_contig_length {
                continue;
            }

            // Get contigs that share reads with this contig
            let mut strain_repeats_contig = FxHashMap::default();
            for (read_id,_) in contig.read_indices_ori.iter() {
                if let Some(similar_reads) = overlap_adj_map.adj_map.get(read_id) {
                    for similar_read in similar_reads.iter() {
                        if let Some(contig_id_similar) = read_to_contig_map.get(similar_read) {
                            if contig_id != contig_id_similar {
                                let count = strain_repeats_contig.entry(*contig_id_similar).or_insert(FxHashSet::default());
                                count.insert(similar_read);
                            }
                        }
                    }
                }
            }

            //Try half contain
            for (contig_id_repeats, count) in strain_repeats_contig.iter() {
                let contig = self.nodes.get(contig_id_repeats).unwrap();
                if contig.cut_length() < min_contig_length {
                    continue;
                }
                if count.len() > contig.read_indices_ori.len() / 2 {
                    let contig_set = strain_repeats.entry(*contig_id).or_insert(FxHashSet::default());
                    contig_set.insert((*contig_id_repeats, count.len()));
                }
            }
        }

        //Print for debugging
        for (contig_id, contig_set) in strain_repeats.iter() {
            let contig = self.nodes.get(contig_id).unwrap();
            let mut contig_str = String::new();
            for (contig_id_2, cont_reads) in contig_set.iter() {
                let contig = self.nodes.get(contig_id_2).unwrap();
                contig_str.push_str(&format!("{}, contained: {}, num reads: {} ", contig.node_id, cont_reads, contig.read_indices_ori.len()));
            }
            log::debug!("Contig {} has strain repeats: {}", contig.node_id, contig_str);
        }

        return strain_repeats.into_iter().map(|(k,v)| (k, v.into_iter().map(|(k,_)| k).collect())).collect();
    }

    pub fn random_walk_over_graph_and_cut<T>(
        &mut self,
        args: &Cli,
        samples: usize,
        temperature: f64,
        steps: usize,
        walk_file: T,
        max_forward: usize,
        max_reads_forward: usize,
        safe_length_back: usize,
        ol_thresh: f64,
        tip_threshold: usize,
        strain_repeat_map: Option<&FxHashMap<NodeIndex, FxHashSet<NodeIndex>>>,
        beam: bool,
    ) where
        T: AsRef<std::path::Path>,
    {
        let mut walk_file = BufWriter::new(std::fs::File::create(walk_file).unwrap());
        let total_edgecount = Mutex::new(FxHashMap::default());
        let keys = self.nodes.keys().collect::<Vec<_>>();
        keys.into_par_iter().for_each(|node_id| {
            let edge_count;
            if beam {
                edge_count =
                    self.beam_search_path_prob(*node_id, temperature, 15, 5, 0.000000001, args.c, steps, steps);
            } else {
                edge_count = self.random_walk_sample(*node_id, samples, temperature, steps);
            }
            for (edge, count) in edge_count {
                let mut lockmap = total_edgecount.lock().unwrap();
                let total_count = lockmap.entry(edge).or_insert(0.);
                *total_count += count;
            }
        });

        let total_edgecount = total_edgecount.into_inner().unwrap();
        let mut min_by_nodeid_edgemap = FxHashMap::default();

        for edge_id in 0..self.edges.len() {
            if self.edges[edge_id].is_none() {
                continue;
            }
            let edge = self.edges[edge_id].as_ref().unwrap();

            let node1_val = total_edgecount
                .get(&(edge_id, edge.from_unitig))
                .unwrap_or(&0.);
            let node2_val = total_edgecount
                .get(&(edge_id, edge.to_unitig))
                .unwrap_or(&0.);
            let min_val = node1_val.min(*node2_val);
            min_by_nodeid_edgemap.insert(edge_id, min_val);
        }

        let mut sorted_edges = min_by_nodeid_edgemap.iter().collect::<Vec<_>>();
        sorted_edges.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        for (edge_id, count) in sorted_edges.iter() {
            let edge = self.edges[**edge_id].as_ref().unwrap();
            writeln!(
                walk_file,
                "u{},u{} has count: {}",
                self.nodes[&edge.from_unitig].node_id, self.nodes[&edge.to_unitig].node_id, count
            )
            .unwrap();
        }

        let mut removed_edges = FxHashSet::default();
        let always_cut_ol_thresh = 1.01; // TODO this is a hack that forces safely_cut_edges to remove this edge no matter what.
        for (edge_id, count) in sorted_edges.iter() {
            let edge = self.edges[**edge_id].as_ref().unwrap();
            let from_direction = edge.node_edge_direction(&edge.from_unitig);
            let from_edges = self.nodes[&edge.from_unitig].edges_direction(&from_direction);
            let to_direction = edge.node_edge_direction(&edge.to_unitig);
            let to_edges = self.nodes[&edge.to_unitig].edges_direction(&to_direction);

            for (i, edge_id_adj) in from_edges.iter().chain(to_edges.iter()).enumerate() {
                let mincount_adj = min_by_nodeid_edgemap.get(edge_id_adj).unwrap_or(&&0.);

                if (**count as f64) / (*mincount_adj) < ol_thresh {

                    let max_forward_adj;
                    let max_reads_forward_adj;

                    if **count / *mincount_adj < 0.01 {
                        max_forward_adj = 2 * max_forward;
                        max_reads_forward_adj = 2 * max_reads_forward;
                    } else {
                        max_forward_adj = max_forward;
                        max_reads_forward_adj = max_reads_forward;
                    }

                    self.safely_cut_edge(
                        **edge_id,
                        &mut removed_edges,
                        always_cut_ol_thresh,
                        None,
                        None,
                        strain_repeat_map,
                        false,
                        max_forward_adj,
                        max_reads_forward_adj,
                        safe_length_back,
                        &mut std::io::empty(),
                        args.c
                    );

                    //Check tip condition
                    if !removed_edges.contains(edge_id) {
                        let unitig;
                        if i < from_edges.len() {
                            unitig = &self.nodes[&edge.from_unitig];
                        } else {
                            unitig = &self.nodes[&edge.to_unitig];
                        }
                        let other_unitig = &self.nodes[&self.edges[**edge_id]
                            .as_ref()
                            .unwrap()
                            .other_node(unitig.node_hash_id)];
                        let direction_edge_other = self.edges[**edge_id]
                            .as_ref()
                            .unwrap()
                            .node_edge_direction(&other_unitig.node_hash_id)
                            .reverse();
                        let other_edges = other_unitig.edges_direction(&direction_edge_other);
                        if other_edges.len() == 0
                            && other_unitig.unique_length.unwrap() < tip_threshold
                        {
                            removed_edges.insert(**edge_id);
                        }
                    } else {
                        break;
                    }
                }
            }
        }
        self.remove_edges(removed_edges);
        self.re_unitig();
    }

    fn log_probability_ol_cov(
        &self,
        edge_id: EdgeIndex,
        start_unitig_id: NodeIndex,
        temperature: f64,
        c: usize,
        coverages: &[(MultiCov, usize)],
    ) -> f64 {
        let edge = self.edges[edge_id].as_ref().unwrap();
        let other_node = edge.other_node(start_unitig_id);
        let other_unitig = &self.nodes[&other_node];

        let other_covs = &other_unitig.read_min_depths_multi;
        //let cov_different = quantile_dist(coverages, other_covs);
        let cov_different = log_distribution_distance(coverages, other_covs);

        let start_unitig = &self.nodes[&start_unitig_id];
        let same_dir_edges =
            start_unitig.edges_direction(&edge.node_edge_direction(&start_unitig_id));

        let max_ol = same_dir_edges
            .iter()
            .map(|x| self.edges[*x].as_ref().unwrap().overlap.overlap_len_bases)
            .max()
            .unwrap();

        let ol = edge.overlap.overlap_len_bases;
        let ol_ratio = ol as f64 / max_ol as f64;
        let max_fsv = same_dir_edges
            .iter()
            .map(|x| self.edges[*x].as_ref().unwrap().edge_id_est(c))
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();

        // 99.5 - 100 -> e^(-1) 
        let fsv_diff = 200.0 * (edge.edge_id_est(c) - max_fsv);

        let cov1;
        let cov2;

        if edge.f1 {
            cov1 = self.nodes[&edge.from_unitig]
                .read_min_depths_multi
                .last()
                .unwrap();
        } else {
            cov1 = self.nodes[&edge.from_unitig]
                .read_min_depths_multi
                .first()
                .unwrap();
        }

        if edge.f2 {
            cov2 = self.nodes[&edge.to_unitig]
                .read_min_depths_multi
                .first()
                .unwrap();
        } else {
            cov2 = self.nodes[&edge.to_unitig]
                .read_min_depths_multi
                .last()
                .unwrap();
        }

        let cov_ratio_weight = 1.0
            - PSEUDOCOUNT
                / ((cov1.0[ID_THRESHOLD_ITERS - 1].min(cov2.0[ID_THRESHOLD_ITERS - 1])).max(1.0)
                    + PSEUDOCOUNT);

        let ol_term = cov_ratio_weight * (ol_ratio - 1.0);
        let fsv_term = fsv_diff;
        if let Some(cov_term) = cov_different {
            let edge_prob =
            (-cov_term + fsv_term + ol_term) / temperature;
            let edge_prob = edge_prob.max(-100.0);
            return edge_prob;
        } else {
            return -100.;
        }
    }

    fn beam_search_path_prob(
        &self,
        starting_node_id: NodeIndex,
        temperature: f64,
        _min_num_soln: usize,
        max_num_soln: usize,
        _min_prob_soln: f64,
        c: usize,
        min_reads_seen: usize,
        depth: usize,
    ) -> FxHashMap<(EdgeIndex, NodeIndex), f64> {
        let mut edge_count = FxHashMap::default();
        let starting_weight = (self.nodes[&starting_node_id].read_indices_ori.len() as f64).sqrt();
        let mut seen_nodes = FxHashSet::default();

        for direction in [Direction::Incoming, Direction::Outgoing] {
            let mut edges_to_search = vec![];
            let starting_edges = self.nodes[&starting_node_id].edges_direction(&direction);
            if starting_edges.len() == 0 {
                continue;
            }

            let mut starting_soln = BeamSearchSoln::default();
            starting_soln.coverages = self.nodes[&starting_node_id].read_min_depths_multi.clone();
            starting_soln.path_nodes = vec![starting_node_id];
            seen_nodes.insert(starting_node_id);
            starting_soln.depth = 0;

            edges_to_search.push((starting_soln, starting_edges));
            let mut current_depth = 0;
            let mut solutions = vec![];

            while current_depth < depth || edges_to_search.iter().any(|x| x.0.coverages.len() < min_reads_seen) {
                current_depth += 1;
                let mut candidate_soln = vec![];
                for (soln, edges) in edges_to_search.iter() {
                    for edge_id in edges.iter() {
                        let prob = self.log_probability_ol_cov(
                            *edge_id,
                            *soln.path_nodes.last().unwrap(),
                            temperature,
                            c, 
                            &soln.coverages,
                        );
                        let new_score = prob + soln.score;
                        candidate_soln.push((edge_id, new_score, soln));
                    }

                    if edges.len() == 0 {
                        solutions.push(soln.clone());
                    }
                }
                candidate_soln.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

                if candidate_soln.len() > max_num_soln {
                    candidate_soln.truncate(max_num_soln);
                }

                let mut new_edges_to_search = vec![];

                for (edge_id, new_score, soln) in candidate_soln {
                    let mut new_soln = soln.clone();
                    new_soln.score = new_score;
                    new_soln.path.push(*edge_id);
                    new_soln.depth += 1;

                    let other_node = self.edges[*edge_id]
                        .as_ref()
                        .unwrap()
                        .other_node(*soln.path_nodes.last().unwrap());
                    new_soln.path_nodes.push(other_node);
                    new_soln
                        .coverages
                        .extend(self.nodes[&other_node].read_min_depths_multi.clone());

                    let other_unitig = &self.nodes[&other_node];
                    let edge = &self.edges[*edge_id].as_ref().unwrap();
                    let new_edges = other_unitig
                        .edges_direction(&edge.node_edge_direction(&other_node).reverse());

                    new_edges_to_search.push((new_soln, new_edges));
                }

                edges_to_search = new_edges_to_search;
                if edges_to_search.is_empty(){
                    break;
                }
            }

            for (soln, _) in edges_to_search {
                solutions.push(soln);
            }

            let max_absolute_score = solutions
                .iter()
                .map(|x| x.score)
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap();
            let lse = solutions
                .iter()
                .map(|x| (x.score - max_absolute_score).exp())
                .sum::<f64>()
                .ln();
            let normalized_log_prob = solutions
                .iter()
                .map(|x| (x.score - max_absolute_score) - lse)
                .collect::<Vec<_>>();
            let probabilities = normalized_log_prob
                .iter()
                .map(|x| x.exp())
                .collect::<Vec<_>>();

            for (i, soln) in solutions.iter().enumerate() {
                let weight = starting_weight * probabilities[i];
                for (j, edge_id) in soln.path.iter().enumerate() {
                    let count = edge_count
                        .entry((*edge_id, soln.path_nodes[j]))
                        .or_insert(0.);
                    *count += weight;
                }
            }
        }
        return edge_count;
    }

    fn random_walk_sample(
        &self,
        starting_node_id: NodeIndex,
        samples: usize,
        temperature: f64,
        steps: usize,
    ) -> FxHashMap<(EdgeIndex, NodeIndex), f64> {
        let max_steps = 20;
        let mut edge_count = FxHashMap::default();
        let mut rng_seed = [0; 32];
        rng_seed[0] = (starting_node_id % 256) as u8;
        let mut rng = rand::rngs::StdRng::from_seed(rng_seed);
        let starting_weight = (self.nodes[&starting_node_id].read_indices_ori.len() as f64).sqrt();
        for _ in 0..samples {
            let mut used_edges = vec![];
            let mut coverages = self.nodes[&starting_node_id].read_min_depths_multi.clone();
            let mut left_unitig = Some(starting_node_id);
            let mut right_unitig = Some(starting_node_id);
            let mut left_direction = Direction::Incoming;
            let mut right_direction = Direction::Outgoing;
            let mut seen_edges = FxHashSet::default();
            for _ in 0..steps {
                // Get edge based on random traversal criteria if it exists
                let left_traverse = self.traverse_walk_probabilistic(
                    left_direction,
                    left_unitig,
                    temperature,
                    max_steps,
                    &mut coverages,
                    &mut rng,
                );
                let right_traverse = self.traverse_walk_probabilistic(
                    right_direction,
                    right_unitig,
                    temperature,
                    max_steps,
                    &mut coverages,
                    &mut rng,
                );

                if left_traverse.is_none() && right_traverse.is_none() {
                    break;
                }

                if let Some(left_edge_id) = left_traverse {
                    if seen_edges.contains(&left_edge_id) {
                        break;
                    } else {
                        seen_edges.insert(left_edge_id);
                    }
                    used_edges.push((left_edge_id, left_unitig.unwrap()));
                    let left_edge = self.edges[left_edge_id].as_ref().unwrap();
                    left_unitig = Some(left_edge.other_node(left_unitig.unwrap()));
                    left_direction = left_edge
                        .node_edge_direction(&left_unitig.unwrap())
                        .reverse();
                } else {
                    left_unitig = None;
                }
                if let Some(right_edge_id) = right_traverse {
                    if seen_edges.contains(&right_edge_id) {
                        break;
                    } else {
                        seen_edges.insert(right_edge_id);
                    }
                    used_edges.push((right_edge_id, right_unitig.unwrap()));
                    let right_edge = self.edges[right_edge_id].as_ref().unwrap();
                    right_unitig = Some(right_edge.other_node(right_unitig.unwrap()));
                    right_direction = right_edge
                        .node_edge_direction(&right_unitig.unwrap())
                        .reverse();
                } else {
                    right_unitig = None;
                }
            }

            for edge_and_node in used_edges {
                let count = edge_count.entry(edge_and_node).or_insert(0.);
                //*count += median_weight as f64 * starting_weight;
                *count += starting_weight;
            }
        }

        return edge_count;
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

            for i in 0..num_reads - 1 {
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

            for _ in 0..num_reads {
                let generic_read = TwinRead {
                    minimizers: vec![],
                    snpmers: vec![],
                    dna_seq: Seq::new(),
                    qual_seq: None,
                    base_length: 2000,
                    k: 21,
                    outer: false,
                    est_id: None,
                    id: "na".to_string(),
                    base_id: "na".to_string(),
                    split_start: 0,
                    min_depth_multi: Some([min_depth, min_depth, min_depth]),
                    median_depth: Some(2. * min_depth),
                    split_chimera: false,
                    snpmer_id_threshold: None,
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
                min_read_depth_multi: Some([min_depth; ID_THRESHOLD_ITERS]),
                median_read_depth: Some(2. * min_depth),
                unique_length: None,
                base_info: BaseInfo::default(),
                mapping_info: MappingInfo::default(),
                read_min_depths_multi: vec![([min_depth; ID_THRESHOLD_ITERS], 1000); num_reads],
                read_median_depths: vec![(min_depth, 1000); num_reads],
                alternate: false,
            };

            self.nodes.insert(node_id, node);
            node_id
        }

        // Add an edge between nodes
        fn add_edge(
            &mut self,
            from: NodeIndex,
            to: NodeIndex,
            overlap_len: usize,
            f1: bool,
            f2: bool,
        ) -> EdgeIndex {
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

            let edge_index = self.edges.len();
            self.edges.push(Some(edge));
            edge_index
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
        assert!(
            graph.nodes[&nodes_to_remove[0]]
                .min_read_depth_multi
                .unwrap()[0]
                == 5.0
        );
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
        assert!(
            graph.nodes[&nodes_to_remove[0]]
                .min_read_depth_multi
                .unwrap()[2]
                == 5.0
        );
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
        assert!(
            graph.nodes[&nodes_to_remove[0]]
                .min_read_depth_multi
                .unwrap()[0]
                == 5.0
        );
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
        assert!(
            graph.nodes[&nodes_to_remove[0]]
                .min_read_depth_multi
                .unwrap()[2]
                == 5.0
        );
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
            &FxHashSet::default(),
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
            &FxHashSet::default(),
        );

        assert!(!result);
    }

    //TODO this fails and it's hard to amend the algorithm to make this work... 
    fn _test_safe_path_detection_circ_tipbubble_small() {
        // Create a graph with a safe path:
        // n1 --- >n2 --> n3
        //   |(safe cut)|
        //       n4
        //       ||
        //       n5
        let mut builder = MockUnitigBuilder::new();

        let n1 = builder.add_node(500, 10.0);
        let n2 = builder.add_node(100, 10.0);
        let n3 = builder.add_node(500, 10.0);
        let n4 = builder.add_node(10, 10.0);
        let n5 = builder.add_node(10, 10.0);

        builder.add_edge(n1, n2, 10000, true, true);
        builder.add_edge(n2, n3, 10000, true, true);
        builder.add_edge(n4, n3, 10000, true, true);
        builder.add_edge(n1, n4, 1000, true, true);
        builder.add_edge(n4, n5, 1000, true, true);
        builder.add_edge(n5, n4, 1000, true, true);

        let (graph, _reads) = builder.build();

        // Test if path through n1->n5 is safe
        let edge_test = [(2, n3),(3, n1)];
        for (edge_id, nid) in edge_test{
            let edge = graph.edges[edge_id].as_ref().unwrap();
            let result = graph.safe_given_forward_back(
                &graph.nodes[&nid],
                edge,
                20000, // max_forward
                5,    // max_reads_forward
                100000, // safe_length_back
                None,
                &FxHashSet::default(),
            );
            assert!(result);
        }
    }

    #[test]
    fn test_safe_path_detection_circ_tip_small() {
        // Create a graph with a safe path:
        // n1 --- >n2
        // | (safe cut)
        // n3
        // (CIRC)
        let mut builder = MockUnitigBuilder::new();

        let n1 = builder.add_node(100, 10.0);
        let n2 = builder.add_node(100, 10.0);
        let n3 = builder.add_node(1, 10.0);

        builder.add_edge(n1, n2, 10000, true, true);
        builder.add_edge(n1, n3, 1000, true, true);
        builder.add_edge(n3, n3, 1000, true, true);

        let (graph, _reads) = builder.build();

        // Test if path through n1->n5 is safe
        let edge = graph.edges[1].as_ref().unwrap();
        let result = graph.safe_given_forward_back(
            &graph.nodes[&n1],
            edge,
            20000, // max_forward
            50,    // max_reads_forward
            10000, // safe_length_back
            None,
            &FxHashSet::default(),
        );

        assert!(result);
    }

    #[test]
    fn safely_cut_edge_test_basic_x() {
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

        for edge_id in edge_order {
            graph.safely_cut_edge(
                edge_id,
                &mut removed_edges,
                0.5,
                None,
                None,
                None,
                false,
                2000,
                5,
                1000,
                &mut unitig_edge_file,
                9
            );
        }

        assert!(removed_edges.contains(&0));
        assert!(!removed_edges.contains(&3));
    }

    #[test]
    fn safely_cut_edge_test_basic_x_fsv_safety() {
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

        let e1 = builder.add_edge(n1, n2, 10000, true, true);
        builder.add_edge(n1, n4, 2000, true, true);
        builder.add_edge(n3, n2, 2000, true, true);
        let e4 = builder.add_edge(n3, n4, 10000, true, true);
        builder.edges[e1].as_mut().unwrap().overlap.diff_snpmers = 10;
        builder.edges[e4].as_mut().unwrap().overlap.diff_snpmers = 10;

        let (graph, _reads) = builder.build();

        let mut removed_edges = FxHashSet::default();
        let mut unitig_edge_file = BufWriter::new(std::io::sink());

        let edge_order = vec![0, 3, 2, 1];

        for edge_id in edge_order.iter() {
            graph.safely_cut_edge(
                *edge_id,
                &mut removed_edges,
                0.5,
                None,
                None,
                None,
                true,
                2000,
                5,
                1000,
                &mut unitig_edge_file,
                9
            );
        }

        assert!(removed_edges.len() == 0);


        for edge_id in edge_order {
            graph.safely_cut_edge(
                edge_id,
                &mut removed_edges,
                0.5,
                None,
                None,
                None,
                false,
                2000,
                5,
                1000,
                &mut unitig_edge_file,
                9
            );
        }

        assert!(removed_edges.len() == 2);
    }

    #[test]
    fn safely_cut_edge_test_basic_x_cov() {
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

        for edge_id in edge_order {
            graph.safely_cut_edge(
                edge_id,
                &mut removed_edges,
                0.5,
                Some(3.),
                None,
                None,
                false,
                2000,
                5,
                1000,
                &mut unitig_edge_file,
                9
            );
        }

        assert!(removed_edges.contains(&0));
        assert!(!removed_edges.contains(&3));

        removed_edges.clear();

        // Edge 1 is cut due to cov, edge 0 is cut due to ol ratio and cov and is compatible
        let edge_order = vec![1, 2, 3, 0];

        for edge_id in edge_order {
            graph.safely_cut_edge(
                edge_id,
                &mut removed_edges,
                0.5,
                Some(3.),
                None,
                None,
                false,
                2000,
                5,
                1000,
                &mut unitig_edge_file,
                9
            );
        }

        dbg!(&removed_edges);

        assert!(removed_edges.contains(&0));
        assert!(removed_edges.contains(&2));
        assert!(!removed_edges.contains(&1));
        assert!(!removed_edges.contains(&3));
    }

    #[test]
    fn safely_cut_edge_test_dont_cut_if_supp_tip() {
        let mut builder = MockUnitigBuilder::new();

        // Make sure not to cut the "safe" edge if it is supported by a tip
        //
        // 1 ---> 2
        //  (3->2)  (small edge)
        // 3 ---> 4 (tip)
        //          10k

        let n1 = builder.add_node(100, 10.0);
        let n2 = builder.add_node(100, 10.0);
        let n3 = builder.add_node(100, 10.0);
        let n4 = builder.add_node(1, 10.0);

        builder.add_edge(n1, n2, 10000, true, true);
        builder.add_edge(n3, n4, 10000, true, true);
        builder.add_edge(n3, n2, 1000, true, true);

        let (graph, _reads) = builder.build();

        let mut removed_edges = FxHashSet::default();
        //let mut unitig_edge_file = BufWriter::new(std::fs::File::create("unitig_edge_file").unwrap());
        //create empty mock file
        let mut unitig_edge_file = BufWriter::new(std::io::sink());

        //Cut 0 (ol 1000) but don't cut 1 (ol 1000)
        let edge_order = vec![2, 0, 1];

        for edge_id in edge_order {
            graph.safely_cut_edge(
                edge_id,
                &mut removed_edges,
                0.5,
                Some(3.),
                None,
                None,
                false,
                2000,
                5,
                2000,
                &mut unitig_edge_file,
                9
            );
        }

        dbg!(&removed_edges);
        assert!(!removed_edges.contains(&2));

        removed_edges.clear();
    }

    #[test]
    fn safely_cut_edge_test_circular_contig_tip() {
        let mut builder = MockUnitigBuilder::new();

        // 1 ---> 1 (circular)
        // |  |
        // 2  3 (large)

        let n1 = builder.add_node(100, 100.0);
        let n2 = builder.add_node(1, 10.0);
        let n3 = builder.add_node(1, 100.0);

        builder.add_edge(n1, n1, 10000, true, true);
        builder.add_edge(n1, n2, 1000, true, true);
        builder.add_edge(n1, n3, 8000, true, true);

        let (graph, _reads) = builder.build();

        let mut removed_edges = FxHashSet::default();
        let mut unitig_edge_file = BufWriter::new(std::io::stdout());

        //Cut 0 (ol 1000) but don't cut 1 (ol 1000)
        let edge_order = vec![1, 0];

        for edge_id in edge_order {
            graph.safely_cut_edge(
                edge_id,
                &mut removed_edges,
                0.5,
                Some(3.),
                None,
                None,
                false,
                2000,
                5,
                1000,
                &mut unitig_edge_file,
                9
            );
        }

        assert!(removed_edges.contains(&1));
        assert!(!removed_edges.contains(&0));
        assert!(!removed_edges.contains(&2));
    }

    #[test]
    fn safely_cut_edge_test_circular_tip() {
        let mut builder = MockUnitigBuilder::new();

        // 1 ----> 2
        // |
        // 3
        // |
        // 4 (-> 4 self loop)
        //

        let n1 = builder.add_node(100, 100.0);
        let n2 = builder.add_node(100, 100.0);
        let n3 = builder.add_node(1, 10.0);
        let n4 = builder.add_node(2, 10.0);

        builder.add_edge(n1, n2, 10000, true, true);
        builder.add_edge(n1, n3, 1000, true, true);
        builder.add_edge(n3, n4, 3000, true, true);
        builder.add_edge(n4, n4, 3000, true, true);

        let (graph, _reads) = builder.build();

        let mut removed_edges = FxHashSet::default();
        let mut unitig_edge_file = BufWriter::new(std::io::stdout());

        //Cut 0 (ol 1000) but don't cut 1 (ol 1000)
        let edge_order = vec![0, 1];

        for edge_id in edge_order {
            graph.safely_cut_edge(
                edge_id,
                &mut removed_edges,
                0.5,
                Some(3.),
                None,
                None,
                false,
                20000,
                20,
                1000,
                &mut unitig_edge_file,
                9
            );
        }

        assert!(removed_edges.contains(&1));
        assert!(removed_edges.len() == 1)
    }

    #[test]
    fn test_strain_repeat_safety() {
        let mut builder = MockUnitigBuilder::new();

        // Create a graph with strain repeats:
        // 1 ---> 2
        //  |     |
        //  v     v
        // 3 ---> 4
        //  \     /
        //   v   v
        //    5

        let n1 = builder.add_node(100, 10.0);
        let n2 = builder.add_node(100, 10.0);
        let n3 = builder.add_node(100, 10.0);
        let n4 = builder.add_node(100, 10.0);
        let n5 = builder.add_node(100, 10.0);

        builder.add_edge(n1, n2, 10000, true, true);
        builder.add_edge(n1, n3, 10000, true, true);
        builder.add_edge(n2, n4, 10000, true, true);
        builder.add_edge(n3, n4, 10000, true, true);
        builder.add_edge(n3, n5, 10000, true, true);
        builder.add_edge(n4, n5, 10000, true, true);

        let (graph, _reads) = builder.build();

        // Create a strain repeat map
        let mut strain_repeat_map = FxHashMap::default();
        let mut repeats = FxHashSet::default();
        repeats.insert(n2);
        strain_repeat_map.insert(n3, repeats);
        let mut repeats = FxHashSet::default();
        repeats.insert(n4);
        strain_repeat_map.insert(n5, repeats);
        
        let mut removed_edges = FxHashSet::default();
        let mut unitig_edge_file = BufWriter::new(std::io::stdout());

        let edges = vec![0, 1, 2, 3, 4, 5];
        // Test safely_cut_edge with strain_repeat_safety
        for edge in edges{
            graph.safely_cut_edge(
                edge, // edge_id
                &mut removed_edges,
                1.01,
                None,
                None,
                Some(&strain_repeat_map),
                false,
                2000,
                5,
                1000000,
                &mut unitig_edge_file,
                9
            );
        }

        dbg!(&removed_edges, &strain_repeat_map);

        assert!(removed_edges.contains(&0));
        assert!(removed_edges.contains(&2));
        assert!(removed_edges.contains(&3));
        assert!(removed_edges.len() == 3);
    }

    #[test]
    fn random_walk_test_simple() {
        let mut builder = MockUnitigBuilder::new();

        // 1 ---> 2
        //  (1->4)
        // 3 ---> 4 (tip)

        let n1 = builder.add_node(10, 10.0);
        let n2 = builder.add_node(10, 10.0);
        let n3 = builder.add_node(10, 50.0);
        let n4 = builder.add_node(10, 50.0);

        builder.add_edge(n1, n2, 10000, true, true);
        builder.add_edge(n3, n4, 10000, true, true);
        builder.add_edge(n1, n4, 10000, true, true);

        let (graph, _reads) = builder.build();

        let temp = 1.;
        let samples = 100;
        let steps = 100;

        //This is fickle-- if RNG changes then we'll have to change this...
        let edge_counts = graph.random_walk_sample(n1, samples, temp, steps);
        assert!(edge_counts[&(0, n1)] > edge_counts[&(2, n1)]);
    }

    #[test]
    fn beam_test_simple() {
        let mut builder = MockUnitigBuilder::new();

        // 1 ---> 2
        //  (1->4)
        // 3 ---> 4 (tip)

        let n1 = builder.add_node(10, 10.0);
        let n2 = builder.add_node(10, 10.0);
        let n3 = builder.add_node(10, 50.0);
        let n4 = builder.add_node(10, 50.0);

        builder.add_edge(n1, n2, 10000, true, true);
        builder.add_edge(n3, n4, 10000, true, true);
        builder.add_edge(n1, n4, 10000, true, true);

        let (graph, _reads) = builder.build();

        let temp = 1.;

        //This is fickle-- if RNG changes then we'll have to change this...
        let edge_counts = graph.beam_search_path_prob(n1, temp, 5, 10, 0.00001, 9, 30, 5);
        dbg!(&edge_counts);
        assert!(edge_counts[&(0, n1)] > edge_counts[&(2, n1)]);
    }

    #[test]
    fn beam_walk_test_digraph() {
        let mut builder = MockUnitigBuilder::new();

        // 1 <---> 2
        //  (1<-<4)
        // 3 <---< 4 (tip)

        let n1 = builder.add_node(10, 10.0);
        let n2 = builder.add_node(10, 10.0);
        let n3 = builder.add_node(10, 50.0);
        let n4 = builder.add_node(10, 50.0);

        builder.add_edge(n1, n2, 10000, false, true);
        builder.add_edge(n3, n4, 10000, false, false);
        builder.add_edge(n1, n4, 10000, false, false);

        let (graph, _reads) = builder.build();

        let temp = 1.;
        let samples = 100;
        let steps = 100;

        let edge_counts = graph.random_walk_sample(n1, samples, temp, steps);
        dbg!(&edge_counts);
        assert!(edge_counts[&(0, n1)] > edge_counts[&(2, n1)]);

        let edge_counts = graph.beam_search_path_prob(n1, temp, 5, 10, 0.00001, 9, 50, 5);
        dbg!(&edge_counts);
        assert!(edge_counts[&(0, n1)] > edge_counts[&(2, n1)]);
    }

    #[test]
    fn random_walk_test_moderate() {
        let mut builder = MockUnitigBuilder::new();

        // 1 >---> 2 > ----> 3 >----> 4
        //    (1->9)           (10->4)
        //         9 >----->10 >----> 11
        //            (6->10)
        // 5 >---> 6 > ----> 7 >----> 8

        let n1 = builder.add_node(12, 10.0);
        let n2 = builder.add_node(10, 12.0);
        let n3 = builder.add_node(13, 10.0);
        let n4 = builder.add_node(15, 15.0);
        let n5 = builder.add_node(10, 50.0);
        let n6 = builder.add_node(10, 60.0);
        let n7 = builder.add_node(13, 40.0);
        let n8 = builder.add_node(15, 35.0);
        let n9 = builder.add_node(10, 1.0);
        let n10 = builder.add_node(10, 2.0);
        let n11 = builder.add_node(13, 25.0);

        let e1 = builder.add_edge(n1, n2, 10000, true, true);
        let e2 = builder.add_edge(n2, n3, 10000, true, true);
        let e3 = builder.add_edge(n3, n4, 10000, true, true);
        let e4 = builder.add_edge(n5, n6, 10000, true, true);
        let e5 = builder.add_edge(n6, n7, 10000, true, true);
        let e6 = builder.add_edge(n7, n8, 10000, true, true);
        let _e7 = builder.add_edge(n9, n10, 10000, true, true);

        let e8 = builder.add_edge(n10, n11, 10000, true, true); // 7
        let e9 = builder.add_edge(n1, n9, 10000, true, true); //8
        let e10 = builder.add_edge(n10, n4, 10000, true, true); // 9
        let e11 = builder.add_edge(n6, n10, 10000, true, true); // 10

        let (graph, _reads) = builder.build();

        let temp = 0.2;
        let samples = 100;
        let steps = 5;

        for j in 0..2 {
            let mut total_edgecounts = FxHashMap::default();
            for nodeid in 0..10 {
                let edge_counts;
                if j == 0 {
                    edge_counts = graph.random_walk_sample(nodeid, samples, temp, steps);
                } else {
                    edge_counts = graph.beam_search_path_prob(nodeid, temp, 5, 10, 0.00001, 9, 50, 5);
                }
                for (edge_id, count) in edge_counts {
                    let total_count = total_edgecounts.entry(edge_id).or_insert(0.);
                    *total_count += count;
                }
            }
            //print by sorted key
            let mut keys: Vec<_> = total_edgecounts.keys().collect();
            keys.sort();
            for key in keys {
                println!("Edge {:?} : {}", key, total_edgecounts[key]);
            }

            let bad_edge1 = (e11, n6);
            let bad_edge2 = (e9, n1);
            let bad_edge3 = (e10, n4);
            let bad_edge4 = (e8, n10);
            let good_edge = (e3, n3);

            for bad_edge in [bad_edge1, bad_edge2, bad_edge3, bad_edge4].iter() {
                assert!(
                    *total_edgecounts.get(bad_edge).unwrap_or(&0.)
                        < *total_edgecounts.get(&good_edge).unwrap_or(&0.) * 0.1
                );
            }

            let good_edges = vec![(e1, n1), (e2, n2), (e4, n5), (e5, n6), (e6, n7)];
            for good_edge in good_edges {
                assert!(*total_edgecounts.get(&good_edge).unwrap_or(&0.) > 0.01);
            }
        }
    }

    #[test]
    fn beam_walk_test_circular() {
        let mut builder = MockUnitigBuilder::new();

        //   e1      e2
        //n0 --> n1 --> n0 (circular)
        // |
        // |  e3
        // v
        // n2


        //Want:  w(e3) ~ w(e1)
        let n0 = builder.add_node(10, 20.0);
        let n1 = builder.add_node(10, 20.0);
        let n2 = builder.add_node(10, 20.0);

        let e1 = builder.add_edge(n0, n1, 5000, true, true);
        let e2 = builder.add_edge(n1, n0, 5000, true, true);
        let e3 = builder.add_edge(n0, n2, 5000, true, false);

        let (graph, _reads) = builder.build();

        let temp = 0.30;
        let samples = 10;
        let steps = 20;

        for j in 1..2 {
            let mut total_edgecounts = FxHashMap::default();
            for nodeid in 0..3 {
                let edge_counts;
                if j == 0 {
                    edge_counts = graph.random_walk_sample(nodeid, samples, temp, steps);
                } else {
                    edge_counts = graph.beam_search_path_prob(nodeid, temp, 5, 10, 0.00001,9, steps, steps);
                }

                for (edge_id, count) in edge_counts {
                    let total_count = total_edgecounts.entry(edge_id).or_insert(0.);
                    *total_count += count;
                }
            }
            //print by sorted key
            let mut keys: Vec<_> = total_edgecounts.keys().collect();
            keys.sort();
            for key in keys {
                println!("Edge {:?} : {}", key, total_edgecounts[key]);
            }

            let good_edge = (e2, n1);
            let opt1 = (e1,n0);
            let opt2 = (e3,n0);

            assert!(total_edgecounts[&good_edge] > total_edgecounts[&opt1]);
            assert!(total_edgecounts[&good_edge] > total_edgecounts[&opt2]);
            //TODO do we want this? Current model biases towards circularization from the topological defintion
            // of k-length paths. This isn't necessarily a _bad_ thing for metagenomics...
            assert!((total_edgecounts[&opt1] - total_edgecounts[&opt2]).abs() < 1.);
        }
    }

    #[test]
    fn cut_z_edge_true() {
        let mut builder = MockUnitigBuilder::new();

        // 1 ---> 2
        //  (1->4)
        // 3 ---> 4 (tip)

        let n1 = builder.add_node(10, 10.0);
        let n2 = builder.add_node(10, 10.0);
        let n3 = builder.add_node(10, 50.0);
        let n4 = builder.add_node(10, 50.0);

        builder.add_edge(n1, n2, 10000, true, true);
        builder.add_edge(n3, n4, 10000, true, true);
        builder.add_edge(n1, n4, 1000, true, true);

        let (graph, _reads) = builder.build();

        let unitig = &graph.nodes[&n1];
        let removed_edges = graph.get_z_edges(unitig, Direction::Outgoing);
        assert!(removed_edges.contains(&2));
        assert!(!removed_edges.contains(&0));
        assert!(!removed_edges.contains(&1));
    }

    //THIS TEST FAILS, BAD Z EDGE ALGO. TODO
    #[test]
    fn cut_z_edge_simple_bidir() {
        // 1 >---> 2
        //  (1>-<4)
        // 3 <---< 4 (tip)
        let mut builder = MockUnitigBuilder::new();
        let n1 = builder.add_node(10, 10.0);
        let n2 = builder.add_node(10, 10.0);
        let n3 = builder.add_node(10, 50.0);
        let n4 = builder.add_node(10, 50.0);

        builder.add_edge(n1, n2, 10000, true, true);
        builder.add_edge(n3, n4, 10000, true, false);
        builder.add_edge(n1, n4, 1000, true, false);

        let (graph, _reads) = builder.build();

        let unitig = &graph.nodes[&n1];
        let removed_edges = graph.get_z_edges(unitig, Direction::Outgoing);
        dbg!(&removed_edges);
        assert!(removed_edges.contains(&2));
        assert!(!removed_edges.contains(&0));
        assert!(!removed_edges.contains(&1));
    }

    // THIS TEST FAILS, BAD Z EDGE ALGO. TODO
    #[test]
    fn cut_z_edge_circular() {
        let mut builder = MockUnitigBuilder::new();

        // 1 --> 1 (circular)
        // |
        // 2

        let n1 = builder.add_node(10, 10.0);
        let n2 = builder.add_node(10, 10.0);

        builder.add_edge(n1, n1, 8000, true, true);
        builder.add_edge(n1, n2, 14000, true, true);

        let (graph, _reads) = builder.build();

        let unitig = &graph.nodes[&n1];
        let removed_edges = graph.get_z_edges(unitig, Direction::Outgoing);
        assert!(!removed_edges.contains(&0));
        assert!(!removed_edges.contains(&1));
    }

    #[test]
    fn cut_z_edge_circular_fp() {
        let mut builder = MockUnitigBuilder::new();

        // 1 --> 2 (circular)
        //   (1->4)
        //   (3->2)
        // 3 --> 4

        let n1 = builder.add_node(10, 10.0);
        let n2 = builder.add_node(10, 10.0);
        let n3 = builder.add_node(10, 10.0);
        let n4 = builder.add_node(10, 10.0);

        builder.add_edge(n1, n2, 14000, true, true);
        builder.add_edge(n1, n4, 5000, true, true);
        builder.add_edge(n3, n4, 14000, true, true);
        builder.add_edge(n3, n2, 5000, true, true);

        let (graph, _reads) = builder.build();

        let unitig = &graph.nodes[&n1];
        let removed_edges = graph.get_z_edges(unitig, Direction::Outgoing);
        assert!(!removed_edges.contains(&3));
        assert!(!removed_edges.contains(&2));
        assert!(!removed_edges.contains(&1));
        assert!(!removed_edges.contains(&0));
    }
}
