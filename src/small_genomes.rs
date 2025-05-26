use crate::unitig;
use crate::types::*;
use crate::cli;
use crate::twin_graph;
use crate::types;
use std::path::PathBuf;
use fxhash::FxHashMap;
use crate::twin_graph::OverlapConfig;


pub fn two_cycle_retrieval(
    unitig_graph: &mut unitig::UnitigGraph,
    twin_reads: &Vec<types::TwinRead>,
    args: &cli::Cli,
    _temp_dir: &PathBuf,
) {
    log::info!("Small circular contig extraction...");
    let node_to_component_char = unitig_graph.get_all_connected_components(true);
    let mut component_to_nodes_map = FxHashMap::default();

    // Extract components that are small
    for (node_id, component_char) in node_to_component_char.iter() {
        let component_max_length = component_char.length;
        let component_num_reads = component_char.num_reads;
        let component_id = component_char.component_num;
        //TODO - make 100_000 depend on read length.
        if component_max_length < 100_000 && component_num_reads < 100 {
            component_to_nodes_map
                .entry(component_id)
                .or_insert_with(Vec::new)
                .push(*node_id);
        }
    }

    let mut nodes_to_remove = vec![];
    let mut good_cycles: Vec<Vec<OverlapConfig>> = vec![];

    for component in component_to_nodes_map.values() {

        let mut read_ids = vec![];

        //Extract reads
        for node_id in component.iter() {
            let node = unitig_graph.nodes.get(node_id).unwrap();
            for (read_id, _) in node.read_indices_ori.iter() {
                read_ids.push(*read_id);
            }
        }

        log::trace!(
            "Component with unitigs {:?} has reads {:?}",
            &component,
            &read_ids
        );

        if read_ids.len() == 1{
            continue;
        }

        //Extract overlaps over the component
       // let file_pathbuf = format!("temp/{}.small_circular", component[0]);
        let mut relaxed_args = args.clone();
        relaxed_args.min_ol = args.min_ol / 2;
        let component_overlaps =
            twin_graph::get_overlaps_outer_reads_twin(&twin_reads, &read_ids, &relaxed_args, None, None);

        // Find 2-cycle overlaps
        let adjacency_map = component_overlaps
            .iter()
            .fold(FxHashMap::default(), |mut acc, x| {
                if x.read_i < x.read_j {
                    acc.entry(x.read_i).or_insert_with(Vec::new).push(x);
                } else {
                    acc.entry(x.read_j).or_insert_with(Vec::new).push(x);
                }
                acc
            });

        let mut two_cycles = vec![];
        for (_, mut overlaps) in adjacency_map.into_iter() {
            overlaps.sort_by(|a, b| (a.read_i, a.read_j).cmp(&(b.read_i, b.read_j)));
            let mut prev: Option<&OverlapConfig> = None;
            for overlap in overlaps {
                if prev.is_some() {
                    let ol1 = prev.unwrap();
                    if (ol1.read_i, ol1.read_j) == (overlap.read_i, overlap.read_j)
                        || (ol1.read_i, ol1.read_j) == (overlap.read_j, overlap.read_i)
                    {
                        log::trace!(
                            "Found 2-cycle overlap between {} and {}",
                            overlap.read_i,
                            overlap.read_j
                        );
                        two_cycles.push((ol1, overlap));
                    }
                }
                prev = Some(overlap);
            }
        }

        if two_cycles.is_empty(){
            continue;
        }

        let mut cycle_stats = vec![];
        for (ol1, ol2) in two_cycles.into_iter() {
            let total_length =
                twin_reads[ol1.read_i].base_length + twin_reads[ol1.read_j].base_length;
            let total_overlap =
                (ol1.overlap1_len + ol1.overlap2_len + ol2.overlap1_len + ol2.overlap2_len) / 2;
            let total_hangs = ol1.hang1 + ol1.hang2 + ol2.hang1 + ol2.hang2;

            let circular_length = total_length - total_overlap - total_hangs;

            let hang1_penalty = ol1.hang1 as i64 - ol1.hang2 as i64;
            let hang2_penalty = ol2.hang1 as i64 - ol2.hang2 as i64;
            let total_mini = ol1.shared_mini + ol2.shared_mini;

            let two_cycle_stat = TwoCycle {
                read_i: ol1.read_i,
                read_j: ol1.read_j,
                circular_length,
                hang_penalty: hang1_penalty.abs() + hang2_penalty.abs(),
                total_mini,
            };

            cycle_stats.push((two_cycle_stat, vec![ol1.clone(), ol2.clone()]));
        }

        cycle_stats.sort_by(|a, b| a.0.circular_length.cmp(&b.0.circular_length));

        let mut prev_cycle_length_opt = None;
        let mut optimal_cycle: Option<&TwoCycle> = None;
        let mut optimal_index: Option<usize> = None;

        for (cycle, _) in cycle_stats.iter() {
            log::trace!(
                "Component {} Cycle {}-{}, Circular length: {}, Hang penalty: {}, Total mini: {}",
                component[0],
                cycle.read_i,
                cycle.read_j,
                cycle.circular_length,
                cycle.hang_penalty,
                cycle.total_mini
            );
        }

        for (i, (cycle, _)) in cycle_stats.iter().enumerate() {
            if let Some(prev_cycle_length) = prev_cycle_length_opt {
                if cycle.circular_length as f64 / prev_cycle_length as f64 > 1.10 {
                    log::trace!(
                        "SMALL CYCLES: Optimal cycle with length {} found",
                        optimal_cycle.unwrap().circular_length
                    );
                    good_cycles.push(cycle_stats[optimal_index.unwrap()].1.clone());
                    optimal_cycle = None;
                    optimal_index = None;
                    prev_cycle_length_opt = None;
                    continue;
                }
            }
            if let Some(optimal_cycle_un) = optimal_cycle {
                if cycle.hang_penalty < optimal_cycle_un.hang_penalty{
                    optimal_cycle = Some(cycle);
                    optimal_index = Some(i);
                }
            }
            else{
                optimal_cycle = Some(cycle);
                optimal_index = Some(i);
            }
            prev_cycle_length_opt = Some(cycle.circular_length);
        }

        if let Some(optimal_cycle) = optimal_cycle {
            log::trace!(
                "SMALL CYCLES: Optimal cycle with length {} found",
                optimal_cycle.circular_length
            );
            good_cycles.push(cycle_stats[optimal_index.unwrap()].1.clone());
        }

        if !good_cycles.is_empty(){
            nodes_to_remove.extend(component.iter().map(|x| *x));
        }
    }

    unitig_graph.remove_nodes(&nodes_to_remove, false);

    log::debug!("Number of small cycles found: {}", good_cycles.len());

    for ols in good_cycles.into_iter(){
        let (new_graph, _) = unitig::UnitigGraph::from_overlaps(&twin_reads, ols, None, args);
        unitig_graph.concatenate(new_graph);
    }

    unitig_graph.get_sequence_info(&twin_reads, &types::GetSequenceInfoConfig::default());
}
