use bincode;
use flexi_logger::style;
use clap::Parser;
use flexi_logger::{DeferredNow, Duplicate, FileSpec, Record};
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use myloasm::cli;
use myloasm::constants::FORWARD_READ_SAFE_SEARCH_CUTOFF;
use myloasm::constants::ID_THRESHOLD_ITERS;
use myloasm::constants::MAX_BUBBLE_UNITIGS_FINAL_STAGE;
use myloasm::constants::POLISHED_CONTIGS_NAME;
use myloasm::constants::TS_DASHES_BLANK_COLONS_DOT_BLANK;
use myloasm::graph::GraphNode;
use myloasm::kmer_comp;
use myloasm::map_processing;
use myloasm::mapping;
use myloasm::polishing_mod;
use myloasm::seq_parse;
use myloasm::twin_graph;
use myloasm::twin_graph::OverlapConfig;
use myloasm::small_genomes;
use myloasm::types;
use myloasm::types::OverlapAdjMap;
use myloasm::unitig;
use myloasm::unitig::NodeSequence;
use myloasm::skani_dereplicate;
use myloasm::utils::*;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;
use tikv_jemallocator::Jemalloc;

#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

fn main() {
    let total_start_time = Instant::now();
    let mut args = cli::Cli::parse();
    let (output_dir, temp_dir) = initialize_setup(&mut args);

    log::info!("Starting assembly...");

    // Step 1: Process k-mers, count k-mers, and get SNPmers
    let mut kmer_info = get_kmers_and_snpmers(&args, &output_dir);
    log_memory_usage(true, "STAGE 1: Obtained SNPmers");

    // Step 2: Get twin reads using k-mer information
    let twin_read_container = get_twin_reads_from_kmer_info(&mut kmer_info, &args, &output_dir);
    let twin_reads = &twin_read_container.twin_reads;
    log_memory_usage(true, "STAGE 2: Obtained twin reads");

    // Step 3: Get overlaps between outer twin reads and construct raw unitig graph
    let overlaps = get_overlaps_from_twin_reads(&twin_read_container, &args, &output_dir);

    // Step 4: Construct raw unitig graph
    let (mut unitig_graph, overlap_adj_map) =
        unitig::UnitigGraph::from_overlaps(&twin_read_container.twin_reads, overlaps, Some(&twin_read_container.outer_indices), &args);

    let out_file = std::path::Path::new(&args.output_dir).join("unitig_graph.gfa");
    unitig_graph.to_gfa(out_file, true, true, &twin_reads, &args);

    // Step 5: First round of cleaning. Progressive cleaning: tips, bubbles, bridged repeats
    light_progressive_cleaning(&mut unitig_graph, &twin_reads, &args, &temp_dir, true);

    
    // Step 6: Second round of progressive cleaning; more aggressive.
    //heavy_cleaning(&mut unitig_graph, &twin_reads, &args, &temp_dir);
    heavy_clean_with_walk(
        &mut unitig_graph,
        &twin_reads,
        &overlap_adj_map,
        &args,
        &temp_dir,
        &output_dir,
    );

    // Step 6.5: Small circular contig retrieval
    small_genomes::two_cycle_retrieval(&mut unitig_graph, &twin_reads, &args, &temp_dir);

    // Step 7: Progressive coverage filter
    let mut contig_graph =
        progressive_coverage_contigs(unitig_graph, &twin_reads, &args, &temp_dir, &output_dir);

    // Step 8: Align reads back to graph. TODO consensus and etc.
    let mut get_seq_config = types::GetSequenceInfoConfig::default();
    get_seq_config.dna_seq_info = true;
    get_seq_config.blunted = false;
    contig_graph.get_sequence_info(&twin_reads, &get_seq_config);

    // Step 8.5: Dereplicate small contigs
    log::info!("Dereplicating spurious contigs...");
    mapping::map_to_dereplicate(&mut contig_graph, &kmer_info, &twin_reads, &args);
    
    log_memory_usage(true, "STAGE 3: Obtained unpolished contigs");
    contig_graph.print_statistics(&args);
    contig_graph.to_fasta(output_dir.join("final_contigs_nopolish.fa"), &args);

    if args.no_polish {
        log::warn!("No polishing requested. This is not recommended.");
        log::info!("Total time elapsed is {:?}", total_start_time.elapsed());
        return;
    }
    
    // Step 9: Align reads back to graph and take consensuses
    log::info!("Beginning final alignment of reads to graph...");
    let start = Instant::now();
    mapping::map_reads_to_unitigs(&mut contig_graph, &kmer_info, &twin_reads, &args);
    log_memory_usage(true, "STAGE 4: Mapped reads to contigs");

    log::info!(
        "Time elapsed for aligning reads to graph is {:?}",
        start.elapsed()
    );

    contig_graph.to_gfa(
        output_dir.join("final_contig_graph_noseq.gfa"),
        true,
        false,
        &twin_reads,
        &args,
    );

    log::info!("Polishing...");
    polishing_mod::polish_assembly(contig_graph, twin_read_container.twin_reads, &args);
    log::info!("Time elapsed for polishing is {:?}", start.elapsed());

    log::info!("Dereplicating polished contigs with skani...");
    skani_dereplicate::dereplicate_with_skani(POLISHED_CONTIGS_NAME, &args);

    log::info!("Total time elapsed is {:?}", total_start_time.elapsed());
}

fn my_own_format_colored(
    w: &mut dyn std::io::Write,
    now: &mut DeferredNow,
    record: &Record,
) -> Result<(), std::io::Error> {
    let mut paintlevel = record.level();
    if paintlevel == log::Level::Info {
        paintlevel = log::Level::Debug;
    }
    write!(
        w,
        "({}) {} [{}] {}",
        now.format(TS_DASHES_BLANK_COLONS_DOT_BLANK),
        style(paintlevel).paint(record.level().to_string()),
        record.module_path().unwrap_or(""),
        &record.args()
    )
}

fn my_own_format(
    w: &mut dyn std::io::Write,
    now: &mut DeferredNow,
    record: &Record,
) -> Result<(), std::io::Error> {
    write!(
        w,
        "({}) {} [{}] {}",
        now.format(TS_DASHES_BLANK_COLONS_DOT_BLANK),
        record.level(),
        record.module_path().unwrap_or(""),
        &record.args()
    )
}

fn initialize_setup(args: &mut cli::Cli) -> (PathBuf, PathBuf) {
    for file in &args.input_files {
        if !Path::new(file).exists() {
            eprintln!(
                "ERROR [myloasm] Input file {} does not exist. Exiting.",
                file
            );
            std::process::exit(1);
        }
    }

    let output_dir = Path::new(args.output_dir.as_str());
    let temp_dir = output_dir.join("temp");

    if !output_dir.exists() {
        std::fs::create_dir_all(output_dir).unwrap();
        std::fs::create_dir_all(output_dir.join("temp")).unwrap();
    } else {
        if !output_dir.is_dir() {
            eprintln!(
                "ERROR [myloasm] Output directory specified by `-o` exists and is not a directory."
            );
            std::process::exit(1);
        }

        if !temp_dir.exists() {
            std::fs::create_dir_all(&temp_dir).unwrap();
        } else {
            if !temp_dir.is_dir() {
                eprintln!(
                    "ERROR [myloasm] Could not create 'temp' directory within output directory."
                );
                std::process::exit(1);
            }
        }
    }

    // Initialize logger with CLI-specified level
    let log_spec = format!("{},skani=info", args.log_level_filter().to_string());
    let filespec = FileSpec::default()
        .directory(output_dir)
        .basename("myloasm");
    flexi_logger::Logger::try_with_str(log_spec)
        .expect("Something went wrong with logging")
        .log_to_file(filespec) // write logs to file
        .duplicate_to_stderr(Duplicate::Info) // print warnings and errors also to the console
        .format(my_own_format_colored) // use a simple colored format
        .format_for_files(my_own_format)
        .start()
        .expect("Something went wrong with creating log file");

    let cli_args: Vec<String> = std::env::args().collect();
    log::info!("COMMAND: {}", cli_args.join(" "));
    log::info!("VERSION: {}", env!("CARGO_PKG_VERSION"));

    // Validate k-mer size
    if args.kmer_size % 2 == 0 {
        log::error!("K-mer size must be odd");
        std::process::exit(1);
    }
    // Initialize thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    if args.hifi {
        args.snpmer_error_rate_lax = 0.001;
        args.snpmer_threshold_lax = 100.;
    }
    if args.r941 {
        args.snpmer_error_rate_lax = 0.05;
        args.contain_subsample_rate = 20;
    }

    return (output_dir.to_path_buf(), temp_dir.to_path_buf());
}

fn get_kmers_and_snpmers(args: &cli::Cli, output_dir: &PathBuf) -> types::KmerGlobalInfo {
    let empty_input = args.input_files.len() == 0;

    let kmer_info;
    if empty_input {
        if !output_dir.join("snpmer_info.bin").exists() {
            log::error!("No input files provided. See --help for usage.");
            std::process::exit(1);
        }
    }

    if empty_input && output_dir.join("snpmer_info.bin").exists() {
        kmer_info = bincode::deserialize_from(BufReader::new(
            File::open(output_dir.join("snpmer_info.bin")).unwrap(),
        ))
        .unwrap();
        log::info!("Loaded snpmer info from file.");
    } else {
        let start = Instant::now();
        let big_kmer_map = seq_parse::read_to_split_kmers(args.kmer_size, args.threads, &args);
        log::info!(
            "Time elapsed in for counting k-mers is: {:?}",
            start.elapsed()
        );

        let start = Instant::now();
        kmer_info = kmer_comp::get_snpmers(big_kmer_map, args.kmer_size, &args);
        log::info!(
            "Time elapsed in for parsing snpmers is: {:?}",
            start.elapsed()
        );
        bincode::serialize_into(
            BufWriter::new(File::create(output_dir.join("snpmer_info.bin")).unwrap()),
            &kmer_info,
        )
        .unwrap();
    }
    return kmer_info;
}

fn get_twin_reads_from_kmer_info(
    kmer_info: &mut types::KmerGlobalInfo,
    args: &cli::Cli,
    output_dir: &PathBuf,
) -> types::TwinReadContainer {
    let empty_input = args.input_files.len() == 0;
    let twin_read_container;

    if empty_input && output_dir.join("twin_reads.bin").exists() {
        twin_read_container = bincode::deserialize_from(BufReader::new(
            File::open(output_dir.join("twin_reads.bin")).unwrap(),
        ))
        .unwrap();
        log::info!("Loaded twin reads from file.");
    } else {
        // (1) read file, get twin reads
        log::info!("Getting twin reads from snpmers...");
        let twin_reads_raw = kmer_comp::twin_reads_from_snpmers(kmer_info, &args);
        let num_reads = twin_reads_raw.len();
        log_memory_usage(true, "STAGE 1.25: Initially obtained twin reads");

        // (2): removed contained reads
        log::info!("Removing contained reads - round 1...");
        let outer_read_indices_raw =
            twin_graph::remove_contained_reads_twin(None, &twin_reads_raw, &args);

        // (3): map all reads to the outer (non-contained) reads
        log::info!("Mapping reads to non-contained (outer) reads - round 1...");
        let start = Instant::now();
        let outer_mapping_info =
            mapping::map_reads_to_outer_reads(&outer_read_indices_raw, &twin_reads_raw, &args);
        log_memory_usage(true, "STAGE 1.5: Finished mapping reads to outer reads");

        // (4): split chimeric twin reads based on mappings to outer reads
        log::info!("Processing mappings...");
        let (mut split_twin_reads, _split_outer_read_indices) =
            map_processing::split_outer_reads(twin_reads_raw, outer_mapping_info, &args);
        log::info!(
            "Gained {} reads after splitting chimeras and mapping to outer reads in {:?}",
            split_twin_reads.len() as i64 - num_reads as i64,
            start.elapsed()
        );

        // (5): the splitted chimeric reads may be contained within the original reads, so we remove contained reads again
        log::info!("Mapping reads to non-contained (outer) reads - round 2...");
        let outer_read_indices = twin_graph::remove_contained_reads_twin(
            //Some(split_outer_read_indices),
            None,
            &split_twin_reads,
            &args,
        );

        // (6): map all reads to the outer new reads, which may be rescued or split chimeric reads
        let remap_indices = outer_read_indices
            .iter()
            .filter(|x| split_twin_reads[**x].min_depth_multi.is_none())
            .map(|x| *x)
            .collect::<Vec<_>>();
        log::info!(
            "Mapping reads to {} candidate outer reads...",
            remap_indices.len()
        );
        let chimeric_mapping_info =
            mapping::map_reads_to_outer_reads(&remap_indices, &split_twin_reads, &args);

        // (7): fix the depths of the chimeric reads
        for mapping_info in chimeric_mapping_info.iter() {
            let chimeric_index = mapping_info.tr_index;
            let split_chimeric_read = &mut split_twin_reads[chimeric_index];
            map_processing::populate_depth_from_map_info(
                split_chimeric_read,
                mapping_info,
                0,
                split_chimeric_read.base_length,
            );
        }

        split_twin_reads.iter_mut().for_each(|x| {
            if x.min_depth_multi.is_some() {
                x.outer = true;
            }
        });

        // (8): remove bad reads that have nothing mapped to them (this can bug for a variety of reasons)
        let outer_and_have_coverage_indices = outer_read_indices
            .iter()
            .filter(|x| split_twin_reads[**x].min_depth_multi.is_some())
            .map(|x| *x)
            .collect::<Vec<_>>();

        twin_read_container = types::TwinReadContainer {
            twin_reads: split_twin_reads,
            outer_indices: outer_and_have_coverage_indices,
            tig_reads: vec![],
        };

        bincode::serialize_into(
            BufWriter::new(File::create(output_dir.join("twin_reads.bin")).unwrap()),
            &twin_read_container,
        )
        .unwrap();
    }
    return twin_read_container;
}

fn get_overlaps_from_twin_reads(
    twin_read_container: &types::TwinReadContainer,
    args: &cli::Cli,
    output_dir: &PathBuf,
) -> Vec<OverlapConfig> {
    let empty_input = args.input_files.len() == 0;
    let twin_reads = &twin_read_container.twin_reads;
    let outer_read_indices = &twin_read_container.outer_indices;

    let overlaps;
    if empty_input && output_dir.join("overlaps.bin").exists() {
        overlaps = bincode::deserialize_from(BufReader::new(
            File::open(output_dir.join("overlaps.bin")).unwrap(),
        ))
        .unwrap();
        log::info!("Loaded overlaps from file.");
    } else {
        log::info!("Getting overlaps between outer reads...");
        let start = Instant::now();
        overlaps = twin_graph::get_overlaps_outer_reads_twin(
            &twin_reads,
            &outer_read_indices,
            &args,
            Some("overlaps.txt"),
        );
        log::info!(
            "Time elapsed for getting overlaps is: {:?}",
            start.elapsed()
        );
        bincode::serialize_into(
            BufWriter::new(File::create(output_dir.join("overlaps.bin")).unwrap()),
            &overlaps,
        )
        .unwrap();
    }
    return overlaps;
}

fn light_progressive_cleaning(
    unitig_graph: &mut unitig::UnitigGraph,
    twin_reads: &Vec<types::TwinRead>,
    args: &cli::Cli,
    temp_dir: &PathBuf,
    output_temp: bool,
) {
    log::info!("Initial light progressive cleaning...");
    let get_seq_config = types::GetSequenceInfoConfig::default();

    let mut iteration = 1;
    let divider = 3;
    let max_attempts = 10;
    let max_dropcut_thresh = 0.50;

    //let safety_edge_cov_score_thresholds = [50., 25., 10.];
    let safety_edge_cov_score_thresholds = [1000000.];
    let bubble_length_cutoff = args.small_bubble_threshold;

    loop {
        log::debug!("LIGHT-CLEAN: Cleaning graph iteration {}", iteration);
        let mut size_graph = unitig_graph.nodes.len();
        let mut counter = 0;
        let tip_length_cutoff = args.tip_length_cutoff;
        let read_cutoff = args.tip_read_cutoff;

        //First iteration, with spurious haplotype edge removal
        unitig_graph.remove_tips(tip_length_cutoff, read_cutoff, false);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.pop_bubbles(bubble_length_cutoff, None, false);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.remove_low_id_haplotype_edges(&args);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.remove_tips(tip_length_cutoff, read_cutoff, false);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);

        // unitig_graph.remove_singleton_lowcov_nodes(&args);
        // unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);

        // Remove tips
        loop {
            unitig_graph.remove_tips(tip_length_cutoff, read_cutoff, false);
            unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
            //unitig_graph.remove_caps();
            //unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
            unitig_graph.pop_bubbles(bubble_length_cutoff, None, false);
            unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);

            if unitig_graph.nodes.len() == size_graph {
                break;
            }

            size_graph = unitig_graph.nodes.len();
            counter += 1;
            if counter == max_attempts {
                break;
            }
        }

        if output_temp {
            unitig_graph.to_gfa(
                temp_dir.join(format!("{}-tip_unitig_graph.gfa", iteration)),
                true,
                true,
                &twin_reads,
                &args,
            );
        }

        // Pop bubbles
        unitig_graph.pop_bubbles(bubble_length_cutoff, None, false);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        if output_temp {
            unitig_graph.to_gfa(
                temp_dir.join(format!("{}-bubble_unitig_graph.gfa", iteration)),
                true,
                false,
                &twin_reads,
                &args,
            );
        }

        //Cut bridged repeats; "drop" cuts
        log::debug!("LIGHT-CLEAN: Resolving bridged repeats...");
        let prebridge_file = temp_dir.join(format!("{}-pre_bridge_cuts.txt", iteration));
        let edge_safe_cov_threshold = safety_edge_cov_score_thresholds
            [(iteration - 1).min(safety_edge_cov_score_thresholds.len() - 1)];
        unitig_graph.resolve_bridged_repeats(
            &args,
            max_dropcut_thresh / (divider as f64) * iteration as f64,
            None,
            Some(edge_safe_cov_threshold),
            prebridge_file,
            FORWARD_READ_SAFE_SEARCH_CUTOFF,
            args.tip_read_cutoff,
            100_000,
        );
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        if output_temp {
            unitig_graph.to_gfa(
                temp_dir.join(format!("{}-resolve_unitig_graph.gfa", iteration)),
                true,
                false,
                &twin_reads,
                &args,
            );
        }

        iteration += 1;
        if iteration == divider + 1 {
            break;
        }
    }

    let mut size_graph = unitig_graph.nodes.len();
    let mut counter = 0;
    loop {
        unitig_graph.remove_tips(args.tip_length_cutoff, args.tip_read_cutoff, false);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.pop_bubbles(args.small_bubble_threshold, None, false);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        if unitig_graph.nodes.len() == size_graph {
            break;
        }
        size_graph = unitig_graph.nodes.len();
        counter += 1;
        if counter == max_attempts {
            break;
        }
    }
    if output_temp {
        unitig_graph.to_gfa(
            temp_dir.join("after-light.gfa"),
            true,
            false,
            &twin_reads,
            &args,
        );
    }
}

fn get_contigs_from_progressive_coverage(
    mut graph_per_coverage: FxHashMap<usize, unitig::UnitigGraph>,
    _args: &cli::Cli,
    _temp_dir: &PathBuf,
) -> unitig::UnitigGraph {
    let mut final_unitig_graph = unitig::UnitigGraph::new();
    let mut used_reads = FxHashSet::default();
    let mut covs = graph_per_coverage.keys().cloned().collect::<Vec<_>>();
    covs.sort();

    // Iterate over indices separately to avoid borrowing conflicts
    for cov in covs.into_iter().rev() {
        let mut contigs_to_keep = FxHashSet::default();

        // First pass: identify contigs to keep
        for contig in graph_per_coverage[&cov].nodes.values() {
            if (contig.min_read_depth_multi.unwrap().iter().sum::<f64>()
                / (ID_THRESHOLD_ITERS as f64))
                >= cov as f64 * 2.0
            {
                let mut unused = true;

                for (read, _) in contig.read_indices_ori.iter() {
                    if used_reads.contains(read) {
                        unused = false;
                        break;
                    }
                }

                if !unused {
                    continue;
                }

                contigs_to_keep.insert(contig.node_hash_id.clone());

                log::trace!(
                    "Keeping contig u{} with covs {:?} of size {} at threshold {}",
                    contig.node_id,
                    contig.min_read_depth_multi.unwrap(),
                    contig.cut_length(),
                    cov
                );

                for (read, _) in contig.read_indices_ori.iter() {
                    used_reads.insert(read.clone());
                }
            }
        }

        // Second pass: remove nodes
        let contigs_to_remove: Vec<_> = graph_per_coverage[&cov]
            .nodes
            .keys()
            .filter(|x| !contigs_to_keep.contains(x))
            .cloned()
            .collect();

        graph_per_coverage
            .get_mut(&cov)
            .unwrap()
            .remove_nodes(&contigs_to_remove, false);

        //Have to do work in order to preserve indices during updates
        let mut new_node_hash_id_map = FxHashMap::default();
        for (count, node) in graph_per_coverage
            .get_mut(&cov)
            .unwrap()
            .nodes
            .values_mut()
            .enumerate()
        {
            for edge_index in node.in_edges_mut().iter_mut() {
                *edge_index += final_unitig_graph.edges.len();
            }
            for edge_index in node.out_edges_mut().iter_mut() {
                *edge_index += final_unitig_graph.edges.len();
            }

            let new_hash_id = count + final_unitig_graph.nodes.len();
            new_node_hash_id_map.insert(node.node_hash_id, new_hash_id);
            node.node_hash_id = new_hash_id;
        }

        for edge_index in 0..graph_per_coverage[&cov].edges.len() {
            let edge_opt = &mut graph_per_coverage.get_mut(&cov).unwrap().edges[edge_index];
            if edge_opt.is_none() {
                continue;
            }
            let edge = &mut edge_opt.as_mut().unwrap();
            let from_unitig_clone = edge.from_unitig.clone();
            let to_unitig_clone = edge.to_unitig.clone();
            edge.from_unitig = new_node_hash_id_map[&from_unitig_clone];
            edge.to_unitig = new_node_hash_id_map[&to_unitig_clone];
        }

        final_unitig_graph.edges.extend(std::mem::take(
            &mut graph_per_coverage.get_mut(&cov).unwrap().edges,
        ));

        for (_, node) in std::mem::take(&mut graph_per_coverage.get_mut(&cov).unwrap().nodes) {
            final_unitig_graph.nodes.insert(node.node_hash_id, node);
        }
    }

    //Invalidate all non-circular edges -- otherwise we get weird joining artefacts.
    let mut remove_edge_set = FxHashSet::default();
    for (id, edge) in final_unitig_graph.edges.iter().enumerate() {
        if edge.is_none() {
            continue;
        }
        let edge = edge.as_ref().unwrap();
        if !final_unitig_graph.nodes[&edge.from_unitig].is_circular()
            || !final_unitig_graph.nodes[&edge.to_unitig].is_circular()
        {
            remove_edge_set.insert(id);
        }
    }

    final_unitig_graph.remove_edges(remove_edge_set);
    final_unitig_graph.re_unitig();
    return final_unitig_graph;
}

fn heavy_clean_with_walk(
    unitig_graph: &mut unitig::UnitigGraph,
    twin_reads: &Vec<types::TwinRead>,
    overlap_adj_map: &OverlapAdjMap,
    args: &cli::Cli,
    temp_dir: &PathBuf,
    out_dir: &PathBuf,
) {
    log::info!("Graph cleaning with walks...");
    let mut save_removed;
    let temperatures = [2., 1.5, 1.0, 0.5];
    let ol_thresholds = [0.125, 0.25, 0.5];
    let aggressive_multipliers = [10, 15, 30];
    let mut global_counter = 0;
    let mut size_graph = unitig_graph.nodes.len();
    let mut special_small;
    let samples = 20;
    let steps = 10;

    let get_seq_config = types::GetSequenceInfoConfig::default();

    //Try switching temp/multiplier
    for temperature in temperatures {
        for multiplier in aggressive_multipliers {
            //TODO 
            if multiplier > 0 {
                save_removed = true;
            }
            else{
                save_removed = false;
            }

            if multiplier > 20 {
                special_small = true;
            } else {
                special_small = false;
            }

            let mut counter = 0;
            loop {
                let tip_length_cutoff_heavy = args.tip_length_cutoff * 5;
                let tip_read_cutoff_heavy = args.tip_read_cutoff * 5;
                let bubble_threshold_heavy =
                    (args.small_bubble_threshold * multiplier).min(1_000_000);
                remove_tips_until_stable(
                    unitig_graph,
                    twin_reads,
                    tip_length_cutoff_heavy,
                    tip_read_cutoff_heavy,
                    bubble_threshold_heavy,
                    usize::MAX,
                    Some(5),
                    temp_dir,
                    save_removed,
                    args,
                );

                let ind = counter.min(ol_thresholds.len() - 1);
                let ol_threshold = ol_thresholds[ind];

                unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
                if counter <= (ol_thresholds.len() - 1) && counter != 0 {
                    unitig_graph.to_gfa(
                        temp_dir.join(format!(
                            "heavy-{}-{}-{}.gfa",
                            multiplier, temperature, ol_threshold
                        )),
                        true,
                        false,
                        &twin_reads,
                        &args,
                    );
                }

                let strain_repeats = unitig_graph.get_strain_repeats(&overlap_adj_map, &args);

                let length_before_cut = unitig_graph.nodes.len();
                unitig_graph.random_walk_over_graph_and_cut(
                    &args,
                    samples,
                    temperature,
                    steps,
                    temp_dir.join(format!("walk-edge-{}.txt", global_counter)),
                    args.tip_length_cutoff * multiplier,
                    args.tip_read_cutoff * multiplier,
                    300_000,
                    ol_threshold,
                    args.tip_length_cutoff * multiplier,
                    Some(&strain_repeats),
                    true,
                    special_small
                );

                let length_after_cut = unitig_graph.nodes.len();
                log::debug!(
                    "WALK-CLEAN: Cut {} nodes at multiplier {}, temperature {}, and threshold {}",
                    length_before_cut - length_after_cut,
                    multiplier,
                    temperature,
                    ol_threshold
                );
                unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);

                counter += 1;
                global_counter += 1;

                if unitig_graph.nodes.len() == size_graph && counter >= ol_thresholds.len() {
                    break;
                }
                size_graph = unitig_graph.nodes.len();
            }
        }
    }

    unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
    unitig_graph.to_gfa(
        out_dir.join("after_walk.gfa"),
        true,
        true,
        &twin_reads,
        &args,
    );
}

fn _heavy_cleaning(
    unitig_graph: &mut unitig::UnitigGraph,
    twin_reads: &Vec<types::TwinRead>,
    args: &cli::Cli,
    temp_dir: &PathBuf,
) {
    log::info!("Heavy graph cleaning...");
    let get_seq_config = types::GetSequenceInfoConfig::default();
    let mut size_graph = unitig_graph.nodes.len();
    let mut counter = 0;
    //let cov_score_thresholds = [20., 10., 5., 3., 2.];
    let cov_score_thresholds = [20., 10.];
    //let ratio_length_thresholds = [0.25, 0.5, 0.75];
    let ratio_length_thresholds = [0.25, 0.5];
    let safety_edge_cov_score_thresholds = [100000000.];

    let save_tips = false;
    loop {
        let ind = counter.min(cov_score_thresholds.len() - 1);
        let ind_ratio = counter.min(ratio_length_thresholds.len() - 1);
        let ind_edge_safety_ratio = counter.min(safety_edge_cov_score_thresholds.len() - 1);
        let tip_length_cutoff_heavy = args.tip_length_cutoff * 5;
        let tip_read_cutoff_heavy = args.tip_read_cutoff;
        let bubble_threshold_heavy = args.max_bubble_threshold;

        remove_tips_until_stable(
            unitig_graph,
            twin_reads,
            tip_length_cutoff_heavy,
            tip_read_cutoff_heavy,
            bubble_threshold_heavy,
            usize::MAX,
            None,
            temp_dir,
            save_tips,
            args,
        );

        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.to_gfa(
            temp_dir.join(format!("heavy-{}-clean_unitig_graph.gfa", counter)),
            true,
            false,
            &twin_reads,
            &args,
        );

        unitig_graph.resolve_bridged_repeats(
            &args,
            ratio_length_thresholds[ind_ratio],
            Some(cov_score_thresholds[ind] as f64),
            Some(safety_edge_cov_score_thresholds[ind_edge_safety_ratio]),
            temp_dir.join(format!("f{}-resolve_unitig_graph.txt", counter)),
            args.tip_length_cutoff * 5,
            args.tip_read_cutoff * 5,
            300_000,
        );
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.to_gfa(
            temp_dir.join(format!("heavy-{}-resolve_unitig_graph.gfa", counter)),
            true,
            false,
            &twin_reads,
            &args,
        );
        if unitig_graph.nodes.len() == size_graph {
            break;
        }
        size_graph = unitig_graph.nodes.len();
        counter += 1;
    }

    unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);

    unitig_graph.to_gfa(
        temp_dir.join("after_heavy.gfa"),
        true,
        true,
        &twin_reads,
        &args,
    );
}

fn remove_tips_until_stable(
    unitig_graph: &mut unitig::UnitigGraph,
    twin_reads: &Vec<types::TwinRead>,
    tip_length_cutoff: usize,
    tip_read_cutoff: usize,
    max_bubble_threshold: usize,
    max_bubble_tigs: usize,
    max_attempts: Option<usize>,
    _temp_dir: &PathBuf,
    save_tips: bool,
    _args: &cli::Cli,
) {
    let get_seq_config = types::GetSequenceInfoConfig::default();
    let mut size_graph = unitig_graph.nodes.len();
    let mut counter = 0;
    loop {
        if let Some(max_attempts) = max_attempts {
            if counter == max_attempts {
                break;
            }
        }
        unitig_graph.remove_tips(tip_length_cutoff, tip_read_cutoff, save_tips);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.pop_bubbles(max_bubble_threshold, Some(max_bubble_tigs), save_tips);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        //unitig_graph.cut_z_edges_circular_only(args);
        //unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        if unitig_graph.nodes.len() == size_graph {
            break;
        }
        size_graph = unitig_graph.nodes.len();
        counter += 1;
    }
}

fn progressive_coverage_contigs(
    unitig_graph: unitig::UnitigGraph,
    twin_reads: &Vec<types::TwinRead>,
    args: &cli::Cli,
    temp_dir: &PathBuf,
    output_dir: &PathBuf,
) -> unitig::UnitigGraph {
    log::info!("Progressive coverage filtering...");
    let prog_dir = temp_dir.join("progressive");

    let max_cov = unitig_graph
        .nodes
        .values()
        .map(|x| x.min_read_depth_multi.unwrap().iter().sum::<f64>() / ID_THRESHOLD_ITERS as f64)
        .max_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap_or(1.);

    let tip_length_cutoff_heavy = args.tip_length_cutoff * 15;
    let tip_read_cutoff_heavy = args.tip_read_cutoff;
    let bubble_threshold_heavy = 6_000_000;

    let mut cov_to_graph_map = FxHashMap::default();
    let mut unitig_graph_with_threshold = unitig_graph.clone();
    remove_tips_until_stable(
        &mut unitig_graph_with_threshold,
        twin_reads,
        tip_length_cutoff_heavy,
        tip_read_cutoff_heavy,
        bubble_threshold_heavy,
        usize::MAX,
        None,
        &temp_dir,
        true,
        &args,
    );

    let up_to_30 = (0..30).collect::<Vec<_>>();
    let after_30 = (30..=max_cov as usize).step_by(2).collect::<Vec<_>>();
    let all_covs = up_to_30.iter().chain(after_30.iter());

    for &cov_thresh in all_covs.into_iter() {
        let cov_thresh = cov_thresh as f64;

        unitig_graph_with_threshold.cut_coverage(cov_thresh);
        unitig_graph_with_threshold
            .get_sequence_info(&twin_reads, &types::GetSequenceInfoConfig::default());

        unitig_graph_with_threshold.cut_z_edges(args);
        unitig_graph_with_threshold
            .get_sequence_info(&twin_reads, &types::GetSequenceInfoConfig::default());

        remove_tips_until_stable(
            &mut unitig_graph_with_threshold,
            twin_reads,
            tip_length_cutoff_heavy,
            tip_read_cutoff_heavy,
            bubble_threshold_heavy,
            MAX_BUBBLE_UNITIGS_FINAL_STAGE,
            None,
            &temp_dir,
            true,
            &args,
        );

        unitig_graph_with_threshold
            .get_sequence_info(&twin_reads, &types::GetSequenceInfoConfig::default());

        if cov_thresh < 50. {
            let prog_dir_cov = prog_dir.join(format!("cov_{}", cov_thresh));
            std::fs::create_dir_all(&prog_dir_cov).unwrap();
            unitig_graph_with_threshold.to_gfa(
                prog_dir_cov.join("filtered_graph.gfa"),
                true,
                true,
                &twin_reads,
                &args,
            );
        }

        // Push contigs at this coverage level
        cov_to_graph_map.insert(cov_thresh as usize, unitig_graph_with_threshold.clone());
    }

    // Includes end and beginning node.
    let mut unitig_graph =
        get_contigs_from_progressive_coverage(cov_to_graph_map, &args, &temp_dir);
    unitig_graph.get_sequence_info(&twin_reads, &types::GetSequenceInfoConfig::default());
    unitig_graph.to_gfa(
        output_dir.join("fnc_nomap.gfa"),
        true,
        true,
        &twin_reads,
        &args,
    );

    return unitig_graph;
}


