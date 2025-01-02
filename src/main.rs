use bincode;
use clap::Parser;
use fxhash::FxHashSet;
use myloasm::cli;
use myloasm::consensus;
use myloasm::constants::FORWARD_READ_SAFE_SEARCH_CUTOFF;
use myloasm::kmer_comp;
use myloasm::map_processing;
use myloasm::mapping;
use myloasm::seq_parse;
use myloasm::twin_graph;
use myloasm::types;
use myloasm::types::TwinOverlap;
use myloasm::unitig;
use std::fs::File;
use std::io::BufReader;
use std::io::BufWriter;
use std::path::Path;
use std::path::PathBuf;
use std::time::Instant;

fn main() {
    let total_start_time = Instant::now();
    let mut args = cli::Cli::parse();
    let (output_dir, temp_dir) = initialize_setup(&mut args);

    log::info!("Starting assembly...");

    // Step 1: Process k-mers, count k-mers, and get SNPmers
    let mut kmer_info = get_kmers_and_snpmers(&args, &output_dir);

    // Step 2: Get twin reads using k-mer information
    let twin_read_container = get_twin_reads_from_kmer_info(&mut kmer_info, &args, &output_dir);
    let twin_reads = &twin_read_container.twin_reads;

    // Step 3: Get overlaps between outer twin reads and construct raw unitig graph
    let overlaps = get_overlaps_from_twin_reads(&twin_read_container, &args, &output_dir);

    // Step 4: Construct raw unitig graph
    let unitig_graph =
        unitig::UnitigGraph::from_overlaps(&twin_read_container.twin_reads, &overlaps, &args);

    // light_progressive_cleaning(&mut unitig_graph, &twin_reads, &args, &temp_dir);

    // heavy_cleaning(&mut unitig_graph, &twin_reads, &args, &temp_dir);
    //     unitig_graph.to_gfa(
    //         output_dir.join("after_heavy.gfa"),
    //         true,
    //         true,
    //         &twin_reads,
    //         &args,
    // );

    let prog_dir = temp_dir.join("progressive");
    let mut contigs_per_coverage = vec![];

    let max_cov = unitig_graph
        .nodes
        .values()
        .map(|x| x.min_read_depth.unwrap())
        .max_by(|x, y| x.partial_cmp(y).unwrap())
        .unwrap_or(1.);

    for cov_thresh in 0..=max_cov as usize {
        log::info!("Processing coverage {}", cov_thresh);
        let prog_dir_cov = prog_dir.join(format!("cov_{}", cov_thresh));
        std::fs::create_dir_all(&prog_dir_cov).unwrap();
        let cov_thresh = cov_thresh as f64;
        let mut unitig_graph_with_threshold = unitig_graph.clone();
        unitig_graph_with_threshold.cut_coverage(cov_thresh);
        unitig_graph_with_threshold
            .get_sequence_info(&twin_reads, &types::GetSequenceInfoConfig::default());

        // Step 5: First round of cleaning. Progressive cleaning: tips, bubbles, bridged repeats
        light_progressive_cleaning(
            &mut unitig_graph_with_threshold,
            &twin_reads,
            &args,
            &prog_dir_cov,
            false
        );

        // Step 6: Second round of cleaning. TODO use coverage information, remove larger tips and bubbles. Imbue graph with sequence information
        heavy_cleaning(
            &mut unitig_graph_with_threshold,
            &twin_reads,
            &args,
            &prog_dir_cov,
        );

        // Push contigs at this coverage level
        unitig_graph_with_threshold.clear_edges();
        let graph_to_contigs: Vec<unitig::UnitigNode> =
            unitig_graph_with_threshold.nodes.into_values().collect();
        contigs_per_coverage.push(graph_to_contigs);
    }

    let mut unitig_graph =
        get_contigs_from_progressive_coverage(&contigs_per_coverage, &args, &temp_dir);

    // Step X: Progressive coverage filter
    //progressive_coverage_filter(&mut unitig_graph, &twin_reads, &args, &temp_dir);

    // Step 7: Align reads back to graph. TODO consensus and etc.
    log::info!("Beginning final alignment of reads to graph...");
    let mut get_seq_config = types::GetSequenceInfoConfig::default();
    get_seq_config.dna_seq_info = true;
    unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);

    let start = Instant::now();
    mapping::map_reads_to_unitigs(&mut unitig_graph, &kmer_info, &twin_reads, &args);

    log::info!(
        "Time elapsed for aligning reads to graph is {:?}",
        start.elapsed()
    );
    unitig_graph.to_gfa(
        output_dir.join("final_contig_graph.gfa"),
        true,
        true,
        &twin_reads,
        &args,
    );
    unitig_graph.to_gfa(
        output_dir.join("final_contig_graph_noseq.gfa"),
        true,
        false,
        &twin_reads,
        &args,
    );
    unitig_graph.to_fasta(output_dir.join("final_contigs.fa"), &args);

    // Step 8: TODO
    if !args.no_minimap2 {
        log::info!("Running minimap2...");
        consensus::write_to_paf(&unitig_graph, &twin_reads, &args);
        log::info!("Time elapsed for minimap2 is {:?}", start.elapsed());
    }

    unitig_graph.print_statistics();
    log::info!("Total time elapsed is {:?}", total_start_time.elapsed());
}

fn initialize_setup(args: &mut cli::Cli) -> (PathBuf, PathBuf) {
    // Initialize logger with CLI-specified level
    simple_logger::SimpleLogger::new()
        .with_level(args.log_level_filter())
        .init()
        .unwrap();

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
        args.snpmer_error_rate = 0.001;
        args.snpmer_threshold = 100.;
    }
    if args.r941 {
        args.snpmer_error_rate = 0.05;
        args.contain_subsample_rate = 20;
    }

    let output_dir = Path::new(args.output_dir.as_str());
    let temp_dir = output_dir.join("temp");

    if !output_dir.exists() {
        std::fs::create_dir_all(output_dir).unwrap();
        std::fs::create_dir_all(output_dir.join("temp")).unwrap();
    } else {
        if !output_dir.is_dir() {
            log::error!("Output directory is not a directory.");
            std::process::exit(1);
        }

        if !temp_dir.exists() {
            std::fs::create_dir_all(&temp_dir).unwrap();
        } else {
            if !temp_dir.is_dir() {
                log::error!("'temp' directory within output directory is not valid.");
                std::process::exit(1);
            }
        }
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
        // First: get twin reads from snpmers
        let twin_reads_raw = kmer_comp::twin_reads_from_snpmers(kmer_info, &args);
        let num_reads = twin_reads_raw.len();

        // Second: removed contained reads
        let outer_read_indices_raw =
            twin_graph::remove_contained_reads_twin(None, &twin_reads_raw, &args);

        // Third: map all reads to the outer (non-contained) reads
        log::info!("Mapping reads to non-contained (outer) reads...");
        let start = Instant::now();
        let outer_mapping_info =
            mapping::map_reads_to_outer_reads(&outer_read_indices_raw, &twin_reads_raw, &args);

        // Fourth: split chimeric twin reads based on mappings to outer reads
        let (split_twin_reads, split_outer_read_indices) =
            map_processing::split_outer_reads(twin_reads_raw, outer_mapping_info, &args);
        log::info!(
            "Gained {} reads after splitting chimeras and mapping to outer reads in {:?}",
            split_twin_reads.len() as i64 - num_reads as i64,
            start.elapsed()
        );

        // Fifth: the splitted chimeric reads may be contained within the original reads, so we remove contained reads again
        let outer_read_indices = twin_graph::remove_contained_reads_twin(
            Some(split_outer_read_indices),
            &split_twin_reads,
            &args,
        );

        twin_read_container = types::TwinReadContainer {
            twin_reads: split_twin_reads,
            outer_indices: outer_read_indices,
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
) -> Vec<TwinOverlap> {
    let empty_input = args.input_files.len() == 0;
    let twin_reads = &twin_read_container.twin_reads;
    let outer_read_indices = &twin_read_container.outer_indices;

    let overlaps;
    if empty_input && output_dir.join("overlaps.bin").exists() {
        overlaps = bincode::deserialize_from(BufReader::new(
            File::open(output_dir.join("overlaps.bin")).unwrap(),
        ))
        .unwrap();
    } else {
        log::info!("Getting overlaps between outer reads...");
        let start = Instant::now();
        overlaps =
            twin_graph::get_overlaps_outer_reads_twin(&twin_reads, &outer_read_indices, &args);
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
    output_temp: bool
) {
    let get_seq_config = types::GetSequenceInfoConfig::default();

    let mut iteration = 1;
    let divider = 3;

    loop {
        log::debug!("Cleaning graph iteration {}", iteration);
        let mut size_graph = unitig_graph.nodes.len();

        // Remove tips
        loop {
            let tip_length_cutoff = args.tip_length_cutoff;
            let read_cutoff = args.tip_read_cutoff;
            unitig_graph.remove_tips(tip_length_cutoff, read_cutoff, false);
            unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
            if output_temp{
                unitig_graph.to_gfa(
                    temp_dir.join(format!("{}-tip_unitig_graph.gfa", iteration)),
                    true,
                    false,
                    &twin_reads,
                    &args,
                );
            }
            if unitig_graph.nodes.len() == size_graph {
                break;
            }
            size_graph = unitig_graph.nodes.len();
        }

        // Pop bubbles
        let bubble_length_cutoff = args.small_bubble_threshold;
        unitig_graph.pop_bubbles(bubble_length_cutoff);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        if output_temp{
            unitig_graph.to_gfa(
                temp_dir.join(format!("{}-bubble_unitig_graph.gfa", iteration)),
                true,
                false,
                &twin_reads,
                &args,
            );
        }

        //Cut bridged repeats; "drop" cuts
        let prebridge_file = temp_dir.join(format!("{}-pre_bridge_cuts.txt", iteration));
        unitig_graph.resolve_bridged_repeats(
            &args,
            0.75 / (divider as f64) * iteration as f64,
            None,
            prebridge_file,
            FORWARD_READ_SAFE_SEARCH_CUTOFF,
            args.tip_read_cutoff,
            100_000,
        );
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        if output_temp{
            unitig_graph.to_gfa(
                temp_dir.join(format!("{}-resolve_unitig_graph.gfa", iteration)),
                true,
                false,
                &twin_reads,
                &args,
            );
        }
        
        iteration += 1;
        if iteration == 4 {
            break;
        }
    }

    let mut size_graph = unitig_graph.nodes.len();
    loop {
        unitig_graph.remove_tips(args.tip_length_cutoff, args.tip_read_cutoff, false);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.pop_bubbles(args.small_bubble_threshold);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        if unitig_graph.nodes.len() == size_graph {
            break;
        }
        size_graph = unitig_graph.nodes.len();
    }
    if output_temp{
        unitig_graph.to_gfa(
            temp_dir.join("pre-final_unitig_graph.gfa"),
            true,
            false,
            &twin_reads,
            &args,
        );
    }
}

fn get_contigs_from_progressive_coverage(
    contigs_per_coverage: &Vec<Vec<unitig::UnitigNode>>,
    _args: &cli::Cli,
    _temp_dir: &PathBuf,
) -> unitig::UnitigGraph {
    let mut final_unitig_graph = unitig::UnitigGraph::new();
    let mut used_reads = FxHashSet::default();
    for (cov, contigs) in contigs_per_coverage.iter().enumerate().rev() {
        for contig in contigs.iter() {
            if contig.min_read_depth.unwrap() as usize >= cov * 2 {
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

                let mut contig_with_new_id = contig.clone();
                contig_with_new_id.node_hash_id = final_unitig_graph.nodes.len();
                final_unitig_graph
                    .nodes
                    .insert(final_unitig_graph.nodes.len(), contig_with_new_id);

                for (read, _) in contig.read_indices_ori.iter() {
                    used_reads.insert(read);
                }
            }
        }
    }

    return final_unitig_graph;
}

fn heavy_cleaning(
    unitig_graph: &mut unitig::UnitigGraph,
    twin_reads: &Vec<types::TwinRead>,
    args: &cli::Cli,
    temp_dir: &PathBuf,
) {
    let get_seq_config = types::GetSequenceInfoConfig::default();
    let mut size_graph = unitig_graph.nodes.len();
    let mut counter = 0;
    let cov_score_thresholds = [10., 7., 5., 3., 2.];

    // TODO cut large tips (but keep them) and remove large bubbles
    loop {
        let ind = counter.min(cov_score_thresholds.len() - 1);
        let tip_length_cutoff_heavy = args.tip_length_cutoff * 5;
        let tip_read_cutoff_heavy = args.tip_read_cutoff;
        let bubble_threshold_heavy = args.max_bubble_threshold;
        let save_tips = true;

        remove_tips_until_stable(
            unitig_graph,
            twin_reads,
            tip_length_cutoff_heavy,
            tip_read_cutoff_heavy,
            bubble_threshold_heavy,
            temp_dir,
            save_tips,
            args,
        );
        unitig_graph.resolve_bridged_repeats(
            &args,
            0.50,
            Some(cov_score_thresholds[ind] as f64),
            temp_dir.join(format!("f{}-resolve_unitig_graph.txt", counter)),
            args.tip_length_cutoff * 5,
            args.tip_read_cutoff * 5,
            300_000,
        );
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
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
    _temp_dir: &PathBuf,
    save_tips: bool,
    _args: &cli::Cli,
) {
    let get_seq_config = types::GetSequenceInfoConfig::default();
    let mut size_graph = unitig_graph.nodes.len();
    loop {
        unitig_graph.remove_tips(tip_length_cutoff, tip_read_cutoff, save_tips);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.pop_bubbles(max_bubble_threshold);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        if unitig_graph.nodes.len() == size_graph {
            break;
        }
        size_graph = unitig_graph.nodes.len();
    }
}

fn _progressive_coverage_filter(
    unitig_graph: &mut unitig::UnitigGraph,
    twin_reads: &Vec<types::TwinRead>,
    args: &cli::Cli,
    temp_dir: &PathBuf,
) {
    let get_seq_config = types::GetSequenceInfoConfig::default();

    let tip_length_cutoff_heavy = args.tip_length_cutoff * 5;
    let tip_read_cutoff_heavy = args.tip_read_cutoff;
    let bubble_threshold_heavy = args.max_bubble_threshold;
    let save_tips = true;
    let mut counter = 0;

    loop {
        remove_tips_until_stable(
            unitig_graph,
            twin_reads,
            tip_length_cutoff_heavy,
            tip_read_cutoff_heavy,
            bubble_threshold_heavy,
            temp_dir,
            save_tips,
            args,
        );

        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        let num_cut = unitig_graph.progressive_cut_lowest();
        log::info!(
            "Progressive coverage filter iteration {}. Cut {} nodes",
            counter,
            num_cut
        );
        unitig_graph.get_sequence_info(twin_reads, &get_seq_config);

        let prog_folder = temp_dir.join("progressive");
        std::fs::create_dir_all(&prog_folder).unwrap();

        unitig_graph.to_gfa(
            //temp_dir.join("progressive-{}-_unitig_graph.gfa"),
            prog_folder.join(format!("progressive-{}-_unitig_graph.gfa", counter)),
            true,
            true,
            &twin_reads,
            &args,
        );
        if num_cut == 0 {
            break;
        }
        counter += 1;
    }
}
