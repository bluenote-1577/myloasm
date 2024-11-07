use phaller::unitig;
use dashmap::DashMap;
use clap::Parser;
use phaller::cli;
use phaller::kmer_comp;
use phaller::twin_graph;

fn main() {
    let args = cli::Cli::parse();

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

    let threads = args.threads;
   // Initialize thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();


    let fastq_files = args.input_files;
    let k = args.kmer_size;
    let c = args.c;
    let start = std::time::Instant::now();
    let hpc = false;
    let mut all_maps = vec![];
    for fastq_file in fastq_files.iter() {
        let test_map = kmer_comp::read_to_split_kmers(fastq_file, k, threads, hpc);
        all_maps.push(test_map);
    }
    let duration = start.elapsed();
    log::info!("Time elapsed in for parsing split-mers is: {:?}", duration);

    let mut big_kmer_map = std::mem::take(&mut all_maps[0]);
    for i in 1..all_maps.len() {
        big_kmer_map.extend(std::mem::take(&mut all_maps[i]));
    }

    let start = std::time::Instant::now();
    let snpmer_info = kmer_comp::get_snpmers(big_kmer_map, k);
    let duration = start.elapsed();
    log::info!("Time elapsed in for parsing snpmers is: {:?}", duration);

    let start = std::time::Instant::now();
    let twin_reads = kmer_comp::twin_reads_from_snpmers(&snpmer_info, &fastq_files, k, c, threads,hpc);
    let duration = start.elapsed();
    log::info!("Time elapsed for obtaining twin reads is: {:?}", duration);

    let start = std::time::Instant::now();
    let outer_read_indices = twin_graph::remove_contained_reads_twin(&twin_reads, k.try_into().unwrap(), c.try_into().unwrap());
    let duration = start.elapsed();

    log::info!("Time elapsed for removing contained reads is: {:?}", duration);

    let start = std::time::Instant::now();
    let overlaps = twin_graph::get_overlaps_outer_reads_twin(&twin_reads, &outer_read_indices, k.try_into().unwrap(), c.try_into().unwrap());
    let duration = start.elapsed();
    log::info!("Time elapsed for getting overlaps is: {:?}", duration);

    let start = std::time::Instant::now();
    let mut graph = twin_graph::read_graph_from_overlaps_twin(twin_reads, &overlaps, c.try_into().unwrap());
    let duration = start.elapsed();

    log::info!("Time elapsed for obtaining overlap graphs is: {:?}", duration);

    let start = std::time::Instant::now();
    twin_graph::print_graph_stdout(&graph, "edges.tsv");
    graph.transitive_reduction();
    let duration = start.elapsed();
    log::info!("Time elapsed for transitive reduction is: {:?}", duration);

    twin_graph::print_graph_stdout(&graph, "edges_reduced.tsv");

    let mut unitig_graph = unitig::UnitigGraph::from_overlap_graph(&graph);
    if cfg!(debug_assertions) {
        unitig_graph.test_consistent_left_right_edges();
    }
    unitig_graph.cut_overlap_boundaries();
    unitig_graph.to_gfa("unitig_graph.gfa", true, &graph);

    //graph.remove_tips();
    //twin_graph::print_graph_stdout(&graph, "tips_removed.tsv");
    //graph.remove_bubbles();
    //twin_graph::print_graph_stdout(&graph, "bubbles_removed.tsv");


//    let outer_reads = graph::remove_contained_reads(&all_reads_cat, &disjoint_set, &used_ids);
//    log::info!("Finished removing contained reads");
//    let overlaps = graph::get_overlaps_outer_reads(&all_reads_cat, &disjoint_set, &used_ids, &outer_reads);
//    let mut graph = graph::read_graph_from_overlaps(all_reads_cat, &overlaps);
//    graph::print_graph_stdout(&graph, "edges.tsv");
//
//    graph.transitive_reduction();
//    graph::print_graph_stdout(&graph, "edges_reduced.tsv");
}


