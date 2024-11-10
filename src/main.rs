use phaller::unitig;
use std::io::BufReader;
use std::io::BufWriter;
use bincode;
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
    let overlaps;
    let twin_reads; 

    if fastq_files[0].contains("overlaps.bin") && fastq_files[1].contains("twin_reads.bin"){
        overlaps = bincode::deserialize_from(BufReader::new(std::fs::File::open(&fastq_files[0]).unwrap())).unwrap();
        twin_reads = bincode::deserialize_from(BufReader::new(
            std::fs::File::open(&fastq_files[1]).unwrap())).unwrap();
        }
    else{
        for fastq_file in fastq_files.iter() {
            let test_map = kmer_comp::read_to_split_kmers(fastq_file, k, threads, hpc);
            all_maps.push(test_map);
        }
        log::info!("Time elapsed in for parsing split-mers is: {:?}", start.elapsed());

        let mut big_kmer_map = std::mem::take(&mut all_maps[0]);
        for i in 1..all_maps.len() {
            big_kmer_map.extend(std::mem::take(&mut all_maps[i]));
        }

        let start = std::time::Instant::now();
        let snpmer_info = kmer_comp::get_snpmers(big_kmer_map, k);
        log::info!("Time elapsed in for parsing snpmers is: {:?}", start.elapsed());

        let start = std::time::Instant::now();
        twin_reads = kmer_comp::twin_reads_from_snpmers(&snpmer_info, &fastq_files, k, c, threads,hpc);
        log::info!("Time elapsed for obtaining twin reads is: {:?}", start.elapsed());

        let start = std::time::Instant::now();
        let outer_read_indices = twin_graph::remove_contained_reads_twin(&twin_reads, k.try_into().unwrap(), c.try_into().unwrap());
        log::info!("Time elapsed for removing contained reads is: {:?}", start.elapsed());

        let start = std::time::Instant::now();
        overlaps = twin_graph::get_overlaps_outer_reads_twin(&twin_reads, &outer_read_indices, k.try_into().unwrap(), c.try_into().unwrap());
        log::info!("Time elapsed for getting overlaps is: {:?}", start.elapsed());

        bincode::serialize_into(BufWriter::new(std::fs::File::create("overlaps.bin").unwrap()), &overlaps).unwrap();
        bincode::serialize_into(BufWriter::new(std::fs::File::create("twin_reads.bin").unwrap()), &twin_reads).unwrap();

    }

    let start = std::time::Instant::now();
    let mut graph = twin_graph::read_graph_from_overlaps_twin(&twin_reads, &overlaps, c.try_into().unwrap());
    log::info!("Time elapsed for obtaining overlap graphs is: {:?}", start.elapsed());

    let start = std::time::Instant::now();
    twin_graph::print_graph_stdout(&graph, "edges.tsv");
    graph.transitive_reduction();
    log::info!("Time elapsed for transitive reduction is: {:?}", start.elapsed());

    twin_graph::print_graph_stdout(&graph, "edges_reduced.tsv");

    let start = std::time::Instant::now();
    let mut unitig_graph = unitig::UnitigGraph::from_overlap_graph(&graph, &twin_reads);
    if cfg!(debug_assertions) {
        unitig_graph.test_consistent_left_right_edges();
    }
    unitig_graph.cut_overlap_boundaries();
    log::info!("Time elapsed for initial unitig construction is: {:?}", start.elapsed());
    unitig_graph.to_gfa("unitig_graph.gfa", true, &twin_reads);

    let start = std::time::Instant::now();
    unitig_graph.remove_tips(args.tip_length_cutoff, args.tip_read_cutoff);
    unitig_graph.re_unitig();
    unitig_graph.cut_overlap_boundaries();
    unitig_graph.to_gfa("tip_unitig_graph.gfa", true, &twin_reads);
    log::info!("Time elapsed for unitig construction is: {:?}", start.elapsed());
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


