use phaller::unitig;
use std::path::Path;
use phaller::mapping;
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


    let fastq_files = args.input_files.clone();
    let k = args.kmer_size;
    let start = std::time::Instant::now();
    let hpc = false;
    let overlaps;
    let mut twin_reads; 
    let kmer_info;
    let output_dir = Path::new(args.output_dir.as_str());
    if !output_dir.exists(){
        std::fs::create_dir_all(output_dir).unwrap();
    }
    else{
        if !output_dir.is_dir(){
            log::error!("Output directory is not a directory.");
            std::process::exit(1);
        }
    }

    if fastq_files.len() == 0{
        if !output_dir.join("overlaps.bin").exists(){
            log::error!("No input files provided. See --help for usage.");
            std::process::exit(1);
        }
        overlaps = bincode::deserialize_from(BufReader::new(std::fs::File::open(output_dir.join("overlaps.bin")).unwrap())).unwrap();
        twin_reads = bincode::deserialize_from(BufReader::new(std::fs::File::open(output_dir.join("twin_reads.bin")).unwrap())).unwrap();
        kmer_info = bincode::deserialize_from(BufReader::new(std::fs::File::open(output_dir.join("snpmer_info.bin")).unwrap())).unwrap();
    }
    else{
        log::info!("Starting assembly...");
        let big_kmer_map = kmer_comp::read_to_split_kmers(k, threads, &args);
        log::info!("Time elapsed in for counting k-mers is: {:?}", start.elapsed());

        let start = std::time::Instant::now();
        kmer_info = kmer_comp::get_snpmers(big_kmer_map, k, &args);
        log::info!("Time elapsed in for parsing snpmers is: {:?}", start.elapsed());

        let start = std::time::Instant::now();
        twin_reads = kmer_comp::twin_reads_from_snpmers(&kmer_info, &fastq_files, &args);
        log::info!("Time elapsed for obtaining twin reads is: {:?}", start.elapsed());

        // let start = std::time::Instant::now();
        // let num_reads = twin_reads.len();
        // twin_reads = kmer_comp::kmer_chimera_filter(twin_reads, &big_kmer_map, args.kmer_size as u64);
        // log::info!("Gained # of reads from splitting chimeras: {} in {:?}", twin_reads.len() - num_reads, start.elapsed());

        let start = std::time::Instant::now();
        let mut outer_read_indices = twin_graph::remove_contained_reads_twin(None, &twin_reads, &args);
        log::info!("Time elapsed for removing contained reads is: {:?}", start.elapsed());

        let start = std::time::Instant::now();
        let num_reads = twin_reads.len();
        let outer_mapping_info = mapping::map_reads_to_outer_reads(&outer_read_indices, &twin_reads, &args);
        let tup = kmer_comp::split_outer_reads(twin_reads, outer_mapping_info, &args);
        twin_reads = tup.0;
        outer_read_indices = tup.1;
        outer_read_indices = twin_graph::remove_contained_reads_twin(Some(outer_read_indices), &twin_reads, &args);
        log::info!("Gained {} reads after splitting chimeras and mapping to outer reads in {:?}", twin_reads.len() - num_reads, start.elapsed());

        let start = std::time::Instant::now();
        overlaps = twin_graph::get_overlaps_outer_reads_twin(&twin_reads, &outer_read_indices, &args);
        log::info!("Time elapsed for getting overlaps is: {:?}", start.elapsed());

        bincode::serialize_into(BufWriter::new(std::fs::File::create(output_dir.join("overlaps.bin")).unwrap()), &overlaps).unwrap();
        bincode::serialize_into(BufWriter::new(std::fs::File::create(output_dir.join("twin_reads.bin")).unwrap()), &twin_reads).unwrap();
        bincode::serialize_into(BufWriter::new(std::fs::File::create(output_dir.join("snpmer_info.bin")).unwrap()), &kmer_info).unwrap();
    }
    let twin_reads = twin_reads;
    
    let start = std::time::Instant::now();
    let mut graph = twin_graph::read_graph_from_overlaps_twin(&twin_reads, &overlaps, &args);
    graph.transitive_reduction();
    twin_graph::print_graph_stdout(&graph, output_dir.join("edges_reduced.tsv"));
    let mut unitig_graph = unitig::UnitigGraph::from_overlap_graph(&graph, &twin_reads);
    log::info!("Time elapsed for obtaining overlap graphs is: {:?}", start.elapsed());

    if cfg!(debug_assertions) {
        unitig_graph.test_consistent_left_right_edges();
    }
    unitig_graph.get_sequence_info(&twin_reads, true, false);
    log::info!("Time elapsed for initial unitig construction is: {:?}", start.elapsed());
    unitig_graph.to_gfa(output_dir.join("unitig_graph.gfa"), true, &twin_reads);


    let start = std::time::Instant::now();
    let mut iteration = 1;
    let divider = 3;

    loop{
        log::info!("Cleaning graph iteration {}", iteration);
        unitig_graph.resolve_bridged_repeats(&args, 0.75/(divider as f64) * iteration as f64);
        unitig_graph.get_sequence_info(&twin_reads, true, false);
        unitig_graph.to_gfa(output_dir.join(format!("{}-resolve_unitig_graph.gfa", iteration)), true, &twin_reads);

        let mut size_graph = unitig_graph.nodes.len();
        loop{
            unitig_graph.remove_tips(args.tip_length_cutoff/divider, args.tip_read_cutoff/divider * iteration);
            unitig_graph.get_sequence_info(&twin_reads, true, false);
            unitig_graph.to_gfa(output_dir.join(format!("{}-tip_unitig_graph.gfa", iteration)), true, &twin_reads);
            if unitig_graph.nodes.len() == size_graph{
                break;
            }
            size_graph = unitig_graph.nodes.len();
        }

        unitig_graph.pop_bubbles(75000/divider * iteration);
        unitig_graph.get_sequence_info(&twin_reads, true, false);
        unitig_graph.to_gfa(output_dir.join(format!("{}-bubble_unitig_graph.gfa", iteration)), true, &twin_reads);

        iteration += 1;
        if iteration == 4{
            break;
        }
    }

    let mut size_graph = unitig_graph.nodes.len();
    loop{
        unitig_graph.remove_tips(args.tip_length_cutoff, args.tip_read_cutoff);
        unitig_graph.get_sequence_info(&twin_reads, true, false);
        unitig_graph.pop_bubbles(75000);
        unitig_graph.get_sequence_info(&twin_reads, true, false);
        if unitig_graph.nodes.len() == size_graph{
            break;
        }
        size_graph = unitig_graph.nodes.len();
    }
    unitig_graph.to_gfa(output_dir.join("pre-final_unitig_graph.gfa"), true, &twin_reads);

    // Cut Z-edges 
    unitig_graph.cut_z_edges(&args);
    unitig_graph.get_sequence_info(&twin_reads, true, true);
    log::info!("Time elapsed for cleaning graph is {:?}", start.elapsed());

    log::info!("Aligning reads back to graph...");
    let start = std::time::Instant::now();
    mapping::map_reads_to_unitigs(&mut unitig_graph, &kmer_info, &twin_reads, &args);
    log::info!("Time elapsed for aligning reads to graph is {:?}", start.elapsed());
    unitig_graph.to_gfa(output_dir.join("final_unitig_graph.gfa"), true, &twin_reads);
    unitig_graph.print_statistics();
}


