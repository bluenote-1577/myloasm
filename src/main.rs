use myloasm::constants::FORWARD_READ_SAFE_SEARCH_CUTOFF;
use myloasm::seq_parse;
use myloasm::types::TwinOverlap;
use myloasm::map_processing;
use myloasm::unitig;
use myloasm::consensus;
use std::path::Path;
use std::path::PathBuf;
use std::fs::File;
use std::time::Instant;
use myloasm::mapping;
use myloasm::types;
use std::io::BufReader;
use std::io::BufWriter;
use bincode;
use clap::Parser;
use myloasm::cli;
use myloasm::kmer_comp;
use myloasm::twin_graph;

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
    let mut unitig_graph = unitig::UnitigGraph::from_overlaps(&twin_read_container.twin_reads, &overlaps, &args);

    // Step 5: First round of cleaning. Progressive cleaning: tips, bubbles, bridged repeats
    light_progressive_cleaning(&mut unitig_graph, &twin_reads, &args, &temp_dir);
    
    // Step 6: Second round of cleaning. TODO use coverage information, remove larger tips and bubbles. Imbue graph with sequence information
    heavy_cleaning(&mut unitig_graph, &twin_reads, &args, &temp_dir);

    // Step 7: Align reads back to graph. TODO consensus and etc. 
    log::info!("Aligning reads back to graph...");
    let start = Instant::now();
    mapping::map_reads_to_unitigs(&mut unitig_graph, &kmer_info, &twin_reads, &args);
    log::info!("Time elapsed for aligning reads to graph is {:?}", start.elapsed());
    unitig_graph.to_gfa(output_dir.join("final_contig_graph.gfa"), true, true, &twin_reads, &args);
    unitig_graph.to_gfa(output_dir.join("final_contig_graph_noseq.gfa"), true, false, &twin_reads, &args);
    unitig_graph.to_fasta(output_dir.join("final_contigs.fa"), &args);

    // Step 8: TODO
    if !args.no_minimap2{
        log::info!("Running minimap2...");
        consensus::write_to_paf(&unitig_graph, &twin_reads, &args);
        log::info!("Time elapsed for minimap2 is {:?}", start.elapsed());
    }

    unitig_graph.print_statistics();
    log::info!("Total time elapsed is {:?}", total_start_time.elapsed());
}

fn initialize_setup(args: &mut cli::Cli) -> (PathBuf, PathBuf){

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

    if args.hifi{
        args.snpmer_error_rate = 0.001;
        args.snpmer_threshold = 100.;
    }
    if args.r941{
        args.snpmer_error_rate = 0.05;
        args.contain_subsample_rate = 20;
    }

    let output_dir = Path::new(args.output_dir.as_str());
    let temp_dir = output_dir.join("temp");

    if !output_dir.exists(){
        std::fs::create_dir_all(output_dir).unwrap();
        std::fs::create_dir_all(output_dir.join("temp")).unwrap();
    }
    else{
        if !output_dir.is_dir(){
            log::error!("Output directory is not a directory.");
            std::process::exit(1);
        }

        if !temp_dir.exists(){
            std::fs::create_dir_all(&temp_dir).unwrap();
        }
        else{
            if !temp_dir.is_dir(){
                log::error!("'temp' directory within output directory is not valid.");
                std::process::exit(1);
            }
        }
    }

    return (output_dir.to_path_buf(), temp_dir.to_path_buf());
}

fn get_kmers_and_snpmers(
    args: &cli::Cli, 
    output_dir: &PathBuf, 
) -> types::KmerGlobalInfo{

    let empty_input = args.input_files.len() == 0;

    let kmer_info;
    if empty_input {
        if !output_dir.join("snpmer_info.bin").exists(){
            log::error!("No input files provided. See --help for usage.");
            std::process::exit(1);
        }
    }

    if empty_input && output_dir.join("snpmer_info.bin").exists(){
        kmer_info = bincode::deserialize_from(BufReader::new(File::open(output_dir.join("snpmer_info.bin")).unwrap())).unwrap();
        log::info!("Loaded snpmer info from file.");
    }
    else{
        let start = Instant::now();
        let big_kmer_map = seq_parse::read_to_split_kmers(args.kmer_size, args.threads, &args);
        log::info!("Time elapsed in for counting k-mers is: {:?}", start.elapsed());

        let start = Instant::now();
        kmer_info = kmer_comp::get_snpmers(big_kmer_map, args.kmer_size, &args);
        log::info!("Time elapsed in for parsing snpmers is: {:?}", start.elapsed());
        bincode::serialize_into(BufWriter::new(File::create(output_dir.join("snpmer_info.bin")).unwrap()), &kmer_info).unwrap();
    }
    return kmer_info;
}

fn get_twin_reads_from_kmer_info(
    kmer_info: &mut types::KmerGlobalInfo, 
    args: &cli::Cli, 
    output_dir: &PathBuf, 
) -> types::TwinReadContainer
{
    let empty_input = args.input_files.len() == 0;

    let twin_read_container;
    if empty_input && output_dir.join("twin_reads.bin").exists(){    
        twin_read_container = bincode::deserialize_from(BufReader::new(File::open(output_dir.join("twin_reads.bin")).unwrap())).unwrap();
        log::info!("Loaded twin reads from file.");
    }
    else{
        // First: get twin reads from snpmers
        let twin_reads_raw = kmer_comp::twin_reads_from_snpmers(kmer_info, &args);
        let num_reads = twin_reads_raw.len();

        // Second: removed contained reads 
        let outer_read_indices_raw = twin_graph::remove_contained_reads_twin(None, &twin_reads_raw, &args);

        // Third: map all reads to the outer (non-contained) reads
        log::info!("Mapping reads to non-contained (outer) reads...");
        let start = Instant::now();
        let outer_mapping_info = mapping::map_reads_to_outer_reads(&outer_read_indices_raw, &twin_reads_raw, &args);

        // Fourth: split chimeric twin reads based on mappings to outer reads
        let (split_twin_reads, split_outer_read_indices) = map_processing::split_outer_reads(twin_reads_raw, outer_mapping_info, &args);
        log::info!("Gained {} reads after splitting chimeras and mapping to outer reads in {:?}", split_twin_reads.len() as i64 - num_reads as i64, start.elapsed());

        // Fifth: the splitted chimeric reads may be contained within the original reads, so we remove contained reads again
        let outer_read_indices = twin_graph::remove_contained_reads_twin(Some(split_outer_read_indices), &split_twin_reads, &args);

        twin_read_container = types::TwinReadContainer{
            twin_reads: split_twin_reads,
            outer_indices: outer_read_indices,
            tig_reads: vec![],
        };

        bincode::serialize_into(BufWriter::new(File::create(output_dir.join("twin_reads.bin")).unwrap()), &twin_read_container).unwrap();
    }
    return twin_read_container;
}

fn get_overlaps_from_twin_reads(
    twin_read_container: &types::TwinReadContainer,
    args: &cli::Cli,
    output_dir: &PathBuf,
) -> Vec<TwinOverlap>{
    let empty_input = args.input_files.len() == 0;
    let twin_reads = &twin_read_container.twin_reads;
    let outer_read_indices = &twin_read_container.outer_indices;

    let overlaps;
    if empty_input && output_dir.join("overlaps.bin").exists(){
        overlaps = bincode::deserialize_from(BufReader::new(File::open(output_dir.join("overlaps.bin")).unwrap())).unwrap();
    }
    else{
        log::info!("Getting overlaps between outer reads...");
        let start = Instant::now();
        overlaps = twin_graph::get_overlaps_outer_reads_twin(&twin_reads, &outer_read_indices, &args);
        log::info!("Time elapsed for getting overlaps is: {:?}", start.elapsed());
        bincode::serialize_into(BufWriter::new(File::create(output_dir.join("overlaps.bin")).unwrap()), &overlaps).unwrap();
    }
    return overlaps;

}

fn light_progressive_cleaning(
    unitig_graph: &mut unitig::UnitigGraph,
    twin_reads: &Vec<types::TwinRead>,
    args: &cli::Cli,
    temp_dir: &PathBuf,
){
    let get_seq_config = types::GetSequenceInfoConfig::default();

    let mut iteration = 1;
    let divider = 3;

    loop{
        log::info!("Cleaning graph iteration {}", iteration);
        let mut size_graph = unitig_graph.nodes.len();

        // Remove tips
        loop{
            let tip_length_cutoff = args.tip_length_cutoff;
            let read_cutoff = args.tip_read_cutoff;
            unitig_graph.remove_tips(tip_length_cutoff, read_cutoff, false);
            unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
            unitig_graph.to_gfa(temp_dir.join(format!("{}-tip_unitig_graph.gfa", iteration)), true, false, &twin_reads, &args);
            if unitig_graph.nodes.len() == size_graph{
                break;
            }
            size_graph = unitig_graph.nodes.len();
        }

        // Pop bubbles
        let bubble_length_cutoff = args.max_bubble_threshold/2; 
        unitig_graph.pop_bubbles(bubble_length_cutoff);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.to_gfa(temp_dir.join(format!("{}-bubble_unitig_graph.gfa", iteration)), true, false, &twin_reads, &args);

        //Cut bridged repeats; "drop" cuts
        let prebridge_file = temp_dir.join(format!("{}-pre_bridge_cuts.txt", iteration));
        unitig_graph.resolve_bridged_repeats(&args, 0.75/(divider as f64) * iteration as f64, prebridge_file, FORWARD_READ_SAFE_SEARCH_CUTOFF, args.tip_read_cutoff, 100_000);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.to_gfa(temp_dir.join(format!("{}-resolve_unitig_graph.gfa", iteration)), true, false, &twin_reads, &args);
        
        iteration += 1;
        if iteration == 4{
            break;
        }
    }

    let mut size_graph = unitig_graph.nodes.len();
    loop{
        unitig_graph.remove_tips(args.tip_length_cutoff, args.tip_read_cutoff, false);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.pop_bubbles(args.max_bubble_threshold/2);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        if unitig_graph.nodes.len() == size_graph{
            break;
        }
        size_graph = unitig_graph.nodes.len();
    }
    unitig_graph.to_gfa(temp_dir.join("pre-final_unitig_graph.gfa"), true, false, &twin_reads, &args);
}

fn heavy_cleaning(
    unitig_graph: &mut unitig::UnitigGraph,
    twin_reads: &Vec<types::TwinRead>,
    args: &cli::Cli,
    temp_dir: &PathBuf,
)
{
    let mut get_seq_config = types::GetSequenceInfoConfig::default();
    let mut size_graph = unitig_graph.nodes.len();
    let mut counter = 0;

    // TODO cut large tips (but keep them) and remove large bubbles
    loop{
        unitig_graph.remove_tips(args.tip_length_cutoff * 5, args.tip_read_cutoff, true);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.pop_bubbles(args.max_bubble_threshold);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        unitig_graph.resolve_bridged_repeats(&args, 0.50, temp_dir.join(format!("f{}-resolve_unitig_graph.txt", counter)), args.tip_length_cutoff * 5, args.tip_read_cutoff * 5, 300_000);
        unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
        if unitig_graph.nodes.len() == size_graph{
            break;
        }
        size_graph = unitig_graph.nodes.len();
        counter += 1;
    }

    get_seq_config.dna_seq_info = true;
    unitig_graph.get_sequence_info(&twin_reads, &get_seq_config);
}