use disjoint::DisjointSet;
use fxhash::FxHashSet;
use phaller::types::*;
use phaller::kmer_comp;
use phaller::graph;
use phaller::twin_graph;
use phaller::mapping;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::io::BufWriter;

fn main() {
    simple_logger::SimpleLogger::new()
        .with_level(log::LevelFilter::Info)
        .init()
        .unwrap();

    //    let bam_file = "/home/jshaw/software/phaller/98id-mn03.bam";
    //    let bam_file = "/home/jshaw/software/phaller/ill-mn03.bam";
    //    let bam_file = "/home/jshaw/software/phaller/98-mgh+mn.bam";
    //let bam_file = "/home/jshaw/software/phaller/temp/99.8-mgh-mn-mix.bam";
    //let bam_file = "/home/jshaw/software/phaller/temp/minimap-metaphlan.bam";

    //    let fastq_files = vec!["/home/jshaw/software/phaller/temp/98acc-10xcov-mn.fq", "/home/jshaw/software/phaller/temp/98acc-10xcov-mgh.fq"];
    //    let cuttlefish_file = "/home/jshaw/software/phaller/temp/cuttle.fa";
    //    let partig = "/home/jshaw/software/phaller/temp/ptig_kleb";

    let fastq_files = vec![
        "/home/jshaw/software/phaller/temp/98acc-10x1M-MN-03.fq",
//        "/home/jshaw/software/phaller/temp/98acc-10x1M-MN-03.fq",
        "/home/jshaw/software/phaller/temp/98acc-10x-1M-MGH.fq",
    ];
    let cuttlefish_file = "/home/jshaw/software/phaller/temp/cuttle-1M-kleb.fa";
    let partig = "/home/jshaw/software/phaller/temp/ptig_kleb-1M";

//    let fastq_files = vec!["/home/jshaw/software/phaller/temp/99.8-mgh-mn-mix.fq", "/home/jshaw/software/phaller/temp/99.8-mgh-mn-mix.fq"];
    //let cuttlefish_file = "/home/jshaw/software/phaller/temp/cuttle-hg002-99.8.fa";
    //let partig = "/home/jshaw/software/phaller/temp/ptig_chr21";

    let fastq_files = vec!["/home/jshaw/software/phaller/temp/99.5acc-10x-ALL-MGH.fq", "/home/jshaw/software/phaller/temp/99.5acc-10x-ALL-MN-03.fq"];
    //    let fastq_files = vec!["/home/jshaw/software/phaller/temp/98acc-10x1m-chr11.fq", "/home/jshaw/software/phaller/temp/98acc-10x1m-chr11-pat.fq"];
    //    let cuttlefish_file = "/home/jshaw/software/phaller/temp/cuttle-hg002-chr11.fa";
    //    let partig = "/home/jshaw/software/phaller/temp/ptig_chr11-1M";

    let (edges, num_elts, used_ids) = parse_tsv_file(partig);
    dbg!(edges.len(), num_elts, used_ids.len());
    let mut disjoint_set = DisjointSet::with_len(num_elts as usize);
    for (id1, id2) in edges {
        disjoint_set.join(id1 as usize, id2 as usize);
    }
    log::info!("Finished reading tsv file");

    //#let bam_file = env::args().nth(1).unwrap();

    //time this operation and print the time
    let threads = 5;
    let k = 27;
    let c = 10;
    let start = std::time::Instant::now();
    let (kmer_to_unitig, unitig_vec) = kmer_comp::parse_unitigs_into_table(cuttlefish_file);
    let mut all_reads = vec![];
    let mut all_maps = vec![];
    for fastq_file in fastq_files.iter() {
        let test_map = kmer_comp::read_to_split_kmers(fastq_file, k, threads);
        all_maps.push(test_map);
    }
            //    bam_parsing::get_snvs(bam_file);
    //bam_parsing::get_pileup(&bam_file, &fasta_file);
    let duration = start.elapsed();
    log::info!("Time elapsed in for parsing split-mers is: {:?}", duration);

    for fastq_file in fastq_files.iter() {
        let reads = kmer_comp::read_to_unitig_count_vector(fastq_file, &kmer_to_unitig, &unitig_vec);
        all_reads.push(reads);
    }

    let mut big_kmer_map = std::mem::take(&mut all_maps[0]);
    big_kmer_map.extend(std::mem::take(&mut all_maps[1]));

    let start = std::time::Instant::now();
    let snpmer_info = kmer_comp::get_snpmers(big_kmer_map, k);
    let duration = start.elapsed();
    log::info!("Time elapsed in for parsing snpmers is: {:?}", duration);

    let start = std::time::Instant::now();
    let twin_reads = kmer_comp::twin_reads_from_snpmers(&snpmer_info, &fastq_files, k, c, threads);
    let duration = start.elapsed();
    log::info!("Time elapsed for obtaining twin reads is: {:?}", duration);
    let all_reads_cat = std::mem::take(&mut all_reads[0])
        .into_iter()
        .chain(std::mem::take(&mut all_reads[1]).into_iter())
        .collect::<Vec<_>>();

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
    graph.transitive_reduction();
    let duration = start.elapsed();
    log::info!("Time elapsed for transitive reduction is: {:?}", duration);

    twin_graph::print_graph_stdout(&graph, "edges_reduced.tsv");


//    let outer_reads = graph::remove_contained_reads(&all_reads_cat, &disjoint_set, &used_ids);
//    log::info!("Finished removing contained reads");
//    let overlaps = graph::get_overlaps_outer_reads(&all_reads_cat, &disjoint_set, &used_ids, &outer_reads);
//    let mut graph = graph::read_graph_from_overlaps(all_reads_cat, &overlaps);
//    graph::print_graph_stdout(&graph, "edges.tsv");
//
//    graph.transitive_reduction();
//    graph::print_graph_stdout(&graph, "edges_reduced.tsv");
}




pub fn parse_tsv_file(file: &str) -> (Vec<(u64, u64)>, usize, FxHashSet<u64>) {
    let mut unique_ids = FxHashSet::default();
    let mut used_ids = FxHashSet::default();
    let mut records = vec![];
    // read file with Bufreader
    let file = File::open(file).expect("Could not open file");

    for line in BufReader::new(file).lines() {
        let record = line.unwrap();
        let record: Vec<&str> = record.split("\t").collect();
        if record[0].starts_with("C") {
            unique_ids.insert(record[1][1..].parse::<u64>().unwrap() - 1);
            continue;
        }
        //parse as u64
        let sim = record[7].parse::<f64>().unwrap();
        let num_kmers1 = record[4].parse::<u64>().unwrap();
        let num_kmers2 = record[5].parse::<u64>().unwrap();
        let num_kmers_shared = record[6].parse::<u64>().unwrap();
        if num_kmers_shared == num_kmers1 || num_kmers_shared == num_kmers2 {
            //    continue
        }
        let id1 = record[1][1..].parse::<u64>().unwrap();
        let id2 = record[2][1..].parse::<u64>().unwrap();
        used_ids.insert(id1 - 1);
        used_ids.insert(id2 - 1);
        records.push((id1 - 1, id2 - 1));
    }
    (records, unique_ids.len(), used_ids)
}



