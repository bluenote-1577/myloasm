use csv;
use disjoint::DisjointSet;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use needletail::*;
use phaller::bam_parsing;
use phaller::seeding;
use phaller::types::*;
use rust_htslib::bam;
use rust_htslib::bam::{Read, Reader, Record};
use std::env;
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

    //let fastq_files = vec!["/home/jshaw/software/phaller/temp/99.8acc-10x100kb-mat.fq", "/home/jshaw/software/phaller/temp/99.8acc-10x100kb-pat.fq"];
    //let cuttlefish_file = "/home/jshaw/software/phaller/temp/cuttle-hg002-99.8.fa";
    //let partig = "/home/jshaw/software/phaller/temp/ptig_chr21";

    //let fastq_files = vec!["/home/jshaw/software/phaller/temp/99acc-10x100kb-mat.fq", "/home/jshaw/software/phaller/temp/99acc-10x100kb-pat.fq"];
    //    let fastq_files = vec!["/home/jshaw/software/phaller/temp/98acc-10x1m-chr11.fq", "/home/jshaw/software/phaller/temp/98acc-10x1m-chr11-pat.fq"];
    //    let cuttlefish_file = "/home/jshaw/software/phaller/temp/cuttle-hg002-chr11.fa";
    //    let partig = "/home/jshaw/software/phaller/temp/ptig_chr11-1M";
    //
    let (edges, num_elts, used_ids) = parse_tsv_file(partig);
    dbg!(edges.len(), num_elts, used_ids.len());
    let mut disjoint_set = DisjointSet::with_len(num_elts as usize);
    for (id1, id2) in edges {
        disjoint_set.join(id1 as usize, id2 as usize);
    }
    log::info!("Finished reading tsv file");

    //#let bam_file = env::args().nth(1).unwrap();

    //time this operation and print the time
    let start = std::time::Instant::now();
    let (kmer_to_unitig, unitig_vec) = parse_unitigs_into_table(cuttlefish_file);
    let mut all_reads = vec![];
    for fastq_file in fastq_files {
        let reads = read_to_unitig_count_vector(fastq_file, &kmer_to_unitig, &unitig_vec);
        all_reads.push(reads);
    }
    //    bam_parsing::get_snvs(bam_file);
    //bam_parsing::get_pileup(&bam_file, &fasta_file);
    let duration = start.elapsed();
    log::info!("Time elapsed in parse_bam() is: {:?}", duration);
    let all_reads_cat = std::mem::take(&mut all_reads[0])
        .into_iter()
        .chain(std::mem::take(&mut all_reads[1]).into_iter())
        .collect::<Vec<_>>();


    let outer_reads = remove_contained_reads(&all_reads_cat, &disjoint_set, &used_ids);
    log::info!("Finished removing contained reads");
    let overlaps = get_overlaps_outer_reads(&all_reads_cat, &disjoint_set, &used_ids, &outer_reads);
    let mut graph = read_graph_from_overlaps(all_reads_cat, &overlaps);
    print_graph_stdout(&graph, "edges.tsv");

    graph.transitive_reduction();
    print_graph_stdout(&graph, "edges_reduced.tsv");
}

fn print_graph_stdout(graph: &OverlapGraph, file: &str) {
    let mut bufwriter = BufWriter::new(File::create(file).unwrap());
    let all_reads_cat = &graph.reads;
    for edge in graph.edges.iter() {
        if let Some(edge) = edge {
            let i = edge.node1;
            let j = edge.node2;
            let read = &all_reads_cat[i];
            let read2 = &all_reads_cat[j];
            let frac_shared_variable = edge.variable_tigs as f64 / edge.variable_roots as f64;
            let forward1 = edge.forward1;
            let forward2 = edge.forward2;
            let aln_len = edge.overlap_len_tigs;

            if read.id.contains("junk") || read2.id.contains("junk") || read.id.contains("chimera") || read2.id.contains("chimera"){
                continue;
            }
            let name1 = read.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[0];
            let name2 = read2.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[0];
            let range1 = read.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[2];
            let range2 = read2.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[2];

            writeln!(bufwriter, 
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                i, j, name1, name2, range1, range2, forward1, forward2
            );
        }
    }
}

fn parse_unitigs_into_table(cuttlefish_file: &str) -> (FxHashMap<u64, u32>, Vec<Vec<u8>>) {
    let mut kmer_to_unitig_count: FxHashMap<u64, u32> = fxhash::FxHashMap::default();
    let mut reader = needletail::parse_fastx_file(cuttlefish_file).expect("valid path");
    let mut count = 0;
    let mut unitig_vec = vec![];
    while let Some(record) = reader.next() {
        let rec = record.expect("Error reading record");
        let seq = rec.seq();
        let mut kmers = vec![];
        seeding::fmh_seeds(&seq, &mut kmers, 10, 27);
        if kmers.len() > 0 {
            for kmer in kmers {
                kmer_to_unitig_count.entry(kmer).or_insert(count);
            }
            unitig_vec.push(seq.to_vec());
        }
        count += 1;
    }
    return (kmer_to_unitig_count, unitig_vec);
}

fn read_to_unitig_count_vector(
    fastq_file: &str,
    kmer_to_unitig_count: &FxHashMap<u64, u32>,
    unitig_vec: &Vec<Vec<u8>>,
) -> Vec<TigRead> {
    let mut reads = vec![];
    let mut reader = needletail::parse_fastx_file(fastq_file).expect("valid path");
    while let Some(record) = reader.next() {
        let rec = record.expect("Error reading record");
        let seq = rec.seq();
        let mut kmers = vec![];
        seeding::fmh_seeds(&seq, &mut kmers, 10, 27);
        if kmers.len() > 0 {
            let mut last_count = u32::MAX;
            let mut unitig_string = vec![];
            for kmer in kmers {
                if let Some(unitig_id) = kmer_to_unitig_count.get(&kmer) {
                    if *unitig_id != last_count {
                        unitig_string.push(*unitig_id);
                        last_count = *unitig_id;
                    }
                }
            }
            reads.push(TigRead{
                tig_seq: unitig_string,
                id: String::from_utf8_lossy(rec.id()).to_string(),
            });
        }
    }
    return reads;
}

fn edit_distance(v1: &[u32], v2: &[u32]) -> usize {
    let len1 = v1.len();
    let len2 = v2.len();

    // Create a 2D vector to store edit distances
    let mut dp = vec![vec![0; len2 + 1]; len1 + 1];

    // Initialize the first row and column
    for i in 0..=len1 {
        dp[i][0] = i;
    }
    for j in 0..=len2 {
        dp[0][j] = j;
    }

    // Fill the rest of the dp table
    for i in 1..=len1 {
        for j in 1..=len2 {
            if v1[i - 1] == v2[j - 1] {
                dp[i][j] = dp[i - 1][j - 1];
            } else {
                dp[i][j] = 1 + dp[i - 1][j - 1].min(dp[i][j - 1]).min(dp[i - 1][j]);
            }
        }
    }

    // The edit distance is the value in the bottom-right corner of the dp table
    dp[len1][len2]
}

fn disjoint_distance(
    v1: &[u32],
    v2: &[u32],
    disjoint_set: &DisjointSet,
    used_ids: &FxHashSet<u64>,
    ind1: usize,
    ind2: usize,
) -> TigdexOverlap {
    let v1_as_disjoint = v1
        .iter()
        .map(|&x| (disjoint_set.root_of(x as usize), x))
        .collect::<Vec<(usize, u32)>>();
    let v2_as_disjoint = v2
        .iter()
        .map(|&x| (disjoint_set.root_of(x as usize), x))
        .collect::<Vec<(usize, u32)>>();
    let v1_roots = v1_as_disjoint.iter().map(|x| x.0).collect::<Vec<usize>>();
    let v2_roots = v2_as_disjoint.iter().map(|x| x.0).collect::<Vec<usize>>();

    let chain_info = find_optimal_chain(&v1_roots, &v2_roots);
    let chain = &chain_info.chain;

    let mut shared_x = 0;
    let mut variable_roots = 0;
    let mut variable_shared_x = 0;

    for (start1, start2, _) in chain.iter() {
        if used_ids.contains(&(v1_as_disjoint[*start1].1 as u64)) {
            variable_roots += 1;
            if v1_as_disjoint[*start1].1 == v2_as_disjoint[*start2].1 {
                variable_shared_x += 1;
            }
        }
        if v1_as_disjoint[*start1].1 == v2_as_disjoint[*start2].1 {
            shared_x += 1;
        }
    }

    if shared_x > v1.len() || shared_x > v2.len(){
        dbg!(shared_x, v1.len(), &chain);
        panic!();
    }

    let l1 = chain[0].0;
    let r1 = chain[chain.len() - 1].0;
    let l2 = chain[0].1;
    let r2 = chain[chain.len() - 1].1;
    let start1 = l1.min(r1);
    let end1 = l1.max(r1);
    let start2 = l2.min(r2);
    let end2 = l2.max(r2);
    let frac_shared = variable_shared_x as f64 / variable_roots as f64;
    let tiger = TigdexOverlap {
        tig1: ind1,
        tig2: ind2,
        tig1_start: start1,
        tig1_end: end1,
        tig2_start: start2,
        tig2_end: end2,
        shared_tig: shared_x,
        variable_roots: variable_roots,
        variable_tigs: variable_shared_x,
        chain_reverse: chain_info.reverse,
    };
    return tiger;
}

fn smith_waterman(
    v1: &[u32],
    v2: &[u32],
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) -> (f64, usize, usize) {
    let len1 = v1.len();
    let len2 = v2.len();

    // Create a 2D vector to store scores
    let mut dp = vec![vec![0; len2 + 1]; len1 + 1];
    let mut traceback = vec![vec![0; len2 + 1]; len1 + 1];
    let mut max_score = 0;
    let mut max_index = (0, 0);

    // Fill the dp table
    for i in 1..=len1 {
        for j in 1..=len2 {
            let score_substitute = if v1[i - 1] == v2[j - 1] {
                match_score
            } else {
                mismatch_penalty
            };

            // Calculate possible scores for this cell
            let score_diag = dp[i - 1][j - 1] + score_substitute;
            let score_up = dp[i - 1][j] + gap_penalty;
            let score_left = dp[i][j - 1] + gap_penalty;

            // Cell score is the max of calculated scores or 0 (Smith-Waterman uses zero as a minimum score)
            dp[i][j] = 0.max(score_diag).max(score_up).max(score_left);

            // Update the traceback matrix
            if dp[i][j] == 0 {
                traceback[i][j] = 0;
            } else if dp[i][j] == score_diag {
                traceback[i][j] = 1;
            } else if dp[i][j] == score_up {
                traceback[i][j] = 2;
            } else {
                traceback[i][j] = 3;
            }

            max_index = if dp[i][j] > max_score {
                max_score = dp[i][j];
                (i, j)
            } else {
                max_index
            };
        }
    }

    let mut aln1 = vec![];
    let mut aln2 = vec![];
    while dp[max_index.0][max_index.1] > 0 {
        match traceback[max_index.0][max_index.1] {
            0 => break,
            1 => {
                aln1.push(v1[max_index.0 - 1]);
                aln2.push(v2[max_index.1 - 1]);
                max_index = (max_index.0 - 1, max_index.1 - 1);
            }
            2 => {
                aln1.push(v1[max_index.0 - 1]);
                aln2.push(0);
                max_index = (max_index.0 - 1, max_index.1);
            }
            3 => {
                aln1.push(0);
                aln2.push(v2[max_index.1 - 1]);
                max_index = (max_index.0, max_index.1 - 1);
            }
            _ => panic!("Invalid traceback value"),
        }
    }
    // The maximum score represents the optimal local alignment score
    (max_score as f64, aln1.len(), aln2.len())
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

fn find_exact_matches_indexes(seq1: &[usize], seq2: &[usize]) -> Vec<(usize, usize, usize)> {
    let mut matches = Vec::new();
    let index_seq = &seq2;
    let query_seq = &seq1;
    let mut index_map = FxHashMap::default();

    //Sorted
    for (i, &x) in index_seq.iter().enumerate() {
        index_map.entry(x).or_insert(vec![]).push(i);
    }

    //Sorted
    for (i, s) in query_seq.iter().enumerate() {
        if let Some(indices) = index_map.get(s) {
            if indices.len() > 500{
                continue;
            }
            for &j in indices {
                matches.push((i, j, 1));
            }
        }
    }

    matches
}

fn find_exact_matches_quadratic(seq1: &[usize], seq2: &[usize]) -> Vec<(usize, usize, usize)> {
    let mut matches = Vec::new();
    let len1 = seq1.len();
    let len2 = seq2.len();

    for i in 0..len1 {
        for j in 0..len2 {
            if seq1[i] == seq2[j] {
                matches.push((i, j, 1)); // (start in seq1, start in seq2, length of match = 1)
            }
        }
    }
    matches
}

fn dp_anchors(
    matches: &[(usize, usize, usize)],
    reverse: bool,
) -> (f64, Vec<(usize, usize, usize)>) {
    let mut dp = vec![0.; matches.len()];
    let mut prev = vec![None; matches.len()];
    let mut max_score = 0.;
    let mut max_index = 0;

    for i in 0..matches.len() {
        let (start1, start2, length) = matches[i];
        dp[i] = (5 * length) as f64;
        let back = if i > 50 { i - 50 } else { 0 };
        for j in back..i {
            let (end1, end2, len_prev) = matches[j];
            if reverse {
                if end1 >= start1 || end2 <= start2 {
                    continue;
                }
            } else {
                if end1 >= start1 || end2 >= start2 {
                    continue;
                }
            }
            let gap_penalty = (start1 as i32 - (end1 + len_prev - 1) as i32).abs()
                - (start2 as i32 - (end2 + len_prev - 1) as i32).abs();
            let score = dp[j] + 5. * length as f64 - gap_penalty.abs() as f64;
            if score > dp[i] {
                dp[i] = score;
                prev[i] = Some(j);
                if score > max_score {
                    max_score = score;
                    max_index = i;
                }
            }
        }
    }

    let mut chain = Vec::new();
    let mut i = Some(max_index);
    while let Some(idx) = i {
        chain.push(matches[idx]);
        i = prev[idx];
    }

    chain.reverse();
    (max_score, chain)
}

fn find_optimal_chain(seq1: &[usize], seq2: &[usize]) -> ChainInfo {
    let matches = find_exact_matches_indexes(seq1, seq2);

    if matches.is_empty() {
        return ChainInfo::default();
    }

    let (score_f, chain_f) = dp_anchors(&matches, false);
    let (score_r, chain_r) = dp_anchors(&matches, true);

    if score_f > score_r {
        ChainInfo {
            chain: chain_f,
            reverse: false,
            score: score_f,
        }
    } else {
        ChainInfo {
            chain: chain_r,
            reverse: true,
            score: score_r,
        }
    }
}

fn remove_contained_reads(all_reads_cat: &Vec<TigRead>, disjoint_set: &DisjointSet, used_ids: &FxHashSet<u64>) -> Vec<usize>
{
    let inverted_index_hashmap =
        all_reads_cat
            .iter()
            .enumerate()
            .fold(FxHashMap::default(), |mut acc, (i, x)| {
                for &y in x.tig_seq.iter() {
                    acc.entry(y).or_insert(vec![]).push(i);
                }
                acc
            });

    //open file for writing
    let mut bufwriter = BufWriter::new(File::create("histo.txt").unwrap());
    let mut bufwriter2 = BufWriter::new(File::create("histo.dbg").unwrap());

    //check contained reads
    let mut outer_reads = vec![];
    let mut contained_reads = FxHashSet::default();
    for (i, read) in all_reads_cat.iter().enumerate() {
        let mut contained = false;
        let mut index_count_map = FxHashMap::default();
        for &y in read.tig_seq.iter() {
            if let Some(indices) = inverted_index_hashmap.get(&y) {
                for &index in indices {
                    if index == i {
                        continue;
                    }
                    *index_count_map.entry(index).or_insert(0) += 1;
                }
            }
        }

        //look at top 5 indices and do disjoint_set distance
        let mut top_indices = index_count_map.iter().collect::<Vec<_>>();
        top_indices.sort_by(|a, b| b.1.cmp(a.1));
        for (index, count) in top_indices.iter() {
            if **count < 5 {
                break;
            }
            if contained_reads.contains(*index) {
                continue;
            }
            let read2 = &all_reads_cat[**index];
            dbg!("{} {} {} {}, {}", i, **index, &read.id, &read2.id);
            let tiglap = disjoint_distance(&read.tig_seq, &read2.tig_seq, &disjoint_set, &used_ids, i, **index);
            if read.id.contains("junk") || read2.id.contains("junk") {
                continue;
            }
            let name1 = read.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[0];
            let name2 = read2.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[0];
            let range1 = read.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[2];
            let range2 = read2.id.split_whitespace().collect::<Vec<&str>>()[1]
                .split(",")
                .collect::<Vec<&str>>()[2];
            let aln_len1 = tiglap.tig1_end - tiglap.tig1_start + 1;
            let aln_len2 = tiglap.tig2_end - tiglap.tig2_start + 1;
            let frac_shared = tiglap.shared_tig as f64 / aln_len1.max(aln_len2) as f64;
            let frac_shared_variable = tiglap.variable_tigs as f64 / tiglap.variable_roots as f64;

            let mut within = false;
            if name1 == name2 {
                within = true;
            }

            if !frac_shared_variable.is_nan()
                && aln_len1 as f64 > (read.tig_seq.len() as f64) * 0.4
                && aln_len2 as f64 > (read2.tig_seq.len() as f64) * 0.4
            {
                write!(
                    bufwriter2,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    tiglap.variable_tigs,
                    tiglap.variable_roots,
                    tiglap.shared_tig,
                    aln_len1,
                    within,
                    read.id,
                    read2.id
                )
                .unwrap();
                write!(
                    bufwriter,
                    "{}\t{}\t{}\n",
                    frac_shared, frac_shared_variable, within
                )
                .unwrap();
            }

            let span_read1 = tiglap.tig1_end - tiglap.tig1_start + 1;
            if frac_shared_variable > 0.8 && span_read1 as f64 > (read.tig_seq.len() as f64) * 0.9 {
                contained = true;
                log::trace!("{} {} CONTAINED", read.id, read2.id);
                contained_reads.insert(i);
            }
        }
        if !contained {
            outer_reads.push(i);
        }
    }

    log::info!("OUTER READS: {}", outer_reads.len());
    for i in outer_reads.iter() {
        log::info!("{}", all_reads_cat[*i].id);
    }



    return outer_reads

}

fn get_overlaps_outer_reads(
    all_reads_cat: &Vec<TigRead>,
    disjoint_set: &DisjointSet,
    used_ids: &FxHashSet<u64>,
    outer_reads: &Vec<usize>,
) -> Vec<TigdexOverlap> {
    let mut bufwriter = BufWriter::new(File::create("overlaps.txt").unwrap());
    let inverted_index_hashmap_outer = outer_reads
        .iter()
        .map(|&i| {
            all_reads_cat[i]
                .tig_seq
                .iter()
                .map(|&x| (x, i))
                .collect::<Vec<(u32, usize)>>()
        })
        .flatten()
        .fold(FxHashMap::default(), |mut acc, x| {
            acc.entry(x.0).or_insert(vec![]).push(x.1);
            acc
        });
    let mut overlaps = vec![];

    let mut compared_reads = FxHashSet::default();
    for (i, read) in all_reads_cat.iter().enumerate() {
        let mut index_count_map = FxHashMap::default();
        if !outer_reads.contains(&i) {
            continue;
        }
        for &y in read.tig_seq.iter() {
            if let Some(indices) = inverted_index_hashmap_outer.get(&y) {
                for &index in indices {
                    if index == i {
                        continue;
                    }
                    *index_count_map.entry(index).or_insert(0) += 1;
                }
            }
        }
        //sort
        let mut top_indices = index_count_map.into_iter().collect::<Vec<_>>();
        top_indices.sort_by(|a, b| b.1.cmp(&a.1));

        for (index, count) in top_indices.iter() {
            let sorted_readpair = if i < *index { (i, *index) } else { (*index, i) };
            if compared_reads.contains(&sorted_readpair) {
                continue;
            }
            let read2 = &all_reads_cat[*index];
            let tiglap = disjoint_distance(&read.tig_seq, &read2.tig_seq, &disjoint_set, &used_ids, i, *index);
            let frac_shared_variable = tiglap.variable_tigs as f64 / tiglap.variable_roots as f64;
            let span_read1 = tiglap.tig1_end - tiglap.tig1_start + 1;
            writeln!(bufwriter,
                "intersect tigs: {} i: {} j: {} fsv: {} leni: {} {}-{} lenj: {} {}-{} shared_tigs: {} REVERSE:{}",
                count,
                i,
                *index,
                frac_shared_variable,
                read.tig_seq.len(),
                tiglap.tig1_start,
                tiglap.tig1_end,
                read2.tig_seq.len(),
                tiglap.tig2_start,
                tiglap.tig2_end,
                tiglap.shared_tig,
                tiglap.chain_reverse
            );
            if frac_shared_variable > 0.8 {
                overlaps.push(tiglap);
            }
            compared_reads.insert(sorted_readpair);
        }
    }

    return overlaps;
}

fn read_graph_from_overlaps(all_reads_cat: Vec<TigRead>, overlaps: &Vec<TigdexOverlap>) -> OverlapGraph
{
    let mut nodes = FxHashMap::default();
    let mut edges = vec![];

    for tiglap in overlaps.iter() {
        let i_tig = tiglap.tig1;
        let j_tig = tiglap.tig2;
        let read = &all_reads_cat[tiglap.tig1].tig_seq;
        let read2 = &all_reads_cat[tiglap.tig2].tig_seq;
        let frac_shared_variable = tiglap.variable_tigs as f64 / tiglap.variable_roots as f64;

        //check if end-to-end overlap
        let mut forward1 = false;
        let mut forward2 = false;
        let mut r1_r2 = true;
        let mut ol = false;

        if tiglap.chain_reverse {
            if tiglap.tig1_start < read.len() / 10 && tiglap.tig2_start < read2.len() / 10 {
                forward1 = false;
                forward2 = true;
                r1_r2 = true;
                ol = true;
            } else if tiglap.tig1_end > read.len() * 9 / 10
                && tiglap.tig2_end > read2.len() * 9 / 10
            {
                forward1 = true;
                forward2 = false;
                r1_r2 = true;
                ol = true;
            }
        } else {
            if tiglap.tig1_start < read.len() / 10 && tiglap.tig2_end > read2.len() * 9 / 10 {
                forward1 = true;
                forward2 = true;
                r1_r2 = false;
                ol = true;
            } else if tiglap.tig2_start < read2.len() / 10
                && tiglap.tig1_end > read.len() * 9 / 10
            {
                forward1 = true;
                forward2 = true;
                r1_r2 = true;
                ol = true;
            }
        }
        let aln_len1 = tiglap.tig1_end - tiglap.tig1_start + 1;
        let aln_len2 = tiglap.tig2_end - tiglap.tig2_start + 1;
        if ol && frac_shared_variable > 0.8 && aln_len1.max(aln_len2) > 20 {
            log::info!(
                "OVERLAP {} {} {} {} {}-{} {} {}-{}, REVERSE: {}",
                i_tig,
                j_tig,
                frac_shared_variable,
                read.len(),
                tiglap.tig1_start,
                tiglap.tig1_end,
                read2.len(),
                tiglap.tig2_start,
                tiglap.tig2_end,
                tiglap.chain_reverse
            );

            let start = if r1_r2 { i_tig } else { j_tig };
            let end = if r1_r2 { j_tig } else { i_tig };
            let new_read_overlap = ReadOverlapEdge {
                node1: start,
                node2: end,
                forward1,
                forward2,
                overlap_len_bases: 0,
                overlap_len_tigs: aln_len1.max(aln_len2),
                shared_tigs: tiglap.shared_tig,
                variable_tigs: tiglap.variable_tigs,
                variable_roots: tiglap.variable_roots,
            };

            edges.push(Some(new_read_overlap));
            let ind = edges.len() - 1;
            {
                let rd1 = nodes.entry(i_tig).or_insert(ReadData::default());
                rd1.index = i_tig;
                if r1_r2 && forward1 {
                    rd1.out_edges.push(ind)
                } else if !r1_r2 && !forward2 {
                    rd1.out_edges.push(ind)
                } else {
                    rd1.in_edges.push(ind)
                }
            }
            {
                let rd2 = nodes.entry(j_tig).or_insert(ReadData::default());
                rd2.index = j_tig;
                if r1_r2 && forward2 {
                    rd2.in_edges.push(ind)
                } else if !r1_r2 && !forward1 {
                    rd2.in_edges.push(ind)
                } else {
                    rd2.out_edges.push(ind)
                }
            }
        }
    }

    let graph = OverlapGraph { reads: all_reads_cat, nodes, edges };

    return graph;
}
