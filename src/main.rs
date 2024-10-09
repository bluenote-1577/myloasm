use phaller::bam_parsing;
use phaller::seeding;
use fxhash::FxHashMap;
use rust_htslib::bam;
use rust_htslib::bam::{Read, Reader, Record};
use needletail::*;
use std::env;

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
    //let fastq_files = vec!["/home/jshaw/software/phaller/temp/98acc-10xcov-mn.fq", "/home/jshaw/software/phaller/temp/98acc-10xcov-mgh.fq"];
    //let fastq_files = vec!["/home/jshaw/software/phaller/temp/99.8acc-10x100kb-mat.fq", "/home/jshaw/software/phaller/temp/99.8acc-10x100kb-pat.fq"];
    //let fastq_files = vec!["/home/jshaw/software/phaller/temp/99acc-10x100kb-mat.fq", "/home/jshaw/software/phaller/temp/99acc-10x100kb-pat.fq"];
    let fastq_files = vec!["/home/jshaw/software/phaller/temp/98acc-10x100kb-mat.fq", "/home/jshaw/software/phaller/temp/98acc-10x100kb-pat.fq"];
    let cuttlefish_file = "/home/jshaw/software/phaller/temp/cuttle-hg002-99.8.fa";
   //#let bam_file = env::args().nth(1).unwrap();

    //time this operation and print the time
    let start = std::time::Instant::now();
    let (kmer_to_unitig, unitig_vec) = parse_unitigs_into_table(cuttlefish_file);
    let mut all_reads = vec![];
    for fastq_file in fastq_files{
        let reads = read_to_unitig_count_vector(fastq_file, &kmer_to_unitig, &unitig_vec);
        all_reads.push(reads);
    }
    //    bam_parsing::get_snvs(bam_file);
    //bam_parsing::get_pileup(&bam_file, &fasta_file);
    let duration = start.elapsed();
    log::info!("Time elapsed in parse_bam() is: {:?}", duration);
    for i in 0..50{
        for j in i..50{
            let read1 = &all_reads[0][i];
            let read2 = &all_reads[1][j];
            let (distance, len1, len2) = smith_waterman(&read1.0, &read2.0, 5, -3, -1);
            if len1 > 10 {
                log::info!("ACROSS");
                log::info!("{}", read1.1);
                log::info!("{}", read2.1);
                let measure = distance as f64 / (len1.max(len2) as f64);
                log::info!("{} {} {}", len1, len2, measure);
                let range1 = read1.1.split_whitespace().collect::<Vec<&str>>()[1].split(",").collect::<Vec<&str>>()[2];
                let range2 = read2.1.split_whitespace().collect::<Vec<&str>>()[1].split(",").collect::<Vec<&str>>()[2];
                println!("ACROSS\t{}\t{}\t{}", measure, range1, range2);
            }
        }
    }

    for i in 0..50{
        for j in i+1..50{
            let read1 = &all_reads[0][i];
            let read2 = &all_reads[0][j];
            let (distance, len1, len2) = smith_waterman(&read1.0, &read2.0, 5, -3, -1);
            if len1 > 10 {
                log::info!("WITHIN");
                log::info!("{}", read1.1);
                log::info!("{}", read2.1);
                let measure = distance as f64 / (len1.max(len2) as f64);
                log::info!("{} {} {}", len1, len2, measure);
                let range1 = read1.1.split_whitespace().collect::<Vec<&str>>()[1].split(",").collect::<Vec<&str>>()[2];
                let range2 = read2.1.split_whitespace().collect::<Vec<&str>>()[1].split(",").collect::<Vec<&str>>()[2];
                println!("WITHIN\t{}\t{}\t{}", measure, range1,range2);
            }
        }
    }
}

fn parse_unitigs_into_table(cuttlefish_file: &str) -> (FxHashMap<u64, u32>, Vec<Vec<u8>>)
{
    let mut kmer_to_unitig_count : FxHashMap<u64, u32> = fxhash::FxHashMap::default();
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
            count += 1;
        }
    }
    return (kmer_to_unitig_count,unitig_vec);
}

fn read_to_unitig_count_vector(fastq_file: &str, kmer_to_unitig_count: &FxHashMap<u64, u32>, unitig_vec: &Vec<Vec<u8>>) -> Vec<(Vec<u32>, String)>
{
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
            reads.push((unitig_string, String::from_utf8_lossy(rec.id()).to_string()));
        }
    }
    return reads
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

fn smith_waterman(v1: &[u32], v2: &[u32], match_score: i32, mismatch_penalty: i32, gap_penalty: i32) -> (i32, usize, usize) {
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
            },
            2 => {
                aln1.push(v1[max_index.0 - 1]);
                aln2.push(0);
                max_index = (max_index.0 - 1, max_index.1);
            },
            3 => {
                aln1.push(0);
                aln2.push(v2[max_index.1 - 1]);
                max_index = (max_index.0, max_index.1 - 1);
            },
            _ => panic!("Invalid traceback value"),
        }
    }
    // The maximum score represents the optimal local alignment score
    (max_score, aln1.len(), aln2.len())
}
