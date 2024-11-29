use dashmap::DashMap;
use smallvec::SmallVec;
use smallvec::smallvec;
use crate::cbloom;
use crate::cli::Cli;
use crate::twin_graph;
use std::path::Path;
use rayon::prelude::*;
use fxhash::FxHashMap;
use crate::mapping;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;
use crate::types::*;
use crate::seeding;
use fishers_exact::fishers_exact;
use fxhash::FxHashSet;
use crate::mapping::*;
use std::io::{BufWriter, Write};
use std::io::BufReader;


fn homopolymer_compression(seq: Vec<u8>) -> Vec<u8> {
    let mut compressed_seq = vec![];
    let mut last_base = seq[0];
    let mut _count = 1;
    for i in 1..seq.len() {
        if seq[i] == last_base {
            _count += 1;
        } else {
            compressed_seq.push(last_base);
            last_base = seq[i];
            _count = 1;
        }
    }
    return compressed_seq;
}

pub fn read_to_split_kmers(
    k: usize,
    threads: usize,
    args: &Cli,
) -> DashMap<u64,[u32;2], MMBuildHasher>{

    let homopolymer_comp = args.homopolymer_compression;
    let bf_size = args.bloom_filter_size * 1_000_000_000.;
    let map = Arc::new(DashMap::with_hasher(MMBuildHasher::default()));

    log::info!("Populating bloom filter of size {} GB and counting kmers.", args.bloom_filter_size);
    let filter = Arc::new(cbloom::Filter::with_size_and_hashers(bf_size as usize, 3));

    if bf_size > 0.{
        for fastq_file in args.input_files.iter(){
            let (mut tx, rx) = spmc::channel();
            let file = fastq_file.to_owned();
            thread::spawn(move || {
                let bufreader = BufReader::new(std::fs::File::open(file).expect("valid path"));
                let mut reader = needletail::parse_fastx_reader(bufreader).expect("valid path");
                while let Some(record) = reader.next() {
                    let rec = record.expect("Error reading record");
                    let seq;
                    if homopolymer_comp {
                        seq = homopolymer_compression(rec.seq().to_vec());
                    } else {
                        seq = rec.seq().to_vec();
                    }
                    tx.send(seq).unwrap();
                }
            });

            let mut handles = Vec::new();
            for _ in 0..threads {
                let filter = Arc::clone(&filter);
                let map = Arc::clone(&map);
                let rx = rx.clone();
                let mut myhashmap = FxHashSet::default();
                handles.push(thread::spawn(move || {
                    loop{
                        match rx.recv() {
                            Ok(msg) => {
                                let split_kmer_info = seeding::split_kmer_mid(msg, k);
                                for kmer_i in split_kmer_info.iter() {
                                    if filter.maybe_insert(kmer_i.full_kmer) {
                                        myhashmap.insert(kmer_i.full_kmer);
                                    }
                                    //Batching might help with repetitive k-mers
                                    if myhashmap.len() > 300_000{
                                        for key in myhashmap.into_iter(){
                                            map.insert(key, [0,0]);
                                        }
                                        myhashmap = FxHashSet::default();
                                    }
                                }
                            }
                            Err(_) => {
                                for key in myhashmap.into_iter(){
                                    map.insert(key, [0,0]);
                                }
                                // When sender is dropped, recv will return an Err, and we can break the loop
                                break;
                            }
                        }
                    }
                }));
            }

            for handle in handles {
                handle.join().unwrap();
            }
        }
    }

    drop(filter);

    log::info!("Bloom filter populated. Counting kmers.");
    //Round 2, actually count
    for fastq_file in args.input_files.iter(){
        let (mut tx, rx) = spmc::channel();
        let file = fastq_file.to_owned();
        thread::spawn(move || {
            let bufreader = BufReader::new(std::fs::File::open(file).expect("valid path"));
            let mut reader = needletail::parse_fastx_reader(bufreader).expect("valid path");
            while let Some(record) = reader.next() {
                let rec = record.expect("Error reading record");
                let seq;
                if homopolymer_comp {
                    seq = homopolymer_compression(rec.seq().to_vec());
                } else {
                    seq = rec.seq().to_vec();
                }
                tx.send(seq).unwrap();
            }
        });

        let mut handles = Vec::new();
        for _ in 0..threads {
            let rx = rx.clone();
            let map = Arc::clone(&map);
            handles.push(thread::spawn(move || {
                loop{
                    match rx.recv() {
                        Ok(msg) => {
                            //Querying the hash table takes > 50x longer than splitting the kmer
                            let _start = std::time::Instant::now();
                            let split_kmer_info = seeding::split_kmer_mid(msg, k);
                            //println!("Split kmer time: {:?}", start.elapsed());
                            let _start = std::time::Instant::now();
                            if bf_size > 0.{
                                for kmer_i in split_kmer_info.iter() {
                                    if let Some(mut counts) = map.get_mut(&kmer_i.full_kmer){
                                        if kmer_i.canonical{
                                            counts[0] += 1;
                                        }
                                        else{
                                            counts[1] += 1;
                                        }
                                    }   
                                }
                            }
                            else{
                                for kmer_i in split_kmer_info.iter(){
                                    let mut counts = map.entry(kmer_i.full_kmer).or_insert([0,0]);
                                    if kmer_i.canonical{
                                        counts[0] += 1;
                                    }
                                    else{
                                        counts[1] += 1;
                                    }
                                }
                            }
                            //println!("Populating table time: {:?}", start.elapsed());
                        }
                        Err(_) => {
                            // When sender is dropped, recv will return an Err, and we can break the loop
                            break;
                        }
                    }
                }
            }));
        }
        for handle in handles {
            handle.join().unwrap();
        }
    }
    
    let map_size_raw = map.len();
    map.retain(|_, v| v[0] > 0 && v[1] > 0);
    let map_size_retain = map.len();
    log::info!("Removed {} kmers with counts < 1 in both strands.", map_size_raw - map_size_retain);
    if map_size_retain < map_size_raw / 1000 {
        log::warn!("Less than 0.1% of kmers have counts > 1 in both strands. This may indicate a problem with the input data or very low coverage.");
    }

    return Arc::try_unwrap(map).unwrap();
}

pub fn twin_reads_from_snpmers(kmer_info: &KmerGlobalInfo, fastq_files: &[String], args: &Cli) -> Vec<TwinRead>{

    let mut snpmer_set = FxHashSet::default();
    for snpmer_i in kmer_info.snpmer_info.iter(){
        let k = snpmer_i.k as usize;
        let snpmer1 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[0] as u64) << (k-1) );
        let snpmer2 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[1] as u64) << (k-1) );
        snpmer_set.insert(snpmer1);
        snpmer_set.insert(snpmer2);
    }

    let snpmer_set = Arc::new(snpmer_set);
    let twin_read_vec = Arc::new(Mutex::new(vec![]));

    let files_owned = fastq_files.iter().map(|x| x.to_string()).collect::<Vec<String>>();
    let hpc = args.homopolymer_compression;

    for fastq_file in files_owned{
        let (mut tx, rx) = spmc::channel();
        thread::spawn(move || {
            let mut reader = needletail::parse_fastx_file(fastq_file).expect("valid path");
            while let Some(record) = reader.next() {
                let rec = record.expect("Error reading record");
                let seq;
                if hpc{
                    seq = homopolymer_compression(rec.seq().to_vec());
                } else {
                    seq = rec.seq().to_vec();
                }
                if seq.len() < 1000{
                    continue;
                }
                let id = String::from_utf8_lossy(rec.id()).to_string();
                if let Some(qualities) = rec.qual(){
                    tx.send((seq, Some(qualities.to_vec()), id)).unwrap();
                }
                else{
                    tx.send((seq, None, id)).unwrap();
                }
            }
        });

        let mut handles = Vec::new();
        let k = args.kmer_size;
        let c = args.c;
        for _ in 0..args.threads{
            let rx = rx.clone();
            let set = Arc::clone(&snpmer_set);
            let twrv = Arc::clone(&twin_read_vec);
            handles.push(thread::spawn(move || {
                loop{
                    match rx.recv() {
                        Ok(msg) => {
                            let seq = msg.0;
                            let qualities = msg.1;
                            let id = msg.2;
                            let twin_read = seeding::get_twin_read(seq, qualities, k, c, set.as_ref(), id);
                            if twin_read.is_some(){
                                let mut vec = twrv.lock().unwrap();
                                vec.push(twin_read.unwrap());
                            }
                        }
                        Err(_) => {
                            // When sender is dropped, recv will return an Err, and we can break the loop
                            break;
                        }
                    }
                }
            }));
        }

        for handle in handles {
            handle.join().unwrap();
        }
    }

    let mut twin_reads = Arc::try_unwrap(twin_read_vec).unwrap().into_inner().unwrap();
    twin_reads.sort_by(|a,b| a.id.cmp(&b.id));
    let solid = &kmer_info.solid_kmers;
    for tr in twin_reads.iter_mut(){
        let mut new_mini = vec![];
        for mini in tr.minimizers.iter(){
            if solid.contains(&mini.1){
                new_mini.push(*mini);
            }
        }
        tr.minimizers = new_mini;
    }

    if log::log_enabled!(log::Level::Trace) {
        for twin_read in twin_reads.iter(){
            let decoded_snpmers = twin_read.snpmers.iter().map(|x| (x.0, decode_kmer(x.1, twin_read.k))).collect::<Vec<_>>();
            log::trace!("{} {:?}", twin_read.id, decoded_snpmers);
        }
    }

    let number_reads_below_threshold = twin_reads.iter().filter(|x| x.est_id.is_some() && x.est_id.unwrap() < args.quality_value_cutoff).count();
    log::info!("Number of valid reads - {}. Number of reads below quality threshold - {}.", twin_reads.len(), number_reads_below_threshold);

    return twin_reads;
}


#[inline]
pub fn split_kmer(kmer: u64, k: usize) -> (u64, u8){
    let split_mask_extract = 3 << (k-1);
    let mid_base = (kmer & split_mask_extract) >> (k-1);
    let masked_kmer = kmer & !split_mask_extract;
    return (masked_kmer, mid_base as u8);
}

pub fn get_snpmers(big_kmer_map: DashMap<Kmer64, [u32;2], MMBuildHasher>, k: usize, _args: &Cli) -> KmerGlobalInfo{

    log::info!("Finding snpmers...");
    let mut new_map_counts_bases : FxHashMap<Kmer64, CountsAndBases> = FxHashMap::default();
    let mut kmer_counts = vec![];
    let mut solid_kmers = FxHashSet::default();

    for pair in big_kmer_map.iter(){
        let counts = pair.value();
        kmer_counts.push(counts[0] + counts[1]);
    }
    kmer_counts.sort_unstable();
    let high_freq_thresh = kmer_counts[kmer_counts.len() - kmer_counts.len() / 50000];
    log::debug!("High frequency threshold: {}", high_freq_thresh);
    drop(kmer_counts);

    for pair in big_kmer_map.into_iter(){
        let kmer = pair.0;
        let (split_kmer, mid_base) = split_kmer(kmer, k);
        let counts = pair.1;
        let count = counts[0] + counts[1];
        if count > 2{
            let v = new_map_counts_bases.entry(split_kmer).or_insert(CountsAndBases{counts: SmallVec::new(), bases: SmallVec::new()});
            v.counts.push(counts);
            v.bases.push(mid_base);
            if count < high_freq_thresh{
                solid_kmers.insert(kmer);
            }
        }
    }


    let potential_snps = Mutex::new(0);
    let snpmers = Mutex::new(vec![]);
    new_map_counts_bases.into_par_iter().for_each(|(split_kmer, c_and_b)|{
        let mut counts = c_and_b.counts;
        let bases = c_and_b.bases;
        if counts.len() > 1{
            counts.sort_unstable_by(|a, b| (b[0] + b[1]).cmp(&(a[0] + a[1])));

            //Errors are differentiated because they will have > 2 alleles
            //and the smallest alleles will have a low count. 
            let n = counts[0][0] + counts[0][1];
            let succ = counts[1][0] + counts[1][1];
            let right_p_val_thresh1 = twin_graph::binomial_test(n as u64 , succ as u64,0.025);
            let right_p_val_thresh2 = twin_graph::binomial_test(n as u64 , succ as u64,0.050);
            let cond1 = right_p_val_thresh1 > 0.05;
            let cond2 = right_p_val_thresh2 > 0.05 && k < 5;
            if cond1 || cond2 {
                log::trace!("NOT SNPMER BINOMIAL c:{:?} c:{:?}",  counts[0], counts[1]);
                return;
            }

            //Add pseudocount... especially when all reads are 
            //already in forward strand (happens when debugging or manipulation via samtools)
            let a = counts[0][0];
            let b = counts[1][0];
            let c = counts[0][1];
            let d = counts[1][1];
            let contingency_table = [
                a.max(c), b.max(d),
                c.min(a), d.min(b)
            ];
            let p_value = fishers_exact(&contingency_table).unwrap().two_tail_pvalue;
            let odds;
            if contingency_table[0] == 0 || contingency_table[1] == 0 || contingency_table[2] == 0 || contingency_table[3] == 0 {
                odds = 0.0;
            } else {
                odds = (contingency_table[0] as f64 * contingency_table[3] as f64) / (contingency_table[1] as f64 * contingency_table[2] as f64);
            }

            if odds == 0.{
                return;
            }

            //Is snpmer
            if p_value > 0.005 || (odds < 1.5 && odds > 1./1.5){
                let mid_bases = bases;
                let snpmer = SnpmerInfo{
                    split_kmer: split_kmer,
                    mid_bases: smallvec![mid_bases[0], mid_bases[1]],
                    counts: smallvec![counts[0][0] + counts[0][1], counts[1][0] + counts[1][1]],
                    k: k as u8,
                };
                snpmers.lock().unwrap().push(snpmer);

                let snpmer1 = split_kmer as u64 | ((mid_bases[0] as u64) << (k-1));
                let snpmer2 = split_kmer as u64 | ((mid_bases[1] as u64) << (k-1));
                log::trace!("{} c:{:?} {} c:{:?}, p:{}, odds:{}", decode_kmer(snpmer1, k as u8), counts[0], decode_kmer(snpmer2, k as u8), counts[1], p_value, odds);
                *potential_snps.lock().unwrap() += 1;

            }
            else{
                log::trace!("NOT SNPMER c:{:?} c:{:?}, p:{}, odds:{}",  counts[0], counts[1], p_value, odds);
            }
        }
    });

    let mut snpmers = snpmers.into_inner().unwrap();
    snpmers.sort();
    log::info!("Number of snpmers: {}. ", potential_snps.into_inner().unwrap());
    return KmerGlobalInfo{
        snpmer_info: snpmers,
        solid_kmers: solid_kmers,
        high_freq_thresh: high_freq_thresh as f64,
    };
}

pub fn parse_unitigs_into_table(cuttlefish_file: &str) -> (FxHashMap<u64, u32>, Vec<Vec<u8>>) {
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


fn split_read(twin_read: TwinRead, mut break_points: Vec<Breakpoints>) -> Vec<TwinRead>{
    let mut new_reads = vec![];
    break_points.push(Breakpoints{pos1: twin_read.base_length, pos2: twin_read.base_length, cov: 0});
    let mut last_break = 0;
    for (i,break_point) in break_points.iter().enumerate(){
        let bp_start = break_point.pos1;
        let bp_end = break_point.pos2;
        if bp_start - last_break > 1000{
            let mut new_read = TwinRead::default();
            new_read.minimizers = twin_read.minimizers.iter().filter(|x| x.0 >= last_break && x.0 < bp_start).copied().map(|x| (x.0 - last_break, x.1)).collect();
            new_read.snpmers = twin_read.snpmers.iter().filter(|x| x.0 >= last_break && x.0 < bp_start).copied().map(|x| (x.0 - last_break, x.1)).collect();
            new_read.id = format!("{}+split{}", &twin_read.id, i);
            log::trace!("Split read {} at {}-{}", &new_read.id, last_break, bp_start);
            new_read.k = twin_read.k;
            new_read.base_length = bp_start - last_break;
            new_read.dna_seq = twin_read.dna_seq[last_break..bp_start].to_owned();
            new_reads.push(new_read);
        }
        last_break = bp_end;
    }
    return new_reads;
}

pub fn split_outer_reads(mut twin_reads: Vec<TwinRead>, tr_map_info: Vec<TwinReadMapping>, args: &Cli)
-> (Vec<TwinRead>, Vec<usize>){
    let tr_map_info_dict = tr_map_info.into_iter().map(|x| (x.tr_index, x)).collect::<FxHashMap<usize, TwinReadMapping>>();
    let mut new_twin_reads = vec![];
    let mut new_outer_indices = vec![];
    let cov_file = Path::new(args.output_dir.as_str()).join("read_coverages.txt");
    let mut writer = BufWriter::new(std::fs::File::create(cov_file).unwrap());
    for (i, twin_read) in twin_reads.iter_mut().enumerate(){
        if tr_map_info_dict.contains_key(&i){
            let map_info = tr_map_info_dict.get(&i).unwrap();
            let twin_read = std::mem::take(twin_read);

            let breakpoints = mapping::cov_mapping_breakpoints(map_info);

            if log::log_enabled!(log::Level::Trace) {
                let depths = map_info.mapping_boundaries().depth().collect::<Vec<_>>();
                for depth in depths{
                    let mut string = format!("{} {}-{} COV:{}, BREAKPOINTS:", twin_read.id, depth.start, depth.stop, depth.val);
                    for breakpoint in breakpoints.iter(){
                        string.push_str(format!("--{} to {}--", breakpoint.pos1, breakpoint.pos2).as_str());
                    }
                    writeln!(writer, "{}", &string).unwrap();
                }
            }
            if breakpoints.len() == 0{
                new_outer_indices.push(new_twin_reads.len());
                new_twin_reads.push(twin_read);
                continue;
            }
            else{
                let splitted_reads = split_read(twin_read, breakpoints);
                for new_read in splitted_reads{
                    new_outer_indices.push(new_twin_reads.len());
                    new_twin_reads.push(new_read);
                }
            }
        }
        else{
            new_twin_reads.push(std::mem::take(twin_read));
        }
    }
    return (new_twin_reads, new_outer_indices);
}