use std::collections::HashSet;
use smallvec::SmallVec;
use smallvec::smallvec;
use crate::cli::Cli;
use crate::constants::MAX_FRACTION_OF_SNPMERS_IN_READ;
use crate::constants::MIN_READ_LENGTH;
use crate::twin_graph;
use rayon::prelude::*;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;
use crate::types::*;
use crate::seeding;
use fishers_exact::fishers_exact;
use std::path::Path;


pub fn homopolymer_compression(seq: Vec<u8>) -> Vec<u8> {
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

// pub fn read_to_split_kmers(
//     k: usize,
//     threads: usize,
//     args: &Cli,
// ) -> Vec<(u64,[u32;2])>{

//     let homopolymer_comp = args.homopolymer_compression;
//     let bf_size = args.bloom_filter_size * 1_000_000_000.;
//     let map = Arc::new(DashMap::with_hasher(MMBuildHasher::default()));

//     log::info!("Populating bloom filter of size {} GB and counting kmers.", args.bloom_filter_size);
//     let filter = Arc::new(cbloom::Filter::with_size_and_hashers(bf_size as usize, 3));

//     if bf_size > 0.{
//         let counter = Arc::new(Mutex::new(0));
//         for fastq_file in args.input_files.iter(){
//             let (mut tx, rx) = spmc::channel();
//             let file = fastq_file.to_owned();
//             thread::spawn(move || {
//                 let bufreader = BufReader::new(std::fs::File::open(file).expect("valid path"));
//                 let mut reader = needletail::parse_fastx_reader(bufreader).expect("valid path");
//                 while let Some(record) = reader.next() {
//                     let rec = record.expect("Error reading record");
//                     let seq;
//                     if homopolymer_comp {
//                         seq = homopolymer_compression(rec.seq().to_vec());
//                     } else {
//                         seq = rec.seq().to_vec();
//                     }
//                     tx.send(seq).unwrap();
//                 }
//             });

//             let mut handles = Vec::new();
//             for _ in 0..threads {
//                 let filter = Arc::clone(&filter);
//                 let map = Arc::clone(&map);
//                 let rx = rx.clone();
//                 let mut myhashmap = FxHashSet::default();
//                 let clone_counter = Arc::clone(&counter);
//                 handles.push(thread::spawn(move || {
//                     loop{
//                         match rx.recv() {
//                             Ok(msg) => {
//                                 let split_kmer_info = seeding::split_kmer_mid(msg, k);
//                                 {
//                                     let mut counter = clone_counter.lock().unwrap();
//                                     *counter += 1;
//                                     if *counter % 10000 == 0{
//                                         log::info!("Processed {} reads.", counter);
//                                     }
//                                 }
//                                 for kmer_i in split_kmer_info.iter() {
//                                     if filter.maybe_insert(kmer_i.full_kmer) {
//                                         myhashmap.insert(kmer_i.full_kmer);
//                                     }
//                                     //Batching might help with repetitive k-mers
//                                     if myhashmap.len() > 300_000{
//                                         for key in myhashmap.into_iter(){
//                                             map.insert(key, [0,0]);
//                                         }
//                                         myhashmap = FxHashSet::default();
//                                     }
//                                 }
//                             }
//                             Err(_) => {
//                                 for key in myhashmap.into_iter(){
//                                     map.insert(key, [0,0]);
//                                 }
//                                 // When sender is dropped, recv will return an Err, and we can break the loop
//                                 break;
//                             }
//                         }
//                     }
//                 }));
//             }

//             for handle in handles {
//                 handle.join().unwrap();
//             }
//         }
//     }

//     log::debug!("Hashmap capacity: {}", map.capacity());
//     log::debug!("Hashmap len: {}", map.len());
//     log::debug!("Memory usage: {:?} MB", memory_stats().unwrap().physical_mem as f32 / 1_000_000.);

//     drop(filter);
//     map.shrink_to_fit();
//     log::debug!("Shrunk Hashmap capacity: {}", map.capacity());
//     log::debug!("Shrunk Hashmap len: {}", map.len());
//     log::debug!("Memory usage: {:?} MB", memory_stats().unwrap().physical_mem as f32 / 1_000_000.);

//     log::info!("Bloom filter populated. Counting kmers.");
//     let counter = Arc::new(Mutex::new(0));
//     //Round 2, actually count
//     for fastq_file in args.input_files.iter(){
//         let (mut tx, rx) = spmc::channel();
//         let file = fastq_file.to_owned();
//         thread::spawn(move || {
//             let bufreader = BufReader::new(std::fs::File::open(file).expect("valid path"));
//             let mut reader = needletail::parse_fastx_reader(bufreader).expect("valid path");
//             while let Some(record) = reader.next() {
//                 let rec = record.expect("Error reading record");
//                 let seq;
//                 if homopolymer_comp {
//                     seq = homopolymer_compression(rec.seq().to_vec());
//                 } else {
//                     seq = rec.seq().to_vec();
//                 }
//                 tx.send(seq).unwrap();
//             }
//         });

//         let mut handles = Vec::new();
//         for _ in 0..threads {
//             let clone_counter = Arc::clone(&counter);
//             let rx = rx.clone();
//             let map = Arc::clone(&map);
//             handles.push(thread::spawn(move || {
//                 loop{
//                     match rx.recv() {
//                         Ok(msg) => {
//                             //Querying the hash table takes > 50x longer than splitting the kmer
//                             let _start = std::time::Instant::now();
//                             let split_kmer_info = seeding::split_kmer_mid(msg, k);
//                             {
//                                 let mut counter = clone_counter.lock().unwrap();
//                                 *counter += 1;
//                                 if *counter % 10000 == 0{
//                                     log::info!("Processed {} reads.", counter);
//                                 }
//                             }
//                             //println!("Split kmer time: {:?}", start.elapsed());
//                             let _start = std::time::Instant::now();
//                             if bf_size > 0.{
//                                 for kmer_i in split_kmer_info.iter() {
//                                     if let Some(mut counts) = map.get_mut(&kmer_i.full_kmer){
//                                         if kmer_i.canonical{
//                                             counts[0] += 1;
//                                         }
//                                         else{
//                                             counts[1] += 1;
//                                         }
//                                     }   
//                                 }
//                             }
//                             else{
//                                 for kmer_i in split_kmer_info.iter(){
//                                     let mut counts = map.entry(kmer_i.full_kmer).or_insert([0,0]);
//                                     if kmer_i.canonical{
//                                         counts[0] += 1;
//                                     }
//                                     else{
//                                         counts[1] += 1;
//                                     }
//                                 }
//                             }
//                             //println!("Populating table time: {:?}", start.elapsed());
//                         }
//                         Err(_) => {
//                             // When sender is dropped, recv will return an Err, and we can break the loop
//                             break;
//                         }
//                     }
//                 }
//             }));
//         }
//         for handle in handles {
//             handle.join().unwrap();
//         }
//     }

//     log::debug!("Final Hashmap capacity: {}", map.capacity());
//     log::debug!("Final Hashmap len: {}", map.len());
//     log::debug!("Final Hash Memory usage: {:?} MB", memory_stats().unwrap().physical_mem as f32 / 1_000_000.);
    
//     let map_size_raw = map.len();
//     let map = Arc::try_unwrap(map).unwrap();
//     let new_map = map.into_iter().filter(|(_, v)| v[0] > 0 && v[1] > 0 && v[0] + v[1] > 2).collect::<Vec<(u64,[u32;2])>>();
//     let map_size_retain = new_map.len();
//     log::info!("Removed {} kmers with counts < 1 in both strands and <= 3 multiplicity.", map_size_raw - map_size_retain);
//     if map_size_retain < map_size_raw / 1000 {
//         log::warn!("Less than 0.1% of kmers have counts > 1 in both strands and > 2 multiplicity. This may indicate a problem with the input data or very low coverage.");
//     }
//     log::debug!("Final Hashmap capacity after vectorization: {}", new_map.capacity());
//     log::debug!("Final Hashmap len after vectorization: {}", new_map.len());
//     log::debug!("Final Hash Memory usage after vectorization : {:?} MB", memory_stats().unwrap().physical_mem as f32 / 1_000_000.);

//     return new_map;
// }

pub fn twin_reads_from_snpmers(kmer_info: &mut KmerGlobalInfo, args: &Cli) -> Vec<TwinRead>{

    let start = std::time::Instant::now();

    let fastq_files = &kmer_info.read_files;
    let mut snpmer_set = HashSet::default();
    for snpmer_i in kmer_info.snpmer_info.iter(){
        let k = snpmer_i.k as usize;
        let snpmer1 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[0] as u64) << (k-1) );
        let snpmer2 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[1] as u64) << (k-1) );
        snpmer_set.insert(snpmer1);
        snpmer_set.insert(snpmer2);
    }

    let snpmer_set = Arc::new(snpmer_set);
    let twin_read_vec = Arc::new(Mutex::new(vec![]));

    let files_owned = fastq_files.clone();
    let hpc = args.homopolymer_compression;
    let solid_kmers_take = std::mem::take(&mut kmer_info.solid_kmers);
    let arc_solid = Arc::new(solid_kmers_take);

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
                if seq.len() < MIN_READ_LENGTH{
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
            let solid = Arc::clone(&arc_solid);
            let twrv = Arc::clone(&twin_read_vec);
            handles.push(thread::spawn(move || {
                loop{
                    match rx.recv() {
                        Ok(msg) => {
                            let seq = msg.0;
                            let seqlen = seq.len();
                            let qualities = msg.1;
                            let id = msg.2;
                            let twin_read = seeding::get_twin_read_syncmer(seq, qualities, k, c, set.as_ref(), id);
                            if twin_read.is_some(){

                                let mut solid_mini_positions = FxHashSet::default();
                                for (i, mini) in twin_read.as_ref().unwrap().minimizer_kmers.iter().enumerate(){
                                    if solid.contains(&mini){
                                        solid_mini_positions.insert(i);
                                    }
                                }
                                //< 5 % of the k-mers are solid; remove. This is usually due to highly repetitive stuff. 
                                if solid_mini_positions.len() < seqlen / c / 20{
                                    continue;
                                }

                                let mut solid_snpmer_positions = FxHashSet::default();
                                for (i, snpmer) in twin_read.as_ref().unwrap().snpmer_kmers.iter().enumerate(){
                                    if solid.contains(&snpmer){
                                        solid_snpmer_positions.insert(i);
                                    }
                                }

                                let mut twin_read = twin_read.unwrap();
                                twin_read.retain_mini_positions(solid_mini_positions);
                                twin_read.retain_snpmer_positions(solid_snpmer_positions);

                                //MinHash top ~ 1/20 * read_length of solid snpmers 
                                minhash_top_snpmers(&mut twin_read, MAX_FRACTION_OF_SNPMERS_IN_READ);

                                let mut vec = twrv.lock().unwrap();
                                vec.push(twin_read);
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

    kmer_info.solid_kmers = Arc::try_unwrap(arc_solid).unwrap();
    let mut twin_reads = Arc::try_unwrap(twin_read_vec).unwrap().into_inner().unwrap();
    twin_reads.sort_by(|a,b| a.id.cmp(&b.id));

    if log::log_enabled!(log::Level::Trace) {
        for twin_read in twin_reads.iter(){
            let decoded_snpmers = twin_read.snpmers().map(|x| (x.0, decode_kmer(x.1, twin_read.k))).collect::<Vec<_>>();
            log::trace!("{} {:?}", twin_read.id, decoded_snpmers);
        }
    }

    let number_reads_below_threshold = twin_reads.iter().filter(|x| x.est_id.is_some() && x.est_id.unwrap() < args.quality_value_cutoff).count();
    log::info!("Number of valid reads with >= 1kb - {}. Number of reads below quality threshold - {}.", twin_reads.len(), number_reads_below_threshold);
    let snpmer_densities = twin_reads.iter().map(|x| x.snpmer_kmers.len() as f64 / x.base_length as f64).collect::<Vec<_>>();
    let mean_snpmer_density = snpmer_densities.iter().sum::<f64>() / snpmer_densities.len() as f64;
    log::info!("Mean SNPmer density: {:.2}%", mean_snpmer_density * 100.);

    log::info!("Time elapsed for obtaining twin reads is: {:?}", start.elapsed());

    return twin_reads;
}

fn minhash_top_snpmers(twin_read: &mut TwinRead, max_fraction: f64){
    let top = (twin_read.base_length as f64 * max_fraction).ceil() as usize;
    if twin_read.snpmer_kmers.len() < top{
        return;
    }

    log::debug!("MinHashing read {}. Top {} snpmers for read of length {} with {} snpmers", &twin_read.id, top, twin_read.base_length, twin_read.snpmer_kmers.len());

    let split_mask = !(3 << (twin_read.k - 1));

    let mut splitmers_hash = twin_read.snpmer_kmers.iter().map(|x| mm_hash_64(x & split_mask)).collect::<Vec<_>>();
    splitmers_hash.sort();

    if top >= splitmers_hash.len(){
        return;
    }

    let hash_cutoff = splitmers_hash[top];
    let retain_positions = twin_read.snpmer_kmers.iter().enumerate().filter(|x| splitmers_hash[x.0] <= hash_cutoff).map(|x| x.0).collect::<FxHashSet<_>>();

    twin_read.snpmer_kmers = twin_read.snpmer_kmers.iter().enumerate().filter(|x| retain_positions.contains(&x.0)).map(|x| *x.1).collect::<Vec<_>>();
    twin_read.snpmer_positions = twin_read.snpmer_positions.iter().enumerate().filter(|x| retain_positions.contains(&x.0)).map(|x| *x.1).collect::<Vec<_>>();
}


#[inline]
pub fn split_kmer(kmer: u64, k: usize) -> (u64, u8){
    let split_mask_extract = 3 << (k-1);
    let mid_base = (kmer & split_mask_extract) >> (k-1);
    let masked_kmer = kmer & !split_mask_extract;
    return (masked_kmer, mid_base as u8);
}

pub fn get_snpmers(big_kmer_map: Vec<(Kmer64, [u32;2])>, k: usize, args: &Cli) -> KmerGlobalInfo{

    log::info!("Number of k-mers passing thresholds: {}", big_kmer_map.len());
    let mut new_map_counts_bases : FxHashMap<Kmer64, CountsAndBases> = FxHashMap::default();
    let mut kmer_counts = vec![];
    let mut solid_kmers = HashSet::default();
    let paths_to_files = args.input_files.iter().map(|x| std::fs::canonicalize(Path::new(x).to_path_buf()).unwrap()).collect::<Vec<_>>();

    for pair in big_kmer_map.iter(){
        let counts = pair.1;
        kmer_counts.push(counts[0] + counts[1]);
    }

    kmer_counts.par_sort_unstable();
    if kmer_counts.len() == 0{
        log::error!("No k-mers found. Exiting.");
        std::process::exit(1);
    }
    let high_freq_thresh = kmer_counts[kmer_counts.len() - (kmer_counts.len() / 50000) - 1].max(100);
    log::info!("High frequency k-mer threshold: {}", high_freq_thresh);
    drop(kmer_counts);

    //Should be able to parallelize this, TODO
    for pair in big_kmer_map.into_iter(){
        let kmer = pair.0;
        let (split_kmer, mid_base) = split_kmer(kmer, k);
        let counts = pair.1;
        if counts[0] > 0 && counts[1] > 0{
            let count = counts[0] + counts[1];
            if count < high_freq_thresh{
                solid_kmers.insert(kmer);
                let v = new_map_counts_bases.entry(split_kmer).or_insert(CountsAndBases{counts: SmallVec::new(), bases: SmallVec::new()});
                v.counts.push(counts);
                v.bases.push(mid_base);

            }
        }
    }

    if args.no_snpmers{
        log::info!("Skipping snpmer detection.");
        return KmerGlobalInfo{
            snpmer_info: vec![],
            solid_kmers: solid_kmers,
            high_freq_thresh: high_freq_thresh as f64,
            read_files: paths_to_files
        };
    }

    log::info!("Finding snpmers...");
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
            let right_p_val_thresh1 = twin_graph::binomial_test(n as u64, succ as u64, 0.025);
            let right_p_val_thresh2 = twin_graph::binomial_test(n as u64, succ as u64, 0.050);
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
    solid_kmers.shrink_to_fit();
    log::info!("Number of snpmers: {}. ", potential_snps.into_inner().unwrap());
    log::info!("Number of solid k-mers: {}.", solid_kmers.len());
    return KmerGlobalInfo{
        snpmer_info: snpmers,
        solid_kmers: solid_kmers,
        high_freq_thresh: high_freq_thresh as f64,
        read_files: paths_to_files
    };
}

pub fn parse_unitigs_into_table(cuttlefish_file: &str) -> (FxHashMap<u64, u32>, Vec<Vec<u8>>) {
    let mut kmer_to_unitig_count: FxHashMap<u64, u32> = FxHashMap::default();
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
