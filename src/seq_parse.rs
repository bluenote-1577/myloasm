use crate::cli::Cli;
use crate::kmc_reader::kmc::KmcReader;
use crate::types::BYTE_TO_SEQ;
use crate::types::Kmer48;
use fastbloom::BloomFilter;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;
use fxhash::FxHashMap;
use crate::seeding;
use crate::utils::*;
use std::io::BufReader;
use crossbeam_channel::bounded;

pub fn read_to_split_kmers(
    k: usize,
    threads: usize,
    args: &Cli,
) -> Vec<(u64,[u32;2])>{
 

    let start = std::time::Instant::now();
    let bf_vec_maps = first_iteration(k, threads, args);
    log::info!("Finished with bloom filter processing in {:?}. Round 2 - Start k-mer counting...", start.elapsed());
    log_memory_usage(true, "Memory usage after bloom filter processing");

    let start = std::time::Instant::now();
    let vec_maps = second_iteration(k, threads, args, bf_vec_maps); 
    let map_size_raw = vec_maps.iter().map(|x| x.len()).sum::<usize>();
    log::info!("Finished with k-mer counting in {:?}. Total kmers after bloom filter: {}", start.elapsed(), map_size_raw);
    log_memory_usage(true, "Memory usage after second round of k-mer counting processing");

    let mut vectorized_map = vec![];
    for map in vec_maps.into_iter(){
        for (kmer, counts) in map.into_iter(){
            if counts[0] > 0 && counts[1] > 0 && counts[0] + counts[1] > 2{
                vectorized_map.push((kmer, counts));
            }
        }
    }

    let map_size_retain = vectorized_map.len();
    log::info!("Removed {} kmers with counts < 1 in both strands and <= 3 multiplicity.", map_size_raw - map_size_retain);
    if map_size_retain < map_size_raw / 1000 {
        log::warn!("Less than 0.1% of kmers have counts > 1 in both strands and > 2 multiplicity. This may indicate a problem with the input data or very low coverage.");
    }
    log::debug!("Final Hashmap len after vectorization: {}", map_size_retain);
    log_memory_usage(false, "Memory usage after second round of k-mer counting processing");

    return vectorized_map;
}

fn estimate_bf_size(
    args: &Cli,
) -> f64{
    let mut est_bases = 0.;
    let mut is_gzipped = false;
    for fq_file in args.input_files.iter(){
        if fq_file.contains(".gz") || fq_file.contains(".gzip") || fq_file.contains(".bz") {
            is_gzipped = true;
        }
        let metadata = std::fs::metadata(fq_file).expect("Unable to read file metadata");
        log::debug!("File: {}, size (Gbytes): {}", fq_file, metadata.len() as f64 / 1_000_000_000.);
        est_bases += metadata.len() as f64 / 1_000_000_000.;
        if is_gzipped{
            est_bases *= 1.5; //rough estimate of compression ratio
        }
        else{
            est_bases /= 2.;
        }
    }
    let bf_size = (est_bases / 2.0).min(200.0).max(2.0);
    return bf_size;
}

fn first_iteration(
    k: usize,
    threads: usize,
    args: &Cli,
) -> Vec<FxHashMap<u64,[u32;2]>>
{
    //Topology is
    //      A-SEND: tx_head, , B-REC: rx_head1, rx_head2...
    // |   |  ... 
    // B   B  ...  B-SEND: txs[0...], txs2[0...],... C-REC: rxs
    // | x | x | ...
    // C   C  ...
    let hm_size = threads;
    let mask = !(1 << 63);
    let bf_size;
    if let Some(bf_size_manual) = args.bloom_filter_size {
        bf_size = bf_size_manual;
    }
     else {
        bf_size = estimate_bf_size(args);
        log::info!("Using automatic bloom filter size: {:.2} Gbytes", bf_size);
    }

    let aggressive_bloom = args.aggressive_bloom;
    let mut bf_vec_maps : Vec<FxHashMap<u64, [u32;2]>> = vec![FxHashMap::default(); hm_size];
    if bf_size > 0.{
        let num_b = threads/10 + 1;
        let counter = Arc::new(Mutex::new(0));
        let read_lengths = Arc::new(Mutex::new(vec![]));
        let mut rxs = vec![];
        let mut txs_vecs = vec![vec![]; num_b];
        for _ in 0..threads {
            //let (tx, rx) = unbounded();
            let (tx, rx) = bounded(500);
            for i in 1..num_b{
                txs_vecs[i].push(tx.clone());
            }
            txs_vecs[0].push(tx);
            rxs.push(rx);
        }

        //let (tx_head, rx_head1) = unbounded();
        let (tx_head, rx_head1) = bounded(500);
        let mut rx_heads = vec![];
        for _ in 1..num_b{
            let rx_head2 = rx_head1.clone();
            rx_heads.push(rx_head2);
        }
        rx_heads.push(rx_head1);

        assert!(txs_vecs.len() == rx_heads.len());

        let fq_files = args.input_files.clone();
        //A: Get k-mers
        thread::spawn(move || {
            for fq_file in fq_files{
                let bufreader = BufReader::new(std::fs::File::open(fq_file).expect("valid path"));
                let mut reader = needletail::parse_fastx_reader(bufreader).expect("valid path");
                while let Some(record) = reader.next() {
                    let rec = record.expect("Error reading record");
                    let seq = rec.seq().to_vec();
                    let qualities = rec.qual().map(Vec::from);
                    tx_head.send((seq,qualities)).unwrap();
                }
            }
            drop(tx_head);
            log::debug!("Finished reading all reads.");
        });
        //B: Process kmers and send to hash maps
        for (rx_head, txs) in rx_heads.into_iter().zip(txs_vecs.into_iter()){
            let clone_counter = Arc::clone(&counter);
            let clone_read_lengths = Arc::clone(&read_lengths);
            thread::spawn(move || {
                loop{
                    match rx_head.recv() {
                        Ok((seq, qualities)) => {

                            {
                                let mut read_lengths = clone_read_lengths.lock().unwrap();
                                read_lengths.push(seq.len());
                            }
                            let split_kmer_info = seeding::split_kmer_mid(seq, qualities, k);
                            let mut vec_and_canon = vec![vec![]; hm_size];
                            for kmer_i_and_canon in split_kmer_info.into_iter() {
                                let kmer = kmer_i_and_canon & mask;
                                let hash = kmer % threads as u64;
                                vec_and_canon[hash as usize].push(kmer_i_and_canon);
                            }

                            for (i, vec) in vec_and_canon.into_iter().enumerate(){
                                txs[i].send(vec).unwrap();
                            }

                            {
                                let mut counter = clone_counter.lock().unwrap();
                                *counter += 1;
                                if *counter % 100000 == 0{
                                    log::info!("Processed {} reads.", counter);
                                }
                                if *counter % 1_000_000 == 0{
                                    log_memory_usage(false, &format!("Processed {} reads for bloom filter stage", *counter));
                                }
                            }
                        }
                        Err(_) => {
                            break;
                        }
                    }
                }
                for tx in txs{
                    drop(tx);
                }
            });
        }

        //C: Update bloom filter
        let mut handles = Vec::new();
        for rx in rxs.into_iter() {
            handles.push(thread::spawn(move || {
                let mut filter_canonical = BloomFilter::with_num_bits((bf_size * 4. * 1_000_000_000. / threads as f64) as usize).seed(&42).expected_items((bf_size * 4. * 1_000_000_000. / 10. / threads as f64) as usize);
                let mut filter_noncanonical = BloomFilter::with_num_bits((bf_size * 4. * 1_000_000_000. / threads as f64) as usize).seed(&42).expected_items((bf_size * 4. * 1_000_000_000. / 10. / threads as f64) as usize);
                let mut map: FxHashMap<u64,[u32;2]> = FxHashMap::default();
                loop{
                    match rx.recv() {
                        Ok(msg) => {
                            let kmer_vecs = msg;
                            for kmer_i_canon in kmer_vecs{
                                let canonical = kmer_i_canon >> 63;
                                let kmer = kmer_i_canon & mask;
                                let kmer_canon = kmer | (1 << 63);
                                if canonical == 1{
                                    let already_present_canon = filter_canonical.insert(&kmer_canon);
                                    let already_present_noncanon = filter_noncanonical.contains(&kmer);
                                    if aggressive_bloom{
                                        if already_present_noncanon && already_present_canon{
                                            map.insert(kmer, [0,0]);
                                        }
                                    }
                                    else{
                                        if already_present_noncanon{
                                            map.insert(kmer, [0,0]);
                                        }
                                    }
                                }
                                else{
                                    let already_present_noncanon = filter_noncanonical.insert(&kmer);
                                    let already_present_canon = filter_canonical.contains(&kmer_canon);
                                    if aggressive_bloom{
                                        if already_present_noncanon && already_present_canon{
                                            map.insert(kmer, [0,0]);
                                        }
                                    }
                                    else{
                                        if already_present_canon{
                                            map.insert(kmer, [0,0]);
                                        }
                                    }
                                }
                            }
                        }
                        Err(_) => {
                            log::trace!("Thread finished.");
                            break;
                        }
                    }
                }
                
                
                map.shrink_to_fit();
                map
            }));
        }

        for (map_ind, handle) in handles.into_iter().enumerate() {
            bf_vec_maps[map_ind] = handle.join().unwrap();
            bf_vec_maps[map_ind].shrink_to_fit();
        };

        let read_lengths = Arc::try_unwrap(read_lengths).unwrap().into_inner().unwrap();
        log::info!("Read lengths - {}", get_nx_from_vec(&read_lengths, &[10,50,90]));
        log::info!("Total bases - {} million bp", (read_lengths.iter().map(|x| *x as usize).sum::<usize>() as f64 / 1_000_000.).round());
    
    }  
    
    return bf_vec_maps;
}

fn second_iteration(
    k: usize,
    threads: usize,
    args: &Cli,
    bf_vec_maps: Vec<FxHashMap<u64, [u32;2]>>) -> Vec<FxHashMap<u64, [u32;2]>>{

    let bf_size = args.bloom_filter_size;
    let mask = !(1 << 63);
    let mut vec_maps : Vec<FxHashMap<u64, [u32;2]>> = vec![FxHashMap::default(); threads];

    let num_b = threads/10 + 1;
    let counter = Arc::new(Mutex::new(0));
    let mut rxs = vec![];
    let mut txs_vecs = vec![vec![]; num_b];
    for _ in 0..threads {
        //let (tx, rx) = unbounded();
        let (tx, rx) = bounded(500);
        for i in 1..num_b{
            txs_vecs[i].push(tx.clone());
        }
        txs_vecs[0].push(tx);
        rxs.push(rx);
    }

    //let (tx_head, rx_head1) = unbounded();
    let (tx_head, rx_head1) = bounded(500);
    let mut rx_heads = vec![];
    for _ in 1..num_b{
        let rx_head2 = rx_head1.clone();
        rx_heads.push(rx_head2);
    }
    rx_heads.push(rx_head1);
    let bf_size = if let Some(bf_size_manual) = bf_size {
        bf_size_manual
    } else {
        estimate_bf_size(args)
    };

    assert!(txs_vecs.len() == rx_heads.len());

    let fq_files = args.input_files.clone();
    thread::spawn(move || {
        for fq_file in fq_files{
            let bufreader = BufReader::new(std::fs::File::open(fq_file).expect("valid path"));
            let mut reader = needletail::parse_fastx_reader(bufreader).expect("valid path");
            while let Some(record) = reader.next() {
                let rec = record.expect("Error reading record");
                let seq = rec.seq().to_vec();
                let qualities = rec.qual().map(Vec::from);
                tx_head.send((seq,qualities)).unwrap();
            }
        }
        drop(tx_head);
        log::debug!("Finished reading all reads.");
    });

    //B: Process kmers and send to hash maps
    for (rx_head, txs) in rx_heads.into_iter().zip(txs_vecs.into_iter()){
        let clone_counter = Arc::clone(&counter);
        thread::spawn(move || {
            loop{
                match rx_head.recv() {
                    Ok((seq,qualities)) => {
                        let split_kmer_info = seeding::split_kmer_mid(seq, qualities, k);
                        let mut vec_and_canon = vec![vec![]; threads];
                        for kmer_i_and_canon in split_kmer_info.into_iter() {
                            let kmer = kmer_i_and_canon & mask;
                            let hash = kmer % threads as u64;
                            vec_and_canon[hash as usize].push(kmer_i_and_canon);
                        }

                        for (i, vec) in vec_and_canon.into_iter().enumerate(){
                            txs[i].send(vec).unwrap();
                        }

                        {
                            let mut counter = clone_counter.lock().unwrap();
                            *counter += 1;
                            if *counter % 100000 == 0{
                                log::debug!("Processed {} reads.", counter);
                            }
                        }
                    }
                    Err(_) => {
                        break;
                    }
                }
            }
            for tx in txs{
                drop(tx);
            }
        });
    }

    let mut handles = Vec::new();
    for (rx, my_map) in rxs.into_iter().zip(bf_vec_maps.into_iter()){
        handles.push(thread::spawn(move || {
            let mut my_map = my_map;
            loop{
                match rx.recv() {
                    Ok(msg) => {
                        let vec_and_canon = msg;
                        if bf_size > 0.{
                            for kmer_and_canon in vec_and_canon.into_iter(){
                                let kmer = kmer_and_canon & mask;
                                let canon = kmer_and_canon >> 63;
                                if let Some(val) = my_map.get_mut(&kmer){
                                    val[canon as usize] += 1;
                                }
                            }
                        }
                        else{
                            for kmer_and_canon in vec_and_canon.into_iter(){
                                let kmer = kmer_and_canon & mask;
                                let canon = kmer_and_canon >> 63;
                                let val = my_map.entry(kmer).or_insert([0,0]);
                                val[canon as usize] += 1
                            }
                        }
                    }
                    Err(_) => {
                        log::trace!("Thread finished.");
                        break;
                    }
                }
            }
            my_map.shrink_to_fit();
            my_map
        }));
    }

    for (i,handle) in handles.into_iter().enumerate() {
        vec_maps[i] = handle.join().unwrap();
    };

    return vec_maps;
}

// Intuitively I think this helps, not sure if we want to use it TODO
pub fn quality_pool(qualities: Vec<u8>) -> Vec<u8>{
    let pool_width = 5;
    let mut pool = Vec::new();
    for i in 0..qualities.len(){
        if i > pool_width/2 && i < qualities.len() - pool_width/2{
            let mut min = 255;
            for j in i-pool_width/2..i+pool_width/2{
                if qualities[j] < min{
                    min = qualities[j];
                }
            }
            pool.push(min);
        }
        else{
            pool.push(qualities[i]);
        }
    }
    return pool;
}

/// Read k-mers from a precomputed KMC database and convert to our internal format.
///
/// KMC stores canonical k-mers (lexicographically smallest between k-mer and its RC).
/// We convert these to our format which uses a different canonicalization based on
/// comparing split k-mers (k-mer with middle base masked out).
///
/// Since KMC doesn't provide strand information per occurrence, we split the count
/// between both strand slots to satisfy downstream filtering requirements.
pub fn read_kmers_from_kmc_db(
    k: usize,
    _threads: usize,
    kmc_db_path: &str,
    _args: &Cli,
) -> Vec<(u64, [u32; 2])> {
    log::info!("Reading k-mers from KMC database: {}", kmc_db_path);

    let reader = KmcReader::open(kmc_db_path).expect("Failed to open KMC database. Ensure the path is correct and the KMC database is valid.");
    let kmc_info = reader.info();

    if kmc_info.kmer_length as usize != k {
        log::error!(
            "KMC database k-mer length ({}) does not match expected k ({})",
            kmc_info.kmer_length, k
        );
        std::process::exit(1);
    }

    log::info!(
        "KMC database info: k={}, total_kmers={}, both_strands={}",
        kmc_info.kmer_length,
        kmc_info.total_kmers,
        kmc_info.both_strands
    );


    // Mask to remove middle base for split comparison
    let split_mask: u64 = !(3u64 << (k - 1));
    // Mask to keep only the k-mer bits
    let marker_mask: u64 = u64::MAX >> (64 - 2 * k);

    let mut kmer_maps: Vec<FxHashMap<Kmer48, [u16; 2]>> = vec![FxHashMap::default(); 50];
    let mut processed_count = 0u64;
    let mut skipped_low_count = 0u64;
    let mut skipped_palindrome = 0u64;


    for bloom in [true, false]{
        processed_count = 0;
        let mut bloom_filter;
        if bloom{
            bloom_filter = BloomFilter::with_num_bits(kmc_info.total_kmers as usize * 8).seed(&42).expected_items(kmc_info.total_kmers as usize);
        }
        else{
            bloom_filter = BloomFilter::with_num_bits(1).seed(&42).expected_items(1); //dummy
        }
        
        if bloom{
            log::info!("Starting bloom filter pass...");
        }
        else{
            log::info!("Starting k-mer counting pass...");
        }
        let mut reader = KmcReader::open(kmc_db_path).expect("Failed to open KMC database");
        while let Some(record) = reader.next_kmer().expect("Error reading k-mer") {
            let kmer_bytes = record.kmer;
            let kmc_count = record.count as u32;

            // Skip if count is too low (will be filtered anyway)
            if kmc_count < 1 {
                skipped_low_count += 1;
                continue;
            }

            // Encode forward k-mer from the KMC canonical form
            let mut forward_kmer: u64 = 0;
            for &byte in kmer_bytes {
                let nuc = BYTE_TO_SEQ[byte as usize] as u64;
                forward_kmer <<= 2;
                forward_kmer |= nuc;
            }
            forward_kmer &= marker_mask;

            // Compute reverse complement
            let mut reverse_kmer: u64 = 0;
            for &byte in kmer_bytes.iter().rev() {
                let nuc = 3 - BYTE_TO_SEQ[byte as usize] as u64;
                reverse_kmer <<= 2;
                reverse_kmer |= nuc;
            }
            reverse_kmer &= marker_mask;

            // Compute split forms (mask out middle base)
            let split_f = forward_kmer & split_mask;
            let split_r = reverse_kmer & split_mask;

            // Skip palindromes (same split form in both directions)
            if split_f == split_r {
                skipped_palindrome += 1;
                continue;
            }

            // Determine our canonical form based on split comparison
            // canonical_marker = true means forward is canonical (bit 63 = 1)
            let canonical = split_f < split_r;
            let canonical_kmer = if canonical {
                forward_kmer
            } else {
                reverse_kmer
            };

            if bloom{
                let already_present = bloom_filter.insert(&canonical_kmer);
                if !already_present{
                    continue;
                }
                else{
                    kmer_maps[canonical_kmer as usize % 50].insert(Kmer48::from_u64(canonical_kmer), [0,0]);
                }
            }

            else{
                if let Some(entry) = kmer_maps[canonical_kmer as usize % 50].get_mut(&Kmer48::from_u64(canonical_kmer)) {
                    let count;
                    if kmc_count > u16::MAX as u32 {
                        count = u16::MAX;
                    }
                    else{
                        count = kmc_count as u16;
                    }
                    if canonical {
                        entry[0] = count;
                    }
                    else{
                        entry[1] = count;
                    }
                }
                else{
                    // This k-mer was not seen in bloom filter pass, skip
                    continue;
                }
            }

            processed_count += 1;
            if processed_count % 10_000_000 == 0 {
                log::info!("Processed {} k-mers from KMC database", processed_count);
            }
        }
        log_memory_usage(true, if bloom{"Memory usage after KMC bloom filter pass"} else{"Memory usage after KMC k-mer counting pass"});
    }

    log::info!(
        "Finished reading KMC database: {} k-mers processed, {} skipped (low count), {} skipped (palindrome)",
        processed_count / 2,
        skipped_low_count / 2,
        skipped_palindrome / 2,
    );

    // Convert to the expected Vec format with filtering
    let mut vectorized_map = vec![];
    for kmer_map in kmer_maps.into_iter() {
        for (kmer, counts) in kmer_map.into_iter() {
            if counts[0] > 0 && counts[1] > 0 && counts[0] + counts[1] > 2 {
                vectorized_map.push((kmer.to_u64(), counts.map(|c| c as u32)));
            }
        }
    }

    log::info!(
        "Retained {} k-mers after filtering (both strands > 0, total > 2)",
        vectorized_map.len()
    );
    log_memory_usage(true, "Memory usage after reading KMC database");

    vectorized_map
}