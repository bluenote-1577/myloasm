use crate::cli::Cli;
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
    let bf_size = args.bloom_filter_size;
    let mut bf_vec_maps : Vec<FxHashMap<u64, [u32;2]>> = vec![FxHashMap::default(); hm_size];
    if bf_size > 0.{
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
            thread::spawn(move || {
                loop{
                    match rx_head.recv() {
                        Ok((seq, qualities)) => {
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
        for rx in rxs.into_iter(){
            handles.push(thread::spawn(move || {
                let mut filter_canonical = BloomFilter::with_num_bits((bf_size * 4. * 1_000_000_000. / threads as f64) as usize).seed(&42).expected_items(25_000_000_000 / threads);
                let mut filter_noncanonical = BloomFilter::with_num_bits((bf_size * 4. * 1_000_000_000. / threads as f64) as usize).seed(&42).expected_items(25_000_000_000 / threads);
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
                                    filter_canonical.insert(&kmer_canon);
                                    let already_present_noncanon = filter_noncanonical.contains(&kmer);
                                    if already_present_noncanon{
                                        map.insert(kmer, [0,0]);
                                    }
                                }
                                else{
                                    filter_noncanonical.insert(&kmer);
                                    let already_present_canon = filter_canonical.contains(&kmer_canon);
                                    if already_present_canon{
                                        map.insert(kmer, [0,0]);
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