use dashmap::DashMap;
use crate::cli::Cli;
use std::path::Path;
use rayon::prelude::*;
use fxhash::FxHashMap;
use std::collections::VecDeque;
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
    fastq_file: &str,
    k: usize,
    threads: usize,
    homopolymer_comp: bool
) -> DashMap<(u64, u8), [u32;2]>{

    let (mut tx, rx) = spmc::channel();

    let file = fastq_file.to_owned();
    thread::spawn(move || {
        let mut reader = needletail::parse_fastx_file(file).expect("valid path");
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

    let map = Arc::new(DashMap::new());
    let mut handles = Vec::new();
    for _ in 0..threads {
        let map = Arc::clone(&map);
        let rx = rx.clone();
        handles.push(thread::spawn(move || {
            loop{
                match rx.recv() {
                    Ok(msg) => {
                        let split_kmer_info = seeding::split_kmer_mid(msg, k);
                        for kmer_i in split_kmer_info.iter() {
                            let mut c = map.entry((kmer_i.split_kmer, kmer_i.mid_base)).or_insert([0,0]);
                            if kmer_i.canonical{
                                c[0] += 1;
                            } else {
                                c[1] += 1;
                            }
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

   return Arc::try_unwrap(map).unwrap();
}

pub fn twin_reads_from_snpmers(snpmer_vec: &Vec<SnpmerInfo>, fastq_files: &[String], args: &Cli) -> Vec<TwinRead>{

    let mut snpmer_set = FxHashSet::default();
    for snpmer_i in snpmer_vec.iter(){
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

    let number_reads_below_threshold = twin_reads.iter().filter(|x| x.est_id.is_some() && x.est_id.unwrap() < args.quality_value_cutoff).count();
    log::info!("Number of valid reads - {}. Number of reads below quality threshold - {}.", twin_reads.len(), number_reads_below_threshold);

    return twin_reads;
}

pub fn get_snpmers(big_kmer_map: &DashMap<(u64, u8), [u32;2]>, k: usize) -> Vec<SnpmerInfo>{

    let mut new_map = FxHashMap::default();
    for pair in big_kmer_map.iter(){
        let kmer = *pair.key();
        let counts = *pair.value();
        let count = counts[0] + counts[1];
        if count > 2{
            let v = new_map.entry(kmer.0).or_insert(vec![]);
            v.push((counts, kmer.1));
        }
    }

    let mut potential_snps = 0;
    let mut snpmers = vec![];
    for (kmer, vec) in new_map.iter_mut(){
        if vec.len() > 1{
            vec.sort_unstable_by(|a, b| (b.0[0] + b.0[1]).cmp(&(a.0[0] + a.0[1])));
            //Add pseudocount... especially when all reads are 
            //already in forward strand (happens when debugging or manipulation via samtools)
            let a = vec[0].0[0] + 1;
            let b = vec[1].0[0] + 1;
            let c = vec[0].0[1] + 1;
            let d = vec[1].0[1] + 1;
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
                continue;
            }

            //Is snpmer
            if p_value > 0.005 || (odds < 1.5 && odds > 1./1.5){
                let snpmer = SnpmerInfo{
                    split_kmer: *kmer,
                    mid_bases: vec![vec[0].1, vec[1].1],
                    counts: vec![vec[0].0[0] + vec[0].0[1], vec[1].0[0] + vec[1].0[1]],
                    k: k as u8,
                };
                snpmers.push(snpmer);

                let snpmer1 = *kmer as u64 | ((vec[0].1 as u64) << k);
                let snpmer2 = *kmer as u64 | ((vec[1].1 as u64) << k);
                log::trace!("{} c:{:?} {} c:{:?}, p:{}, odds:{}", decode_kmer(snpmer1, k as u8), vec[0].0, decode_kmer(snpmer2, k as u8), vec[1].0, p_value, odds);
                potential_snps += 1;

            }
            else{
                log::trace!("NOT SNPMER c:{:?} c:{:?}, p:{}, odds:{}",  vec[0].0, vec[1].0, p_value, odds);
            }
        }
    }

    log::info!("Number of snpmers: {}. ", potential_snps);
    return snpmers;
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

pub fn kmer_chimera_filter(twin_reads: Vec<TwinRead>, kmer_counts: &DashMap<(u64, u8), [u32;2]>, k: u64) -> Vec<TwinRead>{

    let window = 5;
    let gap = 3;
    let mut solid_mini_to_reads = FxHashMap::default();
    let start = std::time::Instant::now();
    for (i,twin_read) in twin_reads.iter().enumerate(){
        let minimizers = &twin_read.minimizers;
        for &(pos, kmer) in minimizers.iter(){
            let mid_base = extract_mid_base(kmer, k);
            let masked_kmer = kmer & !(3 << (k-1));
            let counts = kmer_counts.get(&(masked_kmer, mid_base)).unwrap();
            if counts[0] + counts[1] > 2 {
                let v = solid_mini_to_reads.entry(kmer).or_insert(FxHashMap::default());
                if !v.contains_key(&i){
                    v.insert(i, [pos, pos]);
                }
                else{
                    v.get_mut(&i).unwrap()[1] = pos;
                }
            }
        }
    }
    log::debug!("Time elapsed for solid minimizer extraction: {:?}", start.elapsed());
    let good_twin_reads = Mutex::new(vec![]);
    twin_reads.into_par_iter().for_each(|twin_read| {
        let mut last_sets : VecDeque<FxHashSet<usize>> = VecDeque::new();
        let mut forward_sets: VecDeque<FxHashSet<usize>> = VecDeque::new();
        let mut gap_sets: VecDeque<FxHashSet<usize>> = VecDeque::new();
        let mut solid_forward : FxHashSet<usize> = FxHashSet::default();
        let mut solid_last : FxHashSet<usize> = FxHashSet::default();
        let mut moving_mean_last = 0.;
        let mut moving_mean_forward = 0.;
        let mut last_break = 0;
        let mut break_points = vec![];
        let mut counter = 0;
        let mut gap_positions = VecDeque::new();
        for &(pos,kmer) in twin_read.minimizers.iter(){
            if solid_mini_to_reads.contains_key(&kmer){
                let reads = solid_mini_to_reads.get(&kmer).unwrap();
                let current_set : FxHashSet<usize> = reads.keys().copied().collect();
                if last_sets.len() != window{
                    solid_last.extend(current_set.iter());
                    last_sets.push_back(current_set);
                    continue;
                }
                if gap_sets.len() != gap{
                    gap_positions.push_back(pos);
                    gap_sets.push_back(current_set);
                    continue;
                }
                if forward_sets.len() != window{
                    gap_positions.push_back(pos);
                    solid_forward.extend(current_set.iter());
                    forward_sets.push_back(current_set);
                    continue;
                }

                if moving_mean_forward == 0. && moving_mean_last == 0.{
                    moving_mean_forward = forward_sets.iter().map(|x| x.len()).sum::<usize>() as f64 / window as f64;
                    moving_mean_last = last_sets.iter().map(|x| x.len()).sum::<usize>() as f64 / window as f64;
                }

                let high_cov = solid_forward.len() > 3 && solid_last.len() > 3;

                //If intersection has length == 1, break;
                let mut break_set = true;
                let mut intersect = 0;
                let larger_set;
                let smaller_set;
                if solid_forward.len() < solid_last.len(){
                    larger_set = &solid_last;
                    smaller_set = &solid_forward;
                }
                else{
                    larger_set = &solid_forward;
                    smaller_set = &solid_last;
                }
                for read in smaller_set.iter(){
                    if larger_set.contains(read){
                        intersect += 1;
                        if intersect > (1.max((moving_mean_forward/20.) as usize)).min(3){
                            break_set = false;
                            break;
                        }
                    }
                }

                if break_set && high_cov{
                    if counter - last_break > gap + 1{
                        break_points.push(gap_positions[0]);
                        last_break = counter;
                        log::trace!("{:?}, {:?}, {}, {}", &solid_forward, &solid_last, &intersect, 1.max((moving_mean_forward/15.) as usize).min(2));
                    }
                }

                //Update

                let fset = &current_set;
                let remove_fset = solid_forward.difference(&forward_sets[0]).copied().collect::<Vec<usize>>();
                for read in remove_fset{
                    solid_forward.remove(&read);
                }
                if (fset.len() as f64) < (moving_mean_forward * window  as f64 * 0.5 + 3.){
                    for read in fset.iter(){
                        solid_forward.insert(*read);
                    }
                }

                moving_mean_forward = moving_mean_forward - (forward_sets[0].len() as f64) / window as f64 + (fset.len() as f64) / window as f64;

                let lset = &gap_sets[0];
                let remove_lset = solid_last.difference(&last_sets[0]).copied().collect::<Vec<usize>>();

                for read in remove_lset{
                    solid_last.remove(&read);
                }

                if (lset.len() as f64) < moving_mean_last * window as f64 * 0.5 + 3.{
                    for read in lset.iter(){
                        solid_last.insert(*read);
                    }
                }
                
                moving_mean_last = moving_mean_last - (last_sets[0].len() as f64) / window as f64 + (lset.len() as f64) / window as f64;

                let gap_front = gap_sets.pop_front().unwrap();
                last_sets.push_back(gap_front);
                last_sets.pop_front();
                let forward_front = forward_sets.pop_front().unwrap();
                gap_sets.push_back(forward_front);
                forward_sets.push_back(current_set);
                gap_positions.push_back(pos);
                gap_positions.pop_front();
                counter += 1;
            }
        }
        if break_points.len() > 0{
            log::trace!("Chimera FOUND: {}", twin_read.id);
            let new_reads = split_read(twin_read, break_points);
            for new_read in new_reads{
                good_twin_reads.lock().unwrap().push(new_read);
            }
        }
        else{
            good_twin_reads.lock().unwrap().push(twin_read);
        }
    });
    let mut good_twin_reads = good_twin_reads.into_inner().unwrap();
    good_twin_reads.sort_by(|a,b| a.id.cmp(&b.id));
    return good_twin_reads;
}

fn extract_mid_base(kmer: u64, k: u64) -> u8{
    let split_mask_extract = 3 << (k-1);
    let mid_base = (kmer & split_mask_extract) >> k-1;
    return mid_base as u8;
}


fn split_read(twin_read: TwinRead, mut break_points: Vec<usize>) -> Vec<TwinRead>{
    let mut new_reads = vec![];
    break_points.push(twin_read.base_length);
    let mut last_break = 0;
    for (i,break_point) in break_points.iter().enumerate(){
        if break_point - last_break > 1000{
            let mut new_read = TwinRead::default();
            new_read.minimizers = twin_read.minimizers.iter().filter(|x| x.0 >= last_break && x.0 < *break_point).copied().map(|x| (x.0 - last_break, x.1)).collect();
            new_read.snpmers = twin_read.snpmers.iter().filter(|x| x.0 >= last_break && x.0 < *break_point).copied().map(|x| (x.0 - last_break, x.1)).collect();
            new_read.id = format!("{}+split{}", &twin_read.id, i);
            log::trace!("Split read {} at {}-{}", &new_read.id, last_break, *break_point);
            new_read.k = twin_read.k;
            new_read.base_length = *break_point - last_break;
            new_read.dna_seq = twin_read.dna_seq[last_break..*break_point].to_owned();
            new_reads.push(new_read);
        }
        last_break = *break_point;
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
            let depths = map_info.mapping_boundaries().depth().collect::<Vec<_>>();
            for depth in depths{
                let mut string = format!("{} {}-{} COV:{}, BREAKPOINTS:", twin_read.id, depth.start, depth.stop, depth.val);
                for breakpoint in breakpoints.iter(){
                    string.push_str(format!("--{}--", breakpoint.pos).as_str());
                }
                writeln!(writer, "{}", &string).unwrap();
            }
            if breakpoints.len() == 0{
                new_outer_indices.push(new_twin_reads.len());
                new_twin_reads.push(twin_read);
                continue;
            }
            else{
                let splitted_reads = split_read(twin_read, breakpoints.iter().map(|x| x.pos).collect());
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