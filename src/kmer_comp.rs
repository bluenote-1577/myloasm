use dashmap::DashMap;
use fxhash::FxHashMap;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;
use crate::types::*;
use crate::seeding;
use fishers_exact::fishers_exact;
use fxhash::FxHashSet;

pub fn read_to_split_kmers(
    fastq_file: &str,
    k: usize,
    threads: usize
) -> DashMap<(u64, u8), [u32;2]>{

    let (mut tx, rx) = spmc::channel();

    let file = fastq_file.to_owned();
    thread::spawn(move || {
        let mut reader = needletail::parse_fastx_file(file).expect("valid path");
        while let Some(record) = reader.next() {
            let rec = record.expect("Error reading record");
            let seq = rec.seq().to_vec();
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

pub fn twin_reads_from_snpmers(snpmer_vec: &Vec<SnpmerInfo>, fastq_files: &[&str], k: usize, c: usize, threads: usize) -> Vec<TwinRead>{

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
    for fastq_file in files_owned{
        let (mut tx, rx) = spmc::channel();
        thread::spawn(move || {
            let mut reader = needletail::parse_fastx_file(fastq_file).expect("valid path");
            while let Some(record) = reader.next() {
                let rec = record.expect("Error reading record");
                let seq = rec.seq().to_vec();
                let id = String::from_utf8_lossy(rec.id()).to_string();
                tx.send((seq, id)).unwrap();
            }
        });

        let mut handles = Vec::new();
        for _ in 0..threads {
            let rx = rx.clone();
            let set = Arc::clone(&snpmer_set);
            let twrv = Arc::clone(&twin_read_vec);
            handles.push(thread::spawn(move || {
                loop{
                    match rx.recv() {
                        Ok(msg) => {
                            let seq = msg.0;
                            let id = msg.1;
                            let twin_read = seeding::get_twin_read(seq, k, c, set.as_ref(), id);
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

    let twin_reads = Arc::try_unwrap(twin_read_vec).unwrap().into_inner().unwrap();

    for twin_read in twin_reads.iter(){
        log::trace!("{}\t{}\t{}",twin_read.id, twin_read.minimizers.len(), twin_read.snpmers.len());
        for snpmer in twin_read.snpmers.iter(){
            log::trace!("{}, {}", decode_kmer(snpmer.1, k as u8), snpmer.0);
        }
    }
    return twin_reads;
}

pub fn get_snpmers(big_kmer_map: DashMap<(u64, u8), [u32;2]>, k: usize) -> Vec<SnpmerInfo>{

    dbg!(big_kmer_map.len());

    let mut new_map = FxHashMap::default();
    for (kmer, counts) in big_kmer_map{
        let count = counts[0] + counts[1];
        if count > 2{
            let v = new_map.entry(kmer.0).or_insert(vec![]);
            v.push((counts, kmer.1));
        }
    }

    dbg!(new_map.len());

    let mut potential_snps = 0;
    let mut snpmers = vec![];
    for (kmer, vec) in new_map.iter_mut(){
        if vec.len() > 1{
            vec.sort_unstable_by(|a, b| (b.0[0] + b.0[1]).cmp(&(a.0[0] + a.0[1])));
            let contingency_table = [
                vec[0].0[0], vec[0].0[1],
                vec[1].0[0], vec[1].0[1],
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
                println!("{} c:{:?} {} c:{:?}, p:{}, odds:{}", decode_kmer(snpmer1, k as u8 + 1), vec[0].0, decode_kmer(snpmer2, k as u8 + 1), vec[1].0, p_value, odds);
                potential_snps += 1;

            }
        }
    }

    dbg!(potential_snps);
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



pub fn read_to_unitig_count_vector(
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

