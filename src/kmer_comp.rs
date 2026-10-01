use crate::cli::Cli;
use crate::constants::HUFFMAN_SAMPLE_MAX;
use crate::constants::MAX_KMER_COUNT_IN_READ;
use crate::constants::MIN_READ_LENGTH;
use crate::constants::USE_SOLID_KMERS;
use crate::seeding;
use crate::twin_graph;
use crate::types::*;
use crate::utils;
use fishers_exact::fishers_exact;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use rayon::prelude::*;
use smallvec::SmallVec;
use std::collections::HashSet;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;

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

type TwinReadInput = (Vec<u8>, Option<Vec<u8>>, String);

fn spawn_twin_read_workers(
    rx: spmc::Receiver<TwinReadInput>,
    num_threads: usize,
    k: usize,
    c: usize,
    snpmer_set: Arc<FxHashSet<Kmer64>>,
    solid_kmers: Arc<HashSet<Kmer48>>,
    high_freq_kmers: Arc<HashSet<Kmer48>>,
    twin_reads: Arc<Mutex<Vec<TwinRead>>>,
    num_reads_removed_repetitive: Arc<Mutex<usize>>,
) -> Vec<thread::JoinHandle<()>> {
    let mut handles = Vec::with_capacity(num_threads);
    for _ in 0..num_threads {
        let rx = rx.clone();
        let set = Arc::clone(&snpmer_set);
        let solid = Arc::clone(&solid_kmers);
        let highfreq = Arc::clone(&high_freq_kmers);
        let twrv = Arc::clone(&twin_reads);
        let num_repetitive = Arc::clone(&num_reads_removed_repetitive);
        handles.push(thread::spawn(move || {
            while let Ok((seq, qualities, id)) = rx.recv() {
                let seqlen = seq.len();
                let twin_read = seeding::get_twin_read_syncmer(seq, qualities, k, c, &set, id);
                if let Some(mut twin_read) = twin_read {
                    let mut kmer_counter_map = FxHashMap::default();
                    for mini in twin_read.minimizer_kmers().iter() {
                        *kmer_counter_map.entry(*mini).or_insert(0) += 1;
                    }

                    let mut solid_mini_indices = FxHashSet::default();
                    for (i, mini) in twin_read.minimizer_kmers().iter().enumerate() {
                        if kmer_counter_map[mini] > MAX_KMER_COUNT_IN_READ {
                            continue;
                        }
                        if (USE_SOLID_KMERS && solid.contains(mini))
                            || (!USE_SOLID_KMERS && !highfreq.contains(mini))
                        {
                            solid_mini_indices.insert(i);
                        }
                    }

                    // < 5% of the k-mers are solid; remove. This is usually due to
                    // highly repetitive sequence.
                    if solid_mini_indices.len() < seqlen / c / 20 {
                        *num_repetitive.lock().unwrap() += 1;
                        continue;
                    }

                    let mut solid_snpmer_indices = FxHashSet::default();
                    for (i, snpmer) in twin_read.snpmer_kmers().iter().enumerate() {
                        if (USE_SOLID_KMERS && solid.contains(snpmer))
                            || (!USE_SOLID_KMERS && !highfreq.contains(snpmer))
                        {
                            solid_snpmer_indices.insert(i);
                        }
                    }

                    twin_read.retain_mini_indices(solid_mini_indices);
                    twin_read.retain_snpmer_indices(solid_snpmer_indices);
                    twrv.lock().unwrap().push(twin_read);
                }
            }
        }));
    }
    handles
}

pub fn twin_reads_from_snpmers(kmer_info: &mut KmerGlobalInfo, args: &Cli) -> Vec<TwinRead> {
    let start = std::time::Instant::now();

    let fastq_files = &kmer_info.read_files;
    let mut snpmer_set = FxHashSet::default();
    for snpmer_i in kmer_info.snpmer_info.iter() {
        let k = args.kmer_size as usize;
        let snpmer1 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[0] as u64) << (k - 1));
        let snpmer2 = snpmer_i.split_kmer as u64 | ((snpmer_i.mid_bases[1] as u64) << (k - 1));
        snpmer_set.insert(snpmer1);
        snpmer_set.insert(snpmer2);
    }

    let snpmer_set = Arc::new(snpmer_set);
    let twin_read_vec = Arc::new(Mutex::new(vec![]));

    let files_owned = fastq_files.clone();
    let solid_kmers_take = std::mem::take(&mut kmer_info.solid_kmers);
    let high_freq_kmers_take = std::mem::take(&mut kmer_info.high_freq_kmers);
    let arc_solid = Arc::new(solid_kmers_take);
    let arc_high_freq = Arc::new(high_freq_kmers_take);
    let num_reads_removed_repetitive = Arc::new(Mutex::new(0));

    // The input producer pauses after a deterministic prefix. All prefix work is
    // drained before the Huffman tables are built and atomically published; only
    // then may later reads be encoded. This is a real phase boundary rather than a
    // worker-scheduling-dependent counter.
    let (prefix_tx, prefix_rx) = spmc::channel::<TwinReadInput>();
    let (mut suffix_tx, suffix_rx) = spmc::channel::<TwinReadInput>();
    let (phase_tx, phase_rx) = std::sync::mpsc::channel::<bool>();
    let (resume_tx, resume_rx) = std::sync::mpsc::channel::<()>();
    let input_handle = thread::spawn(move || {
        let mut prefix_tx = Some(prefix_tx);
        let mut prefix_count = 0usize;
        let mut in_prefix = true;

        for fastq_file in files_owned {
            let mut reader = needletail::parse_fastx_file(fastq_file).expect("valid path");
            while let Some(record) = reader.next() {
                let rec = record.expect("Error reading record");
                let seq = rec.seq().to_vec();
                if seq.len() < MIN_READ_LENGTH {
                    continue;
                }
                let id = String::from_utf8_lossy(rec.id()).to_string();
                let qualities = rec.qual().map(|qualities| qualities.to_vec());
                let message = (seq, qualities, id);

                if in_prefix {
                    prefix_tx.as_mut().unwrap().send(message).unwrap();
                    prefix_count += 1;
                    if prefix_count == HUFFMAN_SAMPLE_MAX {
                        drop(prefix_tx.take());
                        phase_tx.send(true).unwrap();
                        resume_rx.recv().unwrap();
                        in_prefix = false;
                    }
                } else {
                    suffix_tx.send(message).unwrap();
                }
            }
        }

        if in_prefix {
            drop(prefix_tx.take());
            phase_tx.send(false).unwrap();
        }
    });

    let prefix_handles = spawn_twin_read_workers(
        prefix_rx,
        args.threads,
        args.kmer_size,
        args.c,
        Arc::clone(&snpmer_set),
        Arc::clone(&arc_solid),
        Arc::clone(&arc_high_freq),
        Arc::clone(&twin_read_vec),
        Arc::clone(&num_reads_removed_repetitive),
    );
    let reached_sample_limit = phase_rx.recv().unwrap();
    for handle in prefix_handles {
        handle.join().unwrap();
    }

    if !huffman_initialized() {
        let mut reads = twin_read_vec.lock().unwrap();
        if reached_sample_limit {
            log::info!(
                "Reached deterministic prefix of {} input reads.",
                HUFFMAN_SAMPLE_MAX
            );
        }
        log::info!("Building Huffman coding from {} reads...", reads.len());
        build_and_set_huffman_coding(&reads, reads.len());
        reencode_reads_huffman(&mut reads);
        log::info!("Huffman coding built and {} reads re-encoded.", reads.len());
    }

    if reached_sample_limit {
        let suffix_handles = spawn_twin_read_workers(
            suffix_rx,
            args.threads,
            args.kmer_size,
            args.c,
            Arc::clone(&snpmer_set),
            Arc::clone(&arc_solid),
            Arc::clone(&arc_high_freq),
            Arc::clone(&twin_read_vec),
            Arc::clone(&num_reads_removed_repetitive),
        );
        resume_tx.send(()).unwrap();
        for handle in suffix_handles {
            handle.join().unwrap();
        }
    }
    input_handle.join().unwrap();

    kmer_info.solid_kmers = Arc::try_unwrap(arc_solid).unwrap();
    kmer_info.high_freq_kmers = Arc::try_unwrap(arc_high_freq).unwrap();
    let mut twin_reads = Arc::try_unwrap(twin_read_vec)
        .unwrap()
        .into_inner()
        .unwrap();
    twin_reads.sort_by(|a, b| a.id.cmp(&b.id));
    twin_reads.par_iter_mut().for_each(|x| x.compact());

    //This is for hefty debugging purposes only; too verbose.
    if log::log_enabled!(log::Level::Trace) && false {
        for twin_read in twin_reads.iter() {
            let decoded_snpmers = twin_read
                .snpmers_vec()
                .into_iter()
                .map(|x| (x.0, decode_kmer48(x.1, twin_read.k)))
                .collect::<Vec<_>>();
            log::trace!("{} {:?}", twin_read.id, decoded_snpmers);
        }
    }

    log::debug!(
        "Number of reads filtered due to repetitiveness: {}",
        num_reads_removed_repetitive.lock().unwrap()
    );
    let number_reads_below_threshold = twin_reads
        .iter()
        .filter(|x| x.est_id.is_some() && x.est_id.unwrap() < args.quality_value_cutoff)
        .count();
    log::info!(
        "Number of valid reads with >= 1kb - {}. Number of reads below quality threshold - {}.",
        twin_reads.len(),
        number_reads_below_threshold
    );
    let snpmer_densities = twin_reads
        .iter()
        .map(|x| x.snpmer_count() as f64 / x.base_length as f64)
        .collect::<Vec<_>>();
    let mean_snpmer_density = snpmer_densities.iter().sum::<f64>() / snpmer_densities.len() as f64;
    log::info!("Mean SNPmer density: {:.2}%", mean_snpmer_density * 100.);
    log::info!(
        "Time elapsed for obtaining twin reads is: {:?}",
        start.elapsed()
    );

    return twin_reads;
}

#[inline]
pub fn retrieve_masked_kmer(kmer: u64, k: usize) -> u64 {
    let split_mask_extract = 3 << (k - 1);
    return kmer & !split_mask_extract;
}

#[inline]
pub fn split_kmer(kmer: u64, k: usize) -> (u64, u8) {
    let split_mask_extract = 3 << (k - 1);
    let mid_base = (kmer & split_mask_extract) >> (k - 1);
    let masked_kmer = kmer & !split_mask_extract;
    return (masked_kmer, mid_base as u8);
}

// Informational only: dump the k-mer count distribution and where the high-frequency
// cutoff landed on it, so the threshold can be sanity-checked against a dataset's actual
// depth. Never fatal -- a failure here just warns.
//
// Note the counts seen here are already post-filtering: k-mers with total count <= 2 or
// with single-strand support only were dropped during counting, so the low-count error
// tail is NOT represented. Bins start wherever the surviving data starts.
fn write_kmer_count_histogram(sorted_counts: &[u32], high_freq_thresh: u32, args: &Cli) {
    let n_total = sorted_counts.len();
    if n_total == 0 {
        return;
    }

    let out_dir = Path::new(&args.output_dir).join("misc");
    if let Err(e) = std::fs::create_dir_all(&out_dir) {
        log::warn!("Could not create {:?} for k-mer histogram: {}", out_dir, e);
        return;
    }
    let path = out_dir.join("kmer_count_histogram.tsv");

    // Exact bins for small counts, doubling bins above that.
    let max_count = *sorted_counts.last().unwrap() as u64;
    let mut edges: Vec<u64> = (1..100).collect();
    let mut edge = 100u64;
    while edge <= max_count {
        edges.push(edge);
        edge *= 2;
    }
    edges.push(max_count + 1);

    // The quantile before the .max(100) floor is applied, so it is visible whether the
    // floor or the quantile actually governed the cutoff.
    let quantile_thresh =
        sorted_counts[n_total - (n_total / args.high_freq_kmer_threshold) - 1];
    // Masking uses a strict >, so ties at the threshold survive.
    let n_masked = n_total - sorted_counts.partition_point(|&c| c <= high_freq_thresh);
    let total_occ: u64 = sorted_counts.iter().map(|&c| c as u64).sum();
    let masked_occ: u64 = sorted_counts[n_total - n_masked..]
        .iter()
        .map(|&c| c as u64)
        .sum();

    match File::create(&path) {
        Ok(f) => {
            let mut w = BufWriter::new(f);
            let _ = writeln!(w, "# n_distinct_kmers\t{}", n_total);
            let _ = writeln!(w, "# min_count\t{}", sorted_counts[0]);
            let _ = writeln!(w, "# max_count\t{}", max_count);
            let _ = writeln!(
                w,
                "# high_freq_kmer_threshold_arg\t{}",
                args.high_freq_kmer_threshold
            );
            let _ = writeln!(w, "# high_freq_thresh_quantile\t{}", quantile_thresh);
            let _ = writeln!(w, "# high_freq_thresh_applied\t{}", high_freq_thresh);
            let _ = writeln!(w, "# floor_100_was_binding\t{}", quantile_thresh < 100);
            let _ = writeln!(w, "# n_kmers_masked\t{}", n_masked);
            let _ = writeln!(
                w,
                "# frac_occurrences_masked\t{:.6}",
                masked_occ as f64 / total_occ as f64
            );
            let _ = writeln!(
                w,
                "count_lo\tcount_hi\tn_kmers\tfrac_kmers\tcum_frac_kmers\toccurrences"
            );

            // Single pass over the sorted counts; empty bins are skipped.
            let mut idx = 0usize;
            let mut cumulative = 0usize;
            for bin in 0..edges.len() - 1 {
                let lo = edges[bin];
                let hi = edges[bin + 1] - 1;
                let mut n = 0usize;
                let mut occurrences = 0u64;
                while idx < n_total && (sorted_counts[idx] as u64) <= hi {
                    n += 1;
                    occurrences += sorted_counts[idx] as u64;
                    idx += 1;
                }
                if n == 0 {
                    continue;
                }
                cumulative += n;
                let _ = writeln!(
                    w,
                    "{}\t{}\t{}\t{:.6}\t{:.6}\t{}",
                    lo,
                    hi,
                    n,
                    n as f64 / n_total as f64,
                    cumulative as f64 / n_total as f64,
                    occurrences
                );
            }
            log::debug!("Wrote k-mer count histogram to {:?}", path);
        }
        Err(e) => log::warn!("Could not write k-mer histogram to {:?}: {}", path, e),
    }
}

pub fn get_snpmers_inplace_sort(
    mut big_kmer_map: Vec<(Kmer64, [u32; 2])>,
    k: usize,
    args: &Cli,
) -> KmerGlobalInfo {
    assert!(!USE_SOLID_KMERS);

    log::debug!(
        "Number of k-mers passing thresholds: {}",
        big_kmer_map.len()
    );
    let mut kmer_counts = vec![];
    let high_freq_kmers = Arc::new(Mutex::new(HashSet::default()));
    let paths_to_files = args
        .input_files
        .iter()
        .map(|x| std::fs::canonicalize(Path::new(x).to_path_buf()).unwrap())
        .collect::<Vec<_>>();

    for pair in big_kmer_map.iter() {
        let counts = pair.1;
        kmer_counts.push(counts[0] + counts[1]);
    }

    kmer_counts.par_sort_unstable();
    if kmer_counts.len() == 0 {
        log::error!("No k-mers found. Exiting.");
        std::process::exit(1);
    }

    let high_freq_thresh = kmer_counts
        [kmer_counts.len() - (kmer_counts.len() / args.high_freq_kmer_threshold) - 1]
        .max(100);
    log::debug!("High frequency k-mer threshold: {}", high_freq_thresh);
    write_kmer_count_histogram(&kmer_counts, high_freq_thresh, args);
    drop(kmer_counts);

    log::info!("Finding snpmers...");
    //big_kmer_map.par_sort_unstable_by_key(|x| retrieve_masked_kmer(x.0, k));
    big_kmer_map.par_sort_unstable_by_key(|x| split_kmer(x.0, k));

    log::trace!("Finished parallel sort");

    let (mut tx, rx) = spmc::channel();
    let high_freq_kmers_arc = Arc::clone(&high_freq_kmers);
    thread::spawn(move || {
        let mut current_split_kmer = None;
        let mut kmer_pairs = vec![];
        for pair in big_kmer_map.into_iter() {
            let counts = pair.1;

            if counts[0] + counts[1] > high_freq_thresh {
                high_freq_kmers_arc
                    .lock()
                    .unwrap()
                    .insert(Kmer48::from_u64(pair.0));
            }

            if counts[0] == 0 || counts[1] == 0 {
                continue;
            }

            let kmer = pair.0;
            let split_kmer = retrieve_masked_kmer(kmer, k);

            if current_split_kmer != Some(split_kmer) {
                if kmer_pairs.len() > 1 {
                    tx.send(kmer_pairs).unwrap();
                }
                kmer_pairs = vec![];
                current_split_kmer = Some(split_kmer);
            }
            kmer_pairs.push(pair);
        }

        if kmer_pairs.len() > 1 {
            tx.send(kmer_pairs).unwrap();
        }
    });

    utils::log_memory_usage(false, "Memory usage during snpmer detection");

    if args.no_snpmers {
        log::info!("Skipping snpmer detection.");
        return KmerGlobalInfo {
            snpmer_info: vec![],
            solid_kmers: HashSet::default(),
            high_freq_kmers: Arc::try_unwrap(high_freq_kmers)
                .unwrap()
                .into_inner()
                .unwrap(),
            use_solid_kmers: USE_SOLID_KMERS,
            high_freq_thresh: high_freq_thresh as f64,
            read_files: paths_to_files,
        };
    }

    let potential_snps = Arc::new(Mutex::new(0));
    let snpmers = Arc::new(Mutex::new(vec![]));

    let mut handles = Vec::new();
    let k = args.kmer_size;

    for _ in 0..args.threads {
        let rx = rx.clone();
        let potential_snps = Arc::clone(&potential_snps);
        let snpmers = Arc::clone(&snpmers);
        handles.push(thread::spawn(move || {
            loop {
                match rx.recv() {
                    Ok(msg) => {
                        // SNPmers (not split) and counts
                        let mut pairsvec = msg;
                        pairsvec.sort_unstable_by(|a, b| (b.1[0] + b.1[1]).cmp(&(a.1[0] + a.1[1])));
                        let n = pairsvec[0].1[0] + pairsvec[0].1[1];
                        let succ = pairsvec[1].1[0] + pairsvec[1].1[1];
                        let right_p_val_thresh1 =
                            twin_graph::binomial_test(n as u64, succ as u64, 0.025);
                        let right_p_val_thresh2 =
                            twin_graph::binomial_test(n as u64, succ as u64, 0.050);
                        let cond1 = right_p_val_thresh1 > 0.05;
                        let cond2 = right_p_val_thresh2 > 0.05 && k < 5;

                        if cond1 || cond2 {
                            if log::log_enabled!(log::Level::Trace) {
                                let snpmer1 = pairsvec[0].0;
                                let snpmer2 = pairsvec[1].0;
                                log::trace!(
                                    "NOT SNPMER BINOMIAL {} {} c:{:?} c:{:?}",
                                    decode_kmer64(snpmer1, k as u8),
                                    decode_kmer64(snpmer2, k as u8),
                                    pairsvec[0].1,
                                    pairsvec[1].1
                                );
                            }
                            continue;
                        }

                        let a = pairsvec[0].1[0];
                        let b = pairsvec[1].1[0];
                        let c = pairsvec[0].1[1];
                        let d = pairsvec[1].1[1];
                        let contingency_table = [a.max(c), b.max(d), c.min(a), d.min(b)];
                        let p_value = fishers_exact(&contingency_table).unwrap().two_tail_pvalue;
                        let odds;
                        if contingency_table[0] == 0
                            || contingency_table[1] == 0
                            || contingency_table[2] == 0
                            || contingency_table[3] == 0
                        {
                            odds = 0.0;
                        } else {
                            odds = (contingency_table[0] as f64 * contingency_table[3] as f64)
                                / (contingency_table[1] as f64 * contingency_table[2] as f64);
                        }
                        if odds == 0. {
                            continue;
                        }

                        //Is snpmer
                        if p_value > 0.005 || (odds < 1.5 && odds > 1. / 1.5) {
                            let splitmer = retrieve_masked_kmer(pairsvec[0].0, k);
                            let mid_bases = pairsvec
                                .iter()
                                .map(|x| split_kmer(x.0, k).1)
                                .collect::<Vec<_>>();
                            let snpmer = SnpmerInfo {
                                split_kmer: splitmer,
                                mid_bases: [mid_bases[0] as u8, mid_bases[1] as u8],
                                counts: [
                                    pairsvec[0].1[0] + pairsvec[0].1[1],
                                    pairsvec[1].1[0] + pairsvec[1].1[1],
                                ],
                            };
                            snpmers.lock().unwrap().push(snpmer);

                            let snpmer1 = splitmer as u64 | ((mid_bases[0] as u64) << (k - 1));
                            let snpmer2 = splitmer as u64 | ((mid_bases[1] as u64) << (k - 1));
                            log::trace!(
                                "{} c:{:?} {} c:{:?}, p:{}, odds:{}",
                                decode_kmer64(snpmer1, k as u8),
                                pairsvec[0].1,
                                decode_kmer64(snpmer2, k as u8),
                                pairsvec[1].1,
                                p_value,
                                odds
                            );
                            *potential_snps.lock().unwrap() += 1;
                        } else {
                            log::trace!(
                                "NOT SNPMER c:{:?} c:{:?}, p:{}, odds:{}",
                                pairsvec[0].1,
                                pairsvec[1].1,
                                p_value,
                                odds
                            );
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

    let mut snpmers = Arc::try_unwrap(snpmers).unwrap().into_inner().unwrap();
    snpmers.sort();
    log::debug!("Number of snpmers: {}. ", potential_snps.lock().unwrap());
    return KmerGlobalInfo {
        snpmer_info: snpmers,
        solid_kmers: HashSet::default(),
        high_freq_kmers: Arc::try_unwrap(high_freq_kmers)
            .unwrap()
            .into_inner()
            .unwrap(),
        use_solid_kmers: USE_SOLID_KMERS,
        high_freq_thresh: high_freq_thresh as f64,
        read_files: paths_to_files,
    };
}

pub fn get_snpmers(big_kmer_map: Vec<(Kmer64, [u32; 2])>, k: usize, args: &Cli) -> KmerGlobalInfo {
    log::debug!(
        "Number of k-mers passing thresholds: {}",
        big_kmer_map.len()
    );
    let mut new_map_counts_bases: FxHashMap<Kmer64, CountsAndBases> = FxHashMap::default();
    let mut kmer_counts = vec![];
    let mut solid_kmers = HashSet::default();
    let mut high_freq_kmers = HashSet::default();
    let paths_to_files = args
        .input_files
        .iter()
        .map(|x| std::fs::canonicalize(Path::new(x).to_path_buf()).unwrap())
        .collect::<Vec<_>>();

    for pair in big_kmer_map.iter() {
        let counts = pair.1;
        kmer_counts.push(counts[0] + counts[1]);
    }

    kmer_counts.par_sort_unstable();
    if kmer_counts.len() == 0 {
        log::error!("No k-mers found. Exiting.");
        std::process::exit(1);
    }
    let high_freq_thresh =
        kmer_counts[kmer_counts.len() - (kmer_counts.len() / 100000) - 1].max(100);
    log::debug!("High frequency k-mer threshold: {}", high_freq_thresh);
    drop(kmer_counts);

    log::info!("Finding snpmers...");
    //Should be able to parallelize this, TODO
    for pair in big_kmer_map.into_iter() {
        let kmer = pair.0;
        let (split_kmer, mid_base) = split_kmer(kmer, k);
        let counts = pair.1;
        if counts[0] > 0 && counts[1] > 0 {
            let count = counts[0] + counts[1];
            if count < high_freq_thresh {
                solid_kmers.insert(Kmer48::from_u64(kmer));
                let v = new_map_counts_bases
                    .entry(split_kmer)
                    .or_insert(CountsAndBases {
                        counts: SmallVec::new(),
                        bases: SmallVec::new(),
                    });
                v.counts.push(counts);
                v.bases.push(mid_base);
            } else {
                high_freq_kmers.insert(Kmer48::from_u64(kmer));
            }
        }
    }

    utils::log_memory_usage(false, "Memory usage during snpmer detection");

    if args.no_snpmers {
        log::info!("Skipping snpmer detection.");
        return KmerGlobalInfo {
            snpmer_info: vec![],
            solid_kmers: solid_kmers,
            high_freq_kmers: high_freq_kmers,
            use_solid_kmers: USE_SOLID_KMERS,
            high_freq_thresh: high_freq_thresh as f64,
            read_files: paths_to_files,
        };
    }

    let potential_snps = Mutex::new(0);
    let snpmers = Mutex::new(vec![]);
    new_map_counts_bases
        .into_par_iter()
        .for_each(|(split_kmer, c_and_b)| {
            let mut counts = c_and_b.counts;
            let bases = c_and_b.bases;
            if counts.len() > 1 {
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
                    if log::log_enabled!(log::Level::Trace) {
                        let mid_bases = bases;
                        let snpmer1 = split_kmer as u64 | ((mid_bases[0] as u64) << (k - 1));
                        let snpmer2 = split_kmer as u64 | ((mid_bases[1] as u64) << (k - 1));
                        log::trace!(
                            "NOT SNPMER BINOMIAL {} {} c:{:?} c:{:?}",
                            decode_kmer64(snpmer1, k as u8),
                            decode_kmer64(snpmer2, k as u8),
                            counts[0],
                            counts[1]
                        );
                    }
                    return;
                }

                //Add pseudocount... especially when all reads are
                //already in forward strand (happens when debugging or manipulation via samtools)
                let a = counts[0][0];
                let b = counts[1][0];
                let c = counts[0][1];
                let d = counts[1][1];
                let contingency_table = [a.max(c), b.max(d), c.min(a), d.min(b)];
                let p_value = fishers_exact(&contingency_table).unwrap().two_tail_pvalue;
                let odds;
                if contingency_table[0] == 0
                    || contingency_table[1] == 0
                    || contingency_table[2] == 0
                    || contingency_table[3] == 0
                {
                    odds = 0.0;
                } else {
                    odds = (contingency_table[0] as f64 * contingency_table[3] as f64)
                        / (contingency_table[1] as f64 * contingency_table[2] as f64);
                }

                if odds == 0. {
                    return;
                }

                //Is snpmer
                if p_value > 0.005 || (odds < 1.5 && odds > 1. / 1.5) {
                    let mid_bases = bases;
                    let snpmer = SnpmerInfo {
                        split_kmer: split_kmer,
                        mid_bases: [mid_bases[0], mid_bases[1]],
                        counts: [counts[0][0] + counts[0][1], counts[1][0] + counts[1][1]],
                    };
                    snpmers.lock().unwrap().push(snpmer);

                    let snpmer1 = split_kmer as u64 | ((mid_bases[0] as u64) << (k - 1));
                    let snpmer2 = split_kmer as u64 | ((mid_bases[1] as u64) << (k - 1));
                    log::trace!(
                        "{} c:{:?} {} c:{:?}, p:{}, odds:{}",
                        decode_kmer64(snpmer1, k as u8),
                        counts[0],
                        decode_kmer64(snpmer2, k as u8),
                        counts[1],
                        p_value,
                        odds
                    );
                    *potential_snps.lock().unwrap() += 1;
                } else {
                    log::trace!(
                        "NOT SNPMER c:{:?} c:{:?}, p:{}, odds:{}",
                        counts[0],
                        counts[1],
                        p_value,
                        odds
                    );
                }
            }
        });

    let mut snpmers = snpmers.into_inner().unwrap();
    snpmers.sort();
    solid_kmers.shrink_to_fit();
    log::debug!(
        "Number of snpmers: {}. ",
        potential_snps.into_inner().unwrap()
    );
    log::debug!("Number of solid k-mers: {}.", solid_kmers.len());
    return KmerGlobalInfo {
        snpmer_info: snpmers,
        solid_kmers: solid_kmers,
        high_freq_kmers: high_freq_kmers,
        use_solid_kmers: USE_SOLID_KMERS,
        high_freq_thresh: high_freq_thresh as f64,
        read_files: paths_to_files,
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
