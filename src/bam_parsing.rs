use rust_htslib::bam::{self, Read, IndexedReader};
use needletail::{parse_fastx_file, Sequence, FastxReader};
use rust_htslib::bam::pileup::Pileup;
use rust_htslib::bcf::{Format, Writer};
use rust_htslib::bcf::header::Header;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib;
use crate::seeding::fmh_seeds_positions;
use std::str;
use crate::types::*;
use fxhash::FxHashMap;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::*;
use statrs::distribution::DiscreteCDF;
use statrs::distribution::{Binomial, Discrete};

pub fn parse_bam(bam_file: &str) {
    let block = 3000;

    let mut bam = Reader::from_path(bam_file).unwrap();
    let header = bam.header().to_owned();
    let target_names = header.target_names();
    let target_ids = target_names
        .iter()
        .map(|x| header.tid(x).unwrap())
        .collect::<Vec<u32>>();
    let target_lens = target_ids
        .iter()
        .map(|x| header.target_len(*x).unwrap())
        .collect::<Vec<u64>>();
    let mut record = Record::new();

    let mut minimizer_index: Vec<FxHashMap<u64, (u32, u64)>> =
        vec![FxHashMap::default(); target_lens[0] as usize / block];
    let mut coverage_index: Vec<f64> = vec![0.0; target_lens[0] as usize / block];
    let mut prefrags = vec![];
    let c = 25;
    let k = 16;

    while let Some(r) = bam.read(&mut record) {
        let prefrag = process_bam_record(
            &record,
            block,
            &mut minimizer_index,
            &mut coverage_index,
            c,
            k,
        );

        if let Some(prefrag) = prefrag {
            prefrags.push(prefrag);
        }
    }

    let mut num_kmers_g1 = 0;
    let mut num_kmers_1 = 0;
    for i in 0..minimizer_index.len() {
        for (_, v) in minimizer_index[i].iter() {
            if v.0 == 1 {
                num_kmers_1 += 1;
            } else if v.0 > 1 {
                num_kmers_g1 += v.0;
            }
        }
        minimizer_index[i].retain(|_, &mut v| v.0 > 2);
    }
    let err_rate = num_kmers_g1 as f64 / (num_kmers_g1 + num_kmers_1) as f64;
    let mut all_varmers = vec![];
    for i in 0..minimizer_index.len() {
        let vars = get_varmers(
            coverage_index[i],
            &minimizer_index[i],
            err_rate,
            k as u8,
            0.005,
        );
        if vars.is_empty() {
            break;
        }
        println!("{:?},{}", vars, err_rate);
        all_varmers.extend(vars);
    }

    all_varmers.sort_by(|a, b| a.pos.cmp(&b.pos));

    let mut var_frags = vec![];

    for pre in prefrags.iter() {
        let mut var_frag = VarmerFragment::default();

        let first_ind = all_varmers
            .binary_search_by(|x| x.pos.cmp(&pre.kmers_with_refpos[0].1))
            .unwrap_or_else(|x| x);
        let last_ind = all_varmers
            .binary_search_by(|x| {
                x.pos
                    .cmp(&pre.kmers_with_refpos[pre.kmers_with_refpos.len() - 1].1)
            })
            .unwrap_or_else(|x| x);

        var_frag.lower = first_ind;
        var_frag.upper = last_ind;
        let mut varmer_indices = vec![];

        let mut track_i = first_ind;
        for (kmer, refpos) in pre.kmers_with_refpos.iter() {
            if track_i >= last_ind {
                break;
            }
            for i in track_i..last_ind {
                let vark = all_varmers[i].kmer;
                if vark == *kmer{
                    println!("{:?},{}", all_varmers[i], refpos);
                    varmer_indices.push(i);
                    track_i = i+1;
                    break;
                }
                else if refpos < &all_varmers[i].pos {
                    track_i = i;
                    break;
                }
            }
        }

        var_frag.varmers = varmer_indices.into_iter().collect();
        var_frag.lower_base = pre.lower_base;
        var_frag.upper_base = pre.upper_base;
        var_frag.id = pre.id.clone();
        var_frags.push(var_frag);
    }

    var_frags.sort_by(|a, b| a.lower.cmp(&b.lower));

    println!("{:?}", all_varmers.len());
    println!("{:?}", var_frags);

    print_graph(var_frags);
}

fn process_bam_record(
    record: &Record,
    block: usize,
    minimizer_index: &mut Vec<FxHashMap<u64, (u32, u64)>>,
    coverage_index: &mut Vec<f64>,
    c: usize,
    k: usize,
) -> Option<PreFragment> {
    let seq = record.seq();
    //check if primary alignment
    if record.is_secondary() || record.is_unmapped() || record.is_duplicate() {
        return None;
    }
    if record.is_supplementary() && record.mapq() < 30 {
        return None;
    }
    if record.mapq() < 15 {
        return None;
    }

    let mut pre_frag = PreFragment::default();
    //convert [u8] to String
    pre_frag.id = String::from_utf8(record.qname().to_vec()).unwrap().to_string();
    pre_frag.lower_base = record.reference_start() as usize;
    pre_frag.upper_base = record.reference_end() as usize;
    //Get start and end alignments
    let start = (record.reference_start() as f64) / block as f64;
    let end = (record.reference_end() as f64) / block as f64;

    // round start down and end up but keep track of how much rounding error there is
    let s = start.floor() as usize;
    let e = end.ceil() as usize;

    let end_overest = e as f64 - end;
    let start_overest = start - s as f64;

    //let start = (record.reference_start() as usize + block/4) / block;
    //let end = (record.reference_end() as usize - block/4) / block;
    for i in s..e {
        if i == s {
            coverage_index[i] += 1. - start_overest;
        } else if i == e - 1 {
            coverage_index[i] += 1. - end_overest;
        } else {
            coverage_index[i as usize] += 1.0;
        }
    }

    // Get minimizer seeds with positions
    let mut minimizers = vec![];
    let mut positions = vec![];
    //minimizer_seeds_positions(&seq.as_bytes(), &mut minimizers, &mut positions, c*2, k);
    fmh_seeds_positions(&seq.as_bytes(), &mut minimizers, &mut positions, c, k);

    if minimizers.len() == 0 {
        return None;
    }

    let aligned_pos = record.aligned_pairs();
    let mut curr_ind = 0;
    let mut curr_pos = positions[0];
    for x in aligned_pos {
        if x[0] < curr_pos as i64 {
            continue;
        } else if x[0] >= curr_pos as i64 {
            let hash = minimizers[curr_ind] as u64;
            let index = x[1] as usize / block;

            let c = minimizer_index[index]
                .entry(hash)
                .or_insert((0, x[1] as u64));
            pre_frag.kmers_with_refpos.push((hash, x[1] as u64));
            c.0 += 1;
        }

        if curr_ind == minimizers.len() - 1 {
            break;
        }
        curr_ind += 1;
        curr_pos = positions[curr_ind];
    }

    return Some(pre_frag);
}

fn get_varmers(
    cov: f64,
    kmer_counts: &FxHashMap<u64, (u32, u64)>,
    err_rate: f64,
    k: u8,
    p: f64,
) -> Vec<Varmer> {
    let binom = Binomial::new(err_rate, cov as u64).unwrap();
    let mut cutoff = 0;
    let mut sum_p = 0.;

    while sum_p < p {
        sum_p += binom.pmf(cutoff);
        cutoff += 1;
    }

    let mut varmers = vec![];
    for (kmer, v) in kmer_counts.iter() {
        if v.0 < cutoff as u32 {
            if binom.cdf(v.0 as u64) < 0.05 {
                let varmer = Varmer {
                    kmer: *kmer,
                    count: v.0,
                    pos: v.1,
                };
                varmers.push(varmer);
            }
        }
    }
    varmers
}

fn print_graph(var_frags: Vec<VarmerFragment>)
{
//print a CSV file representing a graph with varfrags. Edges represent # of shared varmers /
    //overlap between varfrags
    for i in 0..var_frags.len() {
        for j in i + 1..var_frags.len() {
            let mut overlap = (0, 0);
            let x1 = var_frags[i].lower;
            let x2 = var_frags[i].upper;
            let y1 = var_frags[j].lower;
            let y2 = var_frags[j].upper;
            if x1 < y1 {
                if x2 < y1 {
                    continue;
                } else if x2 < y2 {
                    overlap = (y1, x2);
                } else {
                    overlap = (y1, y2);
                }
            } else {
                if y2 < x1 {
                    continue;
                } else if y2 < x2 {
                    overlap = (x1, y2);
                } else {
                    overlap = (x1, x2);
                }
            }
            if overlap.0 == overlap.1{
                continue
            }
            //count the number of varmers in overlap
            let varmer_in_overlap_i = var_frags[i].varmers.iter().filter(|x| {
                **x >= overlap.0 && **x <= overlap.1
            }).count();

            let varmer_in_overlap_j = var_frags[j].varmers.iter().filter(|x| {
                **x >= overlap.0 && **x <= overlap.1
            }).count();
            
            let avg_overlap = (varmer_in_overlap_i + varmer_in_overlap_j) as f64 / 2.0;

            let shared = var_frags[i].varmers.intersection(&var_frags[j].varmers).count();
            let edge = shared as f64;
            println!("{},{},{},{},{},{},{},{},{},{}", i, j, edge, avg_overlap, x1, x2, y1,y2, var_frags[i].id, var_frags[j].id);
        }
    }

}

pub fn get_pileup(bam_file: &str, fasta_file: &str) -> FxHashMap<Vec<u8>, Vec<i64>> {
    let mut bam = IndexedReader::from_path(&bam_file).expect("Error opening BAM file");
    let mut header = bam.header().clone();
    let mut return_map = FxHashMap::default();
    let mut header = Header::new();
    let mut reader = parse_fastx_file(&fasta_file).expect("valid path/file");

    while let Some(r) = reader.next() {
        let record = r.expect("Error reading record");
        let id = record.id().to_owned();
        let id_seq = str::from_utf8(&id).unwrap().split_whitespace().next().unwrap().as_bytes().to_vec();
        //let header_contig_line = r#"##contig=<ID=1,length=10>"#;
        let header_contig_line = format!("##contig=<ID={},length={}>", str::from_utf8(&id_seq).unwrap(), record.seq().len());
        header.push_record(header_contig_line.as_bytes());
        bam.fetch(&id_seq);
        let num_recs = bam.records().count();
        if num_recs == 0 {
            continue;
        }

        let mut pileup_vec_for_contig = vec![];
        for (i,base) in record.seq().iter().enumerate() {
            let new_p = BasePileup{
                ref_pos : i as u64,
                ref_base: base_to_index(*base),
                base_freqs: [0; 4],
            };
            pileup_vec_for_contig.push(new_p);
        }

        bam.fetch(&id_seq);
        let mut printed = false;
        for r in bam.records(){
            let record = r.unwrap();
            if record.is_secondary() || record.is_supplementary() || record.mapq() < 10{
                continue;
            }

            if !printed{
                log::debug!("{:?}", str::from_utf8(&id_seq));
                printed = true;
            }
            for aligned_pairs in record.aligned_pairs(){
                let ref_pos = aligned_pairs[1];
                let pos = aligned_pairs[0];
                let base = BYTE_TO_SEQ[record.seq()[pos as usize] as usize];
                pileup_vec_for_contig[ref_pos as usize].base_freqs[base as usize] += 1;
            }
        }
        let snv_pos = get_snvs(&pileup_vec_for_contig);
        return_map.insert(id_seq, snv_pos);
    }
    let header_gt_line = r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#;
    header.push_record(header_gt_line.as_bytes());
    header.push_sample("test_sample".as_bytes());
    let mut vcf = Writer::from_stdout(&header, true, Format::Vcf).unwrap();
    for (contig, snvs) in return_map.iter(){
        for snv in snvs.iter(){
            let mut record = vcf.empty_record();
            let rid = vcf.header().name2rid(contig).unwrap();
            record.set_rid(Some(rid));
            record.set_pos(snv.0 as i64);
            record.set_alleles(&[&vec![snv.1], &vec![snv.2]]);
            record.push_genotypes(&[GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)]);
            vcf.write(&record).unwrap();
        }
    }
    return FxHashMap::default();
}

pub fn get_snvs(pileups: &Vec<BasePileup> ) -> Vec<(u64, u8, u8)>{

    let mut ret_vec = vec![];
    for p in pileups.iter(){
        let sum = p.base_freqs[0] + p.base_freqs[1] + p.base_freqs[2] + p.base_freqs[3];
        for i in 0..4{
           if i == p.ref_base as usize{
               continue;
           }
           let freq =  p.base_freqs[i] as f32 / sum as f32;
           if freq > 0.10 && freq < 0.80 && p.base_freqs[i] > 5 {
                ret_vec.push((p.ref_pos, index_to_base(p.ref_base), index_to_base(i as u8)));
                log::trace!("{:?}, {:?}", p.ref_pos, p.base_freqs);
                break;
            }
        }
    }
    return ret_vec;

}

fn index_to_base(i: u8) -> u8{
    match i{
        0 => 'A' as u8,
        1 => 'C' as u8,
        2 => 'G' as u8,
        3 => 'T' as u8,
        _ => 0,
    }
}

fn base_to_index(base: u8) -> u8{
    match base{
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}
