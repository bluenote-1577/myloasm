use std::io::BufWriter;
use fxhash::FxHashMap;
use std::path::Path;
use crate::cli::Cli;
use crate::constants::{CIRC_LAX_STRING, CIRC_STRICT_STRING};
use std::fs::File;
use std::io::Write;

use skani::params;
pub fn dereplicate_with_skani(polished_fasta: &str, args: &Cli){

    let initial_fasta_path = Path::new(args.output_dir.as_str()).join(polished_fasta);
    assert!(initial_fasta_path.exists());

    let command_params = params::CommandParams {
        ref_files: vec![initial_fasta_path.to_str().unwrap().to_string()],
        individual_contig_q: true,
        individual_contig_r: true,
        min_aligned_frac: 0.15,
        learned_ani: true,
        rescue_small: false,
        mode: params::Mode::Triangle,
        screen: true,
        screen_val: 0.9,
        ..Default::default()
    };

    let sketch_params = params::SketchParams{
        c: 15,
        k: 15,
        marker_c: 50,
        ..Default::default()
    };

    // Remove (smaller, lower_id) if > 90% aligned and > 99% ANI
    let ani_results = skani::triangle::triangle_return(command_params, sketch_params);
    let mut alternate_indices = FxHashMap::default();
    let mut duplicate_indices = FxHashMap::default();
    let mut contig_name_to_id_map = FxHashMap::default();

    for (ref_id, query_results) in ani_results.iter(){
        for (query_id, ani_result) in query_results.iter(){

            if ani_result.quant_50_contig_len_r > 500_000. && ani_result.quant_50_contig_len_q > 500_000.{
                continue;
            }

            // STEP 1: REMOVE REPETITIVE, misassembled THINGS
            //If something with huge multiplicity consumes a small k-mer multiplicity object,
            // it's likely a duplicated repeat. Remove it. 
            if (ani_result.align_fraction_ref > 0.9 || ani_result.align_fraction_query > 0.9) && ani_result.ani > 0.98{
                let length_id_query = (ani_result.quant_50_contig_len_q, query_id, &ani_result.query_contig);
                let length_id_ref = (ani_result.quant_50_contig_len_r, ref_id, &ani_result.ref_contig);
                let smaller = if length_id_query < length_id_ref {length_id_query} else {length_id_ref};
                let larger = if smaller == length_id_query {length_id_ref} else {length_id_query};
                let kmer_mult_larger_str = larger.2.split("mult=").nth(1).unwrap();
                dbg!(&kmer_mult_larger_str);
                let kmer_mult_larger = kmer_mult_larger_str.parse::<f64>().unwrap_or(1.00);

                if kmer_mult_larger > 1.25 && larger.0 < 1_000_000. {
                    duplicate_indices.insert(larger.1, (smaller.2, ani_result));
                    contig_name_to_id_map.insert(ani_result.query_contig.clone(), query_id);
                    contig_name_to_id_map.insert(ani_result.ref_contig.clone(), ref_id);
                    continue;
                }
            }

            // STEP 2: REMOVE DUPLICATED things
            if (ani_result.align_fraction_ref > 0.99 || ani_result.align_fraction_query > 0.99) && ani_result.ani > 0.999 {
                let length_id_query = (ani_result.quant_50_contig_len_q, query_id, &ani_result.query_contig);
                let length_id_ref = (ani_result.quant_50_contig_len_r, ref_id, &ani_result.ref_contig);
                let mut smaller = if length_id_query < length_id_ref {length_id_query} else {length_id_ref};

                //If they have similar size, keep the circular one
                if length_id_query.0 / length_id_ref.0 > 0.66 && length_id_query.0 / length_id_ref.0 < 1./0.66{
                    if length_id_query.2.contains(CIRC_STRICT_STRING) && !length_id_ref.2.contains(CIRC_STRICT_STRING){
                        smaller = length_id_ref;
                    }
                    else if length_id_ref.2.contains(CIRC_STRICT_STRING) && !length_id_query.2.contains(CIRC_STRICT_STRING){
                        smaller = length_id_query;
                    }
                }
                let larger = if smaller == length_id_query {length_id_ref} else {length_id_query};

                if smaller.0 < 150_000.{
                    duplicate_indices.insert(smaller.1, (larger.2, ani_result));
                    contig_name_to_id_map.insert(ani_result.query_contig.clone(), query_id);
                    contig_name_to_id_map.insert(ani_result.ref_contig.clone(), ref_id);
                }
            }

            // STEP 3: Remove ALTERNATE things
            else if (ani_result.align_fraction_ref > 0.90 || ani_result.align_fraction_query > 0.90) && ani_result.ani > 0.99 {

                let length_id_query = (ani_result.quant_50_contig_len_q, query_id, &ani_result.query_contig);
                let length_id_ref = (ani_result.quant_50_contig_len_r, ref_id, &ani_result.ref_contig);
                let smaller = if length_id_query < length_id_ref {length_id_query} else {length_id_ref};
                let larger = if smaller == length_id_query {length_id_ref} else {length_id_query};

                // Don't allow circular contigs to be considered alternate if they are not duplicates
                // If smaller is the query, remove the reference
                if smaller.0 < 500_000. && !smaller.2.contains(CIRC_STRICT_STRING){
                    alternate_indices.insert(smaller.1, (larger.2, ani_result));
                    contig_name_to_id_map.insert(ani_result.query_contig.clone(), query_id);
                    contig_name_to_id_map.insert(ani_result.ref_contig.clone(), ref_id);
                }
                
                // Really similar small things, doesn't matter if circular, dereplicate. 
                // Small things could arise by sequencing artifacts anyways. 
                else if ani_result.align_fraction_query > 0.90 && ani_result.align_fraction_ref > 0.90  && smaller.0 < 100_000. {
                    alternate_indices.insert(smaller.1, (larger.2, ani_result));
                    contig_name_to_id_map.insert(ani_result.query_contig.clone(), query_id);
                    contig_name_to_id_map.insert(ani_result.ref_contig.clone(), ref_id);
                }
            }
        }
    }

    let mut fasta_records = needletail::parse_fastx_file(&initial_fasta_path).unwrap();
    let alternate_assembly_dir = Path::new(args.output_dir.as_str()).join("alternate_assemblies");
    std::fs::create_dir_all(&alternate_assembly_dir).unwrap();
    let primary_path = Path::new(args.output_dir.as_str()).join("assembly_primary.fa");
    let alternate_path = alternate_assembly_dir.join("assembly_alternate.fa");
    let dup_path = alternate_assembly_dir.join("duplicated_contigs.fa");
    let mut fasta_writer_primary = BufWriter::new(File::create(primary_path).unwrap());
    let mut fasta_writer_alternate = BufWriter::new(File::create(alternate_path).unwrap());
    let mut fasta_writer_duplicated = BufWriter::new(File::create(dup_path).unwrap());
    let mut total_bases_primary = 0;
    let mut number_primary = 0;
    let mut number_alternate = 0;
    let mut number_duplicated = 0;
    let mut largest_contig_size = 0;
    let mut contig_sizes = vec![];
    let mut num_circular_1m = 0;

    while let Some(record) = fasta_records.next() {
        let record = record.unwrap();
        let header = String::from_utf8(record.id().to_vec()).unwrap();
        let seq = record.seq().iter().map(|x| *x as char).collect::<String>();
        let id_opt = contig_name_to_id_map.get(&header);
        if (header.contains(CIRC_STRICT_STRING) || header.contains(CIRC_LAX_STRING)) && seq.len() > 1_000_000{
            num_circular_1m += 1;
        }

        if let Some(id) = id_opt{
            if let Some(val) = duplicate_indices.get(id){
                let af = if val.1.align_fraction_query > val.1.align_fraction_ref {val.1.align_fraction_query} else {val.1.align_fraction_ref};
                let ani = val.1.ani * 100.;
                number_duplicated += 1;
                write!(&mut fasta_writer_duplicated, ">{} secondary_to:|{}| aligned_frac:{:.2} ani:{:.2}\n{}\n", header, val.0, af*100., ani, seq).unwrap();
            }
            else if let Some(val) = alternate_indices.get(id){
                let af = if val.1.align_fraction_query > val.1.align_fraction_ref {val.1.align_fraction_query} else {val.1.align_fraction_ref};
                let ani = val.1.ani * 100.;
                number_alternate += 1;
                write!(&mut fasta_writer_alternate, ">{} secondary_to:|{}| aligned_frac:{:.2} ani:{:.2}\n{}\n", header, val.0, af*100., ani, seq).unwrap();
            }
            else{
                total_bases_primary += seq.len();
                number_primary += 1;
                write!(&mut fasta_writer_primary, ">{}\n{}\n", header, seq).unwrap();
            }
        }
        else{
            total_bases_primary += seq.len();
            number_primary += 1;
            write!(&mut fasta_writer_primary, ">{}\n{}\n", header, seq).unwrap();
        }

        contig_sizes.push(seq.len());
        if seq.len() > largest_contig_size{
            largest_contig_size = seq.len();
        }
    }

    contig_sizes.sort_by(|a, b| b.cmp(a));
    let mut iter_sum_size = 0;
    let mut num_geq_100 = 0;
    let mut num_geq_1000 = 0;
    let mut n50 = None;

    for (_, &size) in contig_sizes.iter().enumerate(){
        if size > 100_000{
            num_geq_100 += 1;
        }
        if size > 1_000_000{
            num_geq_1000 += 1;
        }
        iter_sum_size += size;
        if iter_sum_size > total_bases_primary / 2 && n50.is_none() {
            n50 = Some(size);
        }
    }

    log::info!("-------------- FINAL primary assembly statistics --------------");
    log::info!("N50: {}", n50.unwrap_or(0));
    log::info!("Largest contig has size: {}", largest_contig_size);
    log::info!("Number of primary contigs: {}", number_primary);
    log::info!("Number of alt/dup contigs: {} and {}", number_alternate, number_duplicated);
    log::info!(
        "Number of possibly circular contigs >= 1M: {}",
        num_circular_1m
    );
    log::info!(
        "Number of contigs >= 1M: {}",
        num_geq_1000
    );
    log::info!(
        "Number of contigs >= 100k: {}",
        num_geq_100
    );
    log::info!("Total bases within assembly is {}", total_bases_primary);
    log::info!("-------------------------------------------------");

    // Remove the polished fasta
    let _ = std::fs::remove_file(initial_fasta_path);

}

