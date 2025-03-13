use std::io::BufWriter;
use fxhash::FxHashMap;
use std::path::Path;
use crate::cli::Cli;
use std::fs::File;
use std::io::Write;

use skani::params;
pub fn dereplicate_with_skani(polished_fasta: &str, args: &Cli){

    let path = Path::new(args.output_dir.as_str()).join(polished_fasta);
    assert!(path.exists());

    let command_params = params::CommandParams {
        ref_files: vec![path.to_str().unwrap().to_string()],
        individual_contig_q: true,
        individual_contig_r: true,
        min_aligned_frac: 0.30,
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

            //Too similar - duplicate, not alternate
            if (ani_result.align_fraction_ref > 0.99 || ani_result.align_fraction_query > 0.99) && ani_result.ani > 0.999 {
                let length_id_query = (ani_result.quant_50_contig_len_q, query_id, &ani_result.query_contig);
                let length_id_ref = (ani_result.quant_50_contig_len_r, ref_id, &ani_result.ref_contig);
                let mut smaller = if length_id_query < length_id_ref {length_id_query} else {length_id_ref};

                //If they have similar size, keep the circular one
                if length_id_query.0 / length_id_ref.0 > 0.66 && length_id_query.0 / length_id_ref.0 < 1./0.66{
                    if length_id_query.2.contains("circular") && !length_id_ref.2.contains("circular"){
                        smaller = length_id_ref;
                    }
                    else if length_id_ref.2.contains("circular") && !length_id_query.2.contains("circular"){
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

            else if (ani_result.align_fraction_ref > 0.90 || ani_result.align_fraction_query > 0.90) && ani_result.ani > 0.99 {

                // Keep non-duplicated circular contigs as primary
                if ani_result.query_contig.contains("circular") || ani_result.ref_contig.contains("circular"){
                    continue;
                }

                let length_id_query = (ani_result.quant_50_contig_len_q, query_id, &ani_result.query_contig);
                let length_id_ref = (ani_result.quant_50_contig_len_r, ref_id, &ani_result.ref_contig);
                let smaller = if length_id_query < length_id_ref {length_id_query} else {length_id_ref};
                let larger = if smaller == length_id_query {length_id_ref} else {length_id_query};

                // If smaller is the query, remove the reference
                if smaller.0 < 500_000.{
                    alternate_indices.insert(smaller.1, (larger.2, ani_result));
                    contig_name_to_id_map.insert(ani_result.query_contig.clone(), query_id);
                    contig_name_to_id_map.insert(ani_result.ref_contig.clone(), ref_id);
                }
            }
        }
    }

    let mut fasta_records = needletail::parse_fastx_file(path).unwrap();
    let primary_path = Path::new(args.output_dir.as_str()).join("assembly_primary.fa");
    let alternate_path = Path::new(args.output_dir.as_str()).join("assembly_alternate.fa");
    let dup_path = Path::new(args.output_dir.as_str()).join("duplicated_contigs.fa");
    let mut fasta_writer_primary = BufWriter::new(File::create(primary_path).unwrap());
    let mut fasta_writer_alternate = BufWriter::new(File::create(alternate_path).unwrap());
    let mut fasta_writer_duplicated = BufWriter::new(File::create(dup_path).unwrap());

    while let Some(record) = fasta_records.next() {
        let record = record.unwrap();
        let header = String::from_utf8(record.id().to_vec()).unwrap();
        let seq = record.seq().iter().map(|x| *x as char).collect::<String>();
        let id_opt = contig_name_to_id_map.get(&header);

        if let Some(id) = id_opt{
            if let Some(val) = duplicate_indices.get(id){
                let af = if val.1.align_fraction_query > val.1.align_fraction_ref {val.1.align_fraction_query} else {val.1.align_fraction_ref};
                let ani = val.1.ani * 100.;
                write!(&mut fasta_writer_duplicated, ">{} secondary_to:({}) aligned_frac:{:.2} ani:{:.2}\n{}\n", header, val.0, af*100., ani, seq).unwrap();
            }
            else if let Some(val) = alternate_indices.get(id){
                let af = if val.1.align_fraction_query > val.1.align_fraction_ref {val.1.align_fraction_query} else {val.1.align_fraction_ref};
                let ani = val.1.ani * 100.;
                write!(&mut fasta_writer_alternate, ">{} secondary_to:[{}] aligned_frac:{:.2} ani:{:.2}\n{}\n", header, val.0, af*100., ani, seq).unwrap();
            }
            else{
                write!(&mut fasta_writer_primary, ">{}\n{}\n", header, seq).unwrap();
            }
        }
        else{
            write!(&mut fasta_writer_primary, ">{}\n{}\n", header, seq).unwrap();
        }
    }
}