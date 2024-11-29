use minimap2;
use std::sync::Mutex;
use rayon::prelude::*;
use crate::types::*;
use crate::unitig::*;
use crate::cli::*;
use minimap2::Aligner;
use rust_htslib::bam::{Header, HeaderView};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::Writer;
use bio_seq::prelude::*;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::io::Write;

pub fn outer_consensus(final_graph: &UnitigGraph, reads: &Vec<TwinRead>, args: &Cli) {

    let aligner = Aligner::builder().with_cigar().preset(minimap2::Preset::MapOnt);
    let fasta_out = Path::new(args.output_dir.as_str()).join("reads_out.fasta");
    let fasta_writer = Mutex::new(BufWriter::new(File::create(fasta_out).unwrap()));
    let contig_vec = final_graph.nodes.values().collect::<Vec<&UnitigNode>>();

    //Need to create new fasta for split reads. 
    //Bam file for each contig
    //Merge all and then use racon for now. 
    contig_vec.par_iter().for_each(|contig| {
        let bam_output = Path::new(args.output_dir.as_str()).join("temp").join(format!("contig_{}.bam", contig.node_id));
        let contig_seq: Vec<u8> = contig.base_seq().iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect();
        let aligner_with_index = aligner.clone().with_seq_and_id(&contig_seq, format!("u{}", contig.node_id).as_bytes()).unwrap();
        let mut header = Header::new();
        aligner_with_index.populate_header(&mut header);
        let header_view = HeaderView::from_header(&header);
        let mut contig_records = vec![];
        let mut out_contig = Writer::from_path(&bam_output, &header, rust_htslib::bam::Format::Bam).unwrap();
        for read_id in contig.mapped_indices(){
            let read_seq_u8 = reads[*read_id].dna_seq.iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect::<Vec<u8>>();
            let read_name = format!("r{}", *read_id);
            {
                let mut writer = fasta_writer.lock().unwrap();
                write!(writer, ">r{}\n", *read_id).unwrap();
                write!(writer, "{}\n", std::str::from_utf8(&read_seq_u8).unwrap()).unwrap();
            }
            let records = aligner_with_index.map_to_sam(&read_seq_u8, None, Some(&Vec::from(read_name)), &header_view, None, None).unwrap();
            contig_records.extend(records);
            if contig_records.len() > 1000{
                for record in contig_records.iter(){
                    out_contig.write(&record).unwrap();
                }
                contig_records.clear();
            }
        }
        for record in contig_records{
            out_contig.write(&record).unwrap();
        }
    });

}