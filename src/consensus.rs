use minimap2;
use std::sync::Mutex;
use rayon::prelude::*;
use crate::types::*;
use crate::unitig::*;
use crate::cli::*;
use minimap2::Aligner;
use rust_htslib::bam::{Header, HeaderView};
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::Writer;
use bio_seq::prelude::*;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::io::Write;

fn create_bam_header(sequences: Vec<(String, u32)>) -> Header {
    // Create a new header
    let mut header = Header::new();

    // Add a HD (header) line indicating this is a BAM/SAM file
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", &"1.6")  // SAM/BAM specification version
        .push_tag(b"SO", &"unknown");  // Sorting order
    header.push_record(&hd);

    // Add sequence information (SQ lines)
    for (seq_name, length) in sequences {
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", &seq_name).push_tag(b"LN", &length);   // Sequence length
        header.push_record(&sq);
    }

    header
}

pub fn outer_consensus(final_graph: &UnitigGraph, reads: &Vec<TwinRead>, args: &Cli) {

    let fasta_out = Path::new(args.output_dir.as_str()).join("reads_out.fasta");
    let fasta_writer = Mutex::new(BufWriter::new(File::create(fasta_out).unwrap()));
    let contig_vec = final_graph.nodes.values().collect::<Vec<&UnitigNode>>();
    let bam_output = Path::new(args.output_dir.as_str()).join("map.bam");
    let ids_and_lens_contigs = contig_vec.iter()
    .map(|contig| (format!("u{}", contig.node_id), contig.base_seq().len() as u32))
    .collect::<Vec<(String, u32)>>();
    let header = create_bam_header(ids_and_lens_contigs);
    let mut bam_writer = Writer::from_path(&bam_output, &header, rust_htslib::bam::Format::Bam).unwrap();
    bam_writer.set_threads(args.threads as usize).unwrap();
    let bam_writer_lock = Mutex::new(bam_writer);

    //println!("Size of bam_writer: {}", size_of_val(&bam_writer));
    
    //Need to create new fasta for split reads. 
    //Bam file for each contig
    //Merge all and then use racon for now. 
    contig_vec.iter().for_each(|contig| {
        let contig_seq: Vec<u8> = contig.base_seq().iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect();
        if contig_seq.len() == 0{
            return;
        }
        let aligner = Aligner::builder().with_cigar().preset(minimap2::Preset::MapOnt);
        let aligner_with_index = aligner.with_seq_and_id(&contig_seq, format!("u{}", contig.node_id).as_bytes()).unwrap();
        //let mut header = Header::new();
        //aligner_with_index.populate_header(&mut header);
        contig.mapped_indices().par_iter().for_each(|read_id| {
            let header_view = HeaderView::from_header(&header);
            let read_seq_u8 = reads[*read_id].dna_seq.iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect::<Vec<u8>>();
            if read_seq_u8.len() == 0{
                return;
            }
            let read_name = format!("r{}", *read_id);
            {
                let mut writer = fasta_writer.lock().unwrap();
                write!(writer, ">r{}\n", *read_id).unwrap();
                write!(writer, "{}\n", std::str::from_utf8(&read_seq_u8).unwrap()).unwrap();
            }
            let read_name_u8 = read_name.as_bytes();
            let records = aligner_with_index.map(&read_seq_u8, true, false, None, None, Some(&read_name_u8));
            //let records = aligner_with_index.map_to_sam(&read_seq_u8, None, Some(read_name_u8), &header_view, None, None).unwrap();
            //TODO implement consensus
            return;
        });
    });
}

fn _write_paf(writer: &Mutex<BufWriter<File>>, record: &minimap2::Mapping){
    let mut writer = writer.lock().unwrap();
    write!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        record.query_name.as_ref().unwrap(), record.query_len.unwrap(), record.query_start, record.query_end, 
        record.strand, record.target_name.as_ref().unwrap(), record.target_len, record.target_start, record.target_end,
        record.match_len, record.block_len, record.mapq).unwrap();
}