use crate::polishing::consensus2::PoaConsensusBuilder;
use crate::types::*;
use crate::unitig::*;
use crate::cli::*;
use rust_lapper::Interval;
use rust_htslib::bam::{Header, HeaderView};
use rust_htslib::bam::header::HeaderRecord;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::io::Write;

pub fn polish_assembly(final_graph: UnitigGraph, twin_reads: Vec<TwinRead>, args: &Cli){
    let fasta_out_path = Path::new(args.output_dir.as_str()).join("final_contigs_polished.fa");
    let mut fasta_writer = BufWriter::new(File::create(fasta_out_path).unwrap());
    let mut total_count = 0;
    let mut reset_count = 0;
    let num_passing_nodes = final_graph.nodes.iter().filter(|(_, contig)| UnitigGraph::unitig_pass_filter(contig, args)).count();

    final_graph.nodes.iter().for_each(|(_, contig)| {
        if !UnitigGraph::unitig_pass_filter(contig, args){
            return
        }
        log::trace!("Processing alignments for u{} ...", contig.node_id);
        let mut poa_cons_builder = PoaConsensusBuilder::new(contig.base_seq().len(), format!("u{}", contig.node_id));
        poa_cons_builder.generate_breakpoints(300, 100);
        if contig.mapping_info.max_alignment_boundaries.is_none(){
            return
        }
        let mapping_boundaries = contig.mapping_info.max_alignment_boundaries.as_ref().unwrap().iter().collect::<Vec<&Interval<u32, SmallTwinOl>>>();
        poa_cons_builder.process_mapping_boundaries(&mapping_boundaries, &twin_reads);

        log::trace!("Starting POA consensus for u{} ...", contig.node_id);
        let cons = poa_cons_builder.spoa_blocks();
        let mut final_seq = Vec::new();
        for consensus in cons{
            final_seq.extend(consensus);
        }

        //Output user logging info at 10% intervals
        if reset_count == num_passing_nodes/10{
            log::info!("Polished {:.2}% of contigs...", (total_count as f64/num_passing_nodes as f64)*100.0);
            reset_count = 0;
        }

        write!(&mut fasta_writer, ">u{}\n", contig.node_id).unwrap();
        write!(&mut fasta_writer, "{}\n", std::str::from_utf8(&final_seq).unwrap()).unwrap();
        total_count +=1;
    });
}


fn _create_bam_header(sequences: Vec<(String, u32)>) -> Header {
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