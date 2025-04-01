use crate::constants::POLISHED_CONTIGS_NAME;
use crate::graph::GraphNode;
use crate::polishing::consensus2::PoaConsensusBuilder;
use crate::polishing::consensus2::join_circular_ends;
use crate::types::*;
use crate::unitig::*;
use crate::small_genomes;
use crate::cli::*;
use rust_lapper::Interval;
use rust_htslib::bam::{Header, HeaderView};
use rust_htslib::bam::header::HeaderRecord;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::io::Write;

pub fn polish_assembly(final_graph: UnitigGraph, twin_reads: Vec<TwinRead>, args: &Cli){
    let fasta_out_path = Path::new(args.output_dir.as_str()).join(POLISHED_CONTIGS_NAME);
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

        let mut circ_string = String::new();

        if contig.is_circular(){
            let edge_id = contig.in_edges()[0];
            let edge = final_graph.edges[edge_id].as_ref().unwrap();
            join_circular_ends(&mut final_seq, edge.overlap.overlap_len_bases, edge.overlap.hang1, edge.overlap.hang2, &format!("u{}", &contig.node_id), args); 
            circ_string = "circular".to_string();
        }

        //small_genomes::fix_multimers(&contig, &twin_reads, &mut final_seq, args);

        //Output user logging info at 10% intervals
        if reset_count == num_passing_nodes/10{
            log::info!("Polished {:.0}% of contigs...", 
            //(total_count as f64/num_passing_nodes as f64)*100.0); set to 10,20,30...%
            (total_count as f64/num_passing_nodes as f64)*100.0);
            reset_count = 0;
        }

        write!(&mut fasta_writer, ">u{} len:{} depth:{:?} {}\n", contig.node_id, final_seq.len(), contig.min_read_depth_multi.unwrap(), circ_string).unwrap();
        write!(&mut fasta_writer, "{}\n", std::str::from_utf8(&final_seq).unwrap()).unwrap();
        total_count +=1;
        reset_count += 1;
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

