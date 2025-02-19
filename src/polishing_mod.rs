use minimap2;
//use crate::polishing::consensus::HomopolymerCompressedSeq;
use crate::polishing::consensus2::PoaConsensusBuilder;
use std::sync::Mutex;
use std::thread;
use gzp::{deflate::Gzip, ZBuilder};
use crossbeam_channel::unbounded;
use rayon::prelude::*;
use crate::types::*;
use crate::unitig::*;
use crate::cli::*;
//use crate::polishing::consensus;
use rust_lapper::Interval;
use minimap2::Aligner;
use rust_htslib::bam::{Header, HeaderView};
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::Writer;
use bio_seq::prelude::*;
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
        let mut poa_cons_builder = PoaConsensusBuilder::new(contig.base_seq().len());
        poa_cons_builder.generate_breakpoints(300, 100);
        let mapping_boundaries = contig.mapping_boundaries().iter().collect::<Vec<&Interval<u32, SmallTwinOl>>>();
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

pub fn write_to_paf(final_graph: &UnitigGraph, reads: &Vec<TwinRead>, args: &Cli){
    let fasta_out_path = Path::new(args.output_dir.as_str()).join("reads_out.fa.gz");
    let paf_out_path = Path::new(args.output_dir.as_str()).join("final_mapping.paf.gz");
    let paf_writer = BufWriter::new(File::create(paf_out_path).unwrap());
    let fasta_writer = BufWriter::new(File::create(fasta_out_path).unwrap());
    //let fasta_encoder = Mutex::new(GzEncoder::new(fasta_writer, Compression::default()));
    //let paf_encoder = Mutex::new(GzEncoder::new(BufWriter::new(File::create(paf_out_path).unwrap()), Compression::default()));
    let contig_vec = final_graph.nodes.values().collect::<Vec<&UnitigNode>>();

    let (tx, rx) = unbounded();
    let (tx2, rx2) = unbounded();

    contig_vec.iter().for_each(|contig| {
        let contig_seq: Vec<u8> = contig.base_seq().iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect();
        if contig_seq.len() == 0{
            return;
        }

        let aligner = Aligner::builder().preset(minimap2::Preset::MapPb);
        let aligner_with_index = aligner.with_seq_and_id(&contig_seq, format!("u{}", contig.node_id).as_bytes()).unwrap();

        contig.mapped_indices().par_iter().for_each(|read_id| {
            let read_seq_u8 = reads[*read_id].dna_seq.iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect::<Vec<u8>>();

            if read_seq_u8.len() == 0{
                return;
            }

            let read_name = format!("r{}", *read_id);
            let read_name_u8 = read_name.as_bytes().to_owned();
            let records = aligner_with_index.map(&read_seq_u8, true, false, None, None, Some(&read_name_u8));

            tx.send((read_name, read_seq_u8)).unwrap();
            tx2.send(records).unwrap();
            return;
        });
    });

    drop(tx);
    drop(tx2);
    let num_threads_fasta = args.threads/4;
    let num_threads_paf = args.threads/6;

    thread::spawn(move || {
        let mut parz_fasta = ZBuilder::<Gzip, _>::new().num_threads(num_threads_fasta).from_writer(fasta_writer);
        while let Ok((read_name, read_seq_u8)) = rx.recv() {
            write!(parz_fasta, ">{}\n", read_name).unwrap();
            write!(parz_fasta, "{}\n", std::str::from_utf8(&read_seq_u8).unwrap()).unwrap();
        }
    }).join().unwrap();

    thread::spawn(move || {
        let mut parz_paf = ZBuilder::<Gzip, _>::new().num_threads(num_threads_paf).from_writer(paf_writer);
        while let Ok(records) = rx2.recv() {
            write_paf(&mut parz_paf, records.unwrap());
        }
    }).join().unwrap();
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
    let _bam_writer_lock = Mutex::new(bam_writer);

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
            let _header_view = HeaderView::from_header(&header);
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
            let _records = aligner_with_index.map(&read_seq_u8, true, false, None, None, Some(&read_name_u8));
            //let records = aligner_with_index.map_to_sam(&read_seq_u8, None, Some(read_name_u8), &header_view, None, None).unwrap();
            //TODO implement consensus
            return;
        });
    });
}

fn write_paf<T>(writer: &mut Box<T>, records: Vec<minimap2::Mapping>)
where 
    T: Write + ?Sized{
    for record in records {
        write!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            record.query_name.as_ref().unwrap(), record.query_len.unwrap(), record.query_start, record.query_end, 
            record.strand, record.target_name.as_ref().unwrap(), record.target_len, record.target_start, record.target_end,
            record.match_len, record.block_len, record.mapq).unwrap();
    }
}

