use minimap2;
use std::sync::Arc;
use crate::constants::*;
use crate::polishing::consensus::HomopolymerCompressedSeq;
use std::sync::Mutex;
use std::thread;
use gzp::{deflate::Gzip, ZBuilder};
use crossbeam_channel::unbounded;
use rayon::prelude::*;
use crate::types::*;
use crate::unitig::*;
use crate::cli::*;
use crate::kmer_comp::*;
use crate::polishing::alignment;
use crate::polishing::consensus;
use minimap2::Aligner;
use rust_htslib::bam::{Header, HeaderView};
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::Writer;
use bio_seq::prelude::*;
use fxhash::FxHashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::io::Write;
use std::path::PathBuf;

pub fn read_fastq_and_polish(final_graph: UnitigGraph, mut twin_reads: Vec<TwinRead>, args: &Cli, read_files: &Vec<PathBuf>){
    let fasta_out_path = Path::new(args.output_dir.as_str()).join("final_contigs.fa");
    let mut fasta_writer = BufWriter::new(File::create(fasta_out_path).unwrap());
    let mut read_id_to_map_loc: FxHashMap<String, Vec<(NodeIndex, u32, u32, u32, u32, bool)>> = FxHashMap::default();
    let mut contig_index_to_base_seq: FxHashMap<NodeIndex, Seq<Dna>> = FxHashMap::default();
    let mut contig_index_to_consensus_builders: FxHashMap<NodeIndex, Mutex<consensus::ConsensusBuilder>> = FxHashMap::default();

    log::debug!("Clearing reads...");
    for read in twin_reads.iter_mut(){
        read.clear();
    }

    log::debug!("Setting up alignments...");
    final_graph.nodes.iter().for_each(|(contig_index, contig)| {
        contig.mapping_boundaries().iter().for_each(|interval| {
            let ol = &interval.val;
            let reference_tr = &twin_reads[ol.query_id as usize];
            let read_base_id = reference_tr.base_id.clone();
            let q_start = ol.query_range.0 + reference_tr.split_start;
            let q_stop = ol.query_range.1 + reference_tr.split_start;
            let r_start = interval.start;
            let r_stop = interval.stop;
            read_id_to_map_loc.entry(read_base_id).or_insert(Vec::new()).push((*contig_index, q_start, q_stop, r_start, r_stop,ol.reverse));
            contig_index_to_base_seq.insert(*contig_index, contig.base_seq().clone());
            if contig_index_to_consensus_builders.get(contig_index).is_none(){
                let consensus_builder = consensus::ConsensusBuilder::new(contig.base_seq().len(), false);
                contig_index_to_consensus_builders.insert(*contig_index, Mutex::new(consensus_builder));
            }
        });
    });


    let files_owned = read_files.clone();
    let hpc = args.homopolymer_compression;
    let arc_locations = Arc::new(read_id_to_map_loc);
    let arc_contig_seqs = Arc::new(contig_index_to_base_seq);
    let arc_consensus_builders = Arc::new(contig_index_to_consensus_builders);

    log::debug!("Re-reading reads and aligning...");
    for fastq_file in files_owned{
        let (mut tx, rx) = spmc::channel();
        thread::spawn(move || {
            let mut reader = needletail::parse_fastx_file(fastq_file).expect("valid path");
            while let Some(record) = reader.next() {
                let rec = record.expect("Error reading record");
                let seq;
                if hpc{
                    seq = homopolymer_compression(rec.seq().to_vec());
                } else {
                    seq = rec.seq().to_vec();
                }
                if seq.len() < MIN_READ_LENGTH{
                    continue;
                }
                let id = String::from_utf8_lossy(rec.id()).to_string();
                if let Some(qualities) = rec.qual(){
                    tx.send((seq, Some(qualities.to_vec()), id)).unwrap();
                }
                else{
                    tx.send((seq, None, id)).unwrap();
                }
            }
        });

        let mut handles = Vec::new();
        for _ in 0..args.threads{
            let rx = rx.clone();
            let location_map = Arc::clone(&arc_locations);
            let contig_seqs = Arc::clone(&arc_contig_seqs);
            let consensus_builders = Arc::clone(&arc_consensus_builders);
            handles.push(thread::spawn(move || {
                loop{
                    match rx.recv() {
                        Ok(msg) => {
                            let seq = msg.0;
                            let seqlen = seq.len() as u32;
                            let qualities = msg.1.unwrap();
                            let qualities = quality_pool(qualities);
                            //Take the quality as the minimum over the 3 left and 3 right bases
                            let id = msg.2;
                            let locations = location_map.as_ref().get(&id);
                            if let Some(locations) = locations{
                                for (contig_index, q_start, q_stop, r_start, r_stop, reverse) in locations.iter(){
                                    let contig_seq = contig_seqs.as_ref().get(contig_index).unwrap();
                                    let start_read;
                                    let stop_read;
                                    let query_u8;
                                    let query_qualities_u8;
                                    if *reverse{
                                        //[0,1,2] -> [2,1,0]; seqlen - pos - 1
                                        start_read = seqlen - *q_stop - 1;
                                        stop_read = seqlen - *q_start - 1;
                                        query_u8 = revcomp_u8(&seq);
                                        query_qualities_u8 = qualities.iter().rev().cloned().collect();
                                    }
                                    else{
                                        start_read = *q_start;
                                        stop_read = *q_stop;
                                        //for now... inefficent but it's probably fine
                                        query_u8 = seq.clone();
                                        query_qualities_u8 = qualities.clone();
                                    }
                                    let start_read = start_read as usize;
                                    let stop_read = stop_read as usize;

                                    let contig_slice = &contig_seq[*r_start as usize..*r_stop as usize];
                                    let contig_u8_ref = dna_slice_to_u8(contig_slice);
                                    let query_u8_ref = &query_u8[start_read..stop_read+1];
                                    let query_qualities_ref = query_qualities_u8[start_read..stop_read+1].to_vec();

                                    let cigar = alignment::align_seq_to_ref_slice(&contig_u8_ref, query_u8_ref, &GAPS);
                                    let seq_wrapper = HomopolymerCompressedSeq::new(&query_u8_ref, query_qualities_ref, false);
                                    let mut consensus_builder = consensus_builders.as_ref().get(contig_index).unwrap().lock().unwrap();
                                    consensus_builder.process_alignment(cigar, &seq_wrapper, *r_start as usize, 0);
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
    }

    log::info!("Finished alignments, generating consensus...");

    let consensuses = Arc::try_unwrap(arc_consensus_builders).unwrap();
    for (contig_index, consensus_builder) in consensuses{
        let consensus_builder = consensus_builder.into_inner().unwrap();
        let node_index = final_graph.nodes[&contig_index].node_id;
        let consensus_seq = consensus_builder.produce_consensus();
        write!(&mut fasta_writer, ">u{}\n", node_index).unwrap();
        write!(&mut fasta_writer, "{}\n", std::str::from_utf8(&consensus_seq).unwrap()).unwrap();
    }
}

fn dna_slice_to_u8(slice: &SeqSlice<Dna>) -> Vec<u8>{
    slice.iter().map(|x| x.to_char().to_ascii_uppercase() as u8).collect()
}

fn revcomp_u8(seq: &Vec<u8>) -> Vec<u8>{
    seq.iter().rev().map(|x| match x{
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N'
    }).collect()
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

fn quality_pool(qualities: Vec<u8>) -> Vec<u8>{
    let pool_width = 5;
    let mut pool = Vec::new();
    for i in 0..qualities.len(){
        if i > pool_width/2 && i < qualities.len() - pool_width/2{
            let mut min = 255;
            for j in i-pool_width/2..i+pool_width/2{
                if qualities[j] < min{
                    min = qualities[j];
                }
            }
            pool.push(min);
        }
        else{
            pool.push(qualities[i]);
        }
    }
    return pool;
}