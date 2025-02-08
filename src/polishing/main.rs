// use clap::Parser;
// use std::io::Write;

// use needletail::parse_fastx_file;
// use std::path::PathBuf;
// use block_aligner::{scan_block::*, scores::*};
// use polishing::alignment::*;
// use polishing::consensus::*;
// use polishing::polish_const::*;


// /// CLI tool to parse reference FASTA and read files
// #[derive(Parser, Debug)]
// #[clap(author, version, about)]
// struct Args {
//     /// Path to the reference FASTA file
//     #[clap(short, long)]
//     reference: PathBuf,

//     /// Paths to one or more read files (FASTA/FASTQ)
//     #[clap(short, long)]
//     query: Vec<PathBuf>,

//     #[clap(short, long)]
//     hpc: bool
// }

// fn main() -> Result<(), Box<dyn std::error::Error>> {
//     let args = Args::parse();
//     let mut padded_refs = vec![];
//     let mut refs = vec![];
//     let mut consensuses = vec![];
//     let gaps = Gaps { open: -5, extend: -4 };
//     let hpc = args.hpc;


//     // Parse reference sequence
//     println!("Processing reference file: {:?}", args.reference);
//     let mut reference_reader = parse_fastx_file(&args.reference)?;
//     while let Some(record) = reference_reader.next() {
//         let reference = record?;
//         println!(
//             "Reference sequence: id={}, length={}",
//             String::from_utf8_lossy(reference.id()),
//             reference.seq().len()
//         );
        
//         // Access the sequence if needed:
//         let sequence = reference.seq();
//         let quals;
//         if reference.qual().is_none() {
//             quals = vec![10; sequence.len()];
//         }
//         else{
//             quals = reference.qual().unwrap().to_vec();
//         }
//         let hpc_seq = HomopolymerCompressedSeq::new(&sequence, quals, hpc);
//         let padded_ref = PaddedBytes::from_bytes::<NucMatrix>(&hpc_seq.seq, MAX_BLOCK_SIZE);
//         padded_refs.push(padded_ref);
//         consensuses.push(ConsensusBuilder::new(hpc_seq.seq.len(), hpc));
//         refs.push(hpc_seq.seq);
//     }

//     // Process each read file
//     for read_file in args.query{
//         println!("\nProcessing read file: {:?}", read_file);
//         let mut read_reader = parse_fastx_file(&read_file)?;
        
//         let mut read_count = 0;
//         while let Some(record) = read_reader.next() {
//             let read = record?;
//             read_count += 1;
            
//             // Access read information if needed:
//             // let id = read.id();
//             let seq = read.seq();
//             let quals = read.qual().unwrap();
//             let hpc_seq = HomopolymerCompressedSeq::new_slice(&seq, &quals, hpc);
//             let padded_seq = PaddedBytes::from_bytes::<NucMatrix>(&hpc_seq.seq, MAX_BLOCK_SIZE);
//             // let qual = read.qual(); // Only available for FASTQ
//             //align_seq_to_ref_end_to_end(&padded_refs[0], &padded_seq, &hpc_seq, &gaps, &mut consensuses[0]);
//             align_seq_to_ref_slice(&refs[0], 0, refs[0].len(), &seq, 0, seq.len(), &hpc_seq, &gaps, &mut consensuses[0]);
//         }
//         println!("Total reads processed: {}", read_count);
//     }

//     for consensus in consensuses{
//         for (i,op) in consensus.consensus.iter().enumerate(){
//             println!("POS:{}, INS: {:?}, DEL: {}, A {} , C {}, G {}, T {}, HPC_MATCH: {:?}, PREV_INS: {}, NONPREV_INS: {}", i, &op.insertions, op.deletion_count, op.a_count, op.c_count, op.g_count, op.t_count, op.homopolymer_lengths, op.prev_ins_weight, op.prev_nonins_weight);
//         }
//         let file = std::fs::File::create("consensus.fasta")?;
//         let mut bufwriter = std::io::BufWriter::new(file);
//         let consensus_seq = consensus.produce_consensus();
//         write!(&mut bufwriter, ">test\n{}", String::from_utf8_lossy(&consensus_seq)).unwrap();
//     }

//     Ok(())
// }
