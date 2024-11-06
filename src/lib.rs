pub mod seeding;
pub mod types;
pub mod bam_parsing;
pub mod kmer_comp;
pub mod mapping;
pub mod graph;
pub mod twin_graph;
pub mod unitig;
pub mod cli;

#[cfg(target_arch = "x86_64")]
pub mod avx2_seeding;


