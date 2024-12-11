pub mod seeding;
pub mod types;
pub mod kmer_comp;
pub mod mapping;
pub mod twin_graph;
pub mod graph;
pub mod unitig;
pub mod cli;
pub mod constants;
pub mod cbloom;
pub mod consensus;
pub mod seq_parse;

#[cfg(target_arch = "x86_64")]
pub mod avx2_seeding;
#[cfg(target_arch = "x86_64")]
pub mod avx2_chaining;


