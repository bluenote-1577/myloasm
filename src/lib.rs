pub mod cli;
pub mod constants;
pub mod graph;
pub mod kmc_reader;
pub mod kmer_comp;
pub mod map_processing;
pub mod mapping;
pub mod mphmap;
pub mod polishing;
pub mod polishing_mod;
pub mod seeding;
pub mod seq_parse;
pub mod skani_dereplicate;
pub mod small_genomes;
pub mod twin_graph;
pub mod types;
pub mod unitig;
pub mod unitig_utils;
pub mod utils;

//pub mod cbloom;
//
//#[cfg(target_arch = "x86_64")]
//pub mod avx2_seeding;
//#[cfg(target_arch = "x86_64")]
//pub mod avx2_chaining;

// Use of a mod or pub mod is not actually necessary.
pub mod built_info {
    // The file has been placed there by the build script.
    // include!(concat!(env!("OUT_DIR"), "/built.rs"));
}
