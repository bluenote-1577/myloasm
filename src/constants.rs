use block_aligner;

//pub const ID_CUTOFF: f64 = 0.999;
pub const FORWARD_READ_SAFE_SEARCH_CUTOFF: usize = 10000;
//foward_read_safe_search_cutoff: usize = 10000;
pub const MAX_GAP_CHAINING: usize = 200;
pub const SAMPLING_RATE_COV: usize = 10;
pub const MINIMIZER_END_NTH_COV: usize = 30;
pub const MINIMIZER_END_NTH_OVERLAP: usize = 30;
pub const QUANTILE_UNITIG_WEIGHT: f64 = 0.50;
pub const MID_BASE_THRESHOLD_READ: u8 = 23; // 98%
pub const MID_BASE_THRESHOLD_INITIAL: u8 = 10; // 90%
pub const MAX_BUBBLE_UNITIGS_FINAL_STAGE: usize = 5;
pub const TS_DASHES_BLANK_COLONS_DOT_BLANK: &str = "%Y-%m-%d %H:%M:%S%.3f";
pub const MIN_CHAIN_SCORE_COMPARE: i32 = 150;
pub const MIN_READ_LENGTH: usize = 1000;
pub const ENDPOINT_MAPPING_FUZZ : u32 = 200;
// seed with 42 and 31 0s
pub const RNG_SEED: [u8; 32] = [42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
pub const PSEUDOCOUNT: f64 = 3.;
pub const ID_THRESHOLD_ITERS: usize = 3;
//pub const IDENTITY_THRESHOLDS: [f64; ID_THRESHOLD_ITERS] = [0.995, 0.9975, 1.0];
pub const IDENTITY_THRESHOLDS: [f64; ID_THRESHOLD_ITERS] = [0.99, 0.9975, 1.0];
//pub const COV_MULTI_WEIGHTS: [f64; ID_THRESHOLD_ITERS] = [0.0, 0.0, 1.0];
pub const COV_MULTI_WEIGHTS: [f64; ID_THRESHOLD_ITERS] = [0.333, 0.333, 0.333];
pub const MIN_COV_READ: usize = 5;
pub const MIN_COV_READ_SMALL: usize = 3;
pub const SUB_MATRIX: block_aligner::scores::NucMatrix = block_aligner::scores::NucMatrix::new_simple(3, -2);
pub const GAPS: block_aligner::scores::Gaps = block_aligner::scores::Gaps { open: -3, extend: -2 };
pub const MIN_BLOCK_SIZE: usize = 32;
pub const MAX_BLOCK_SIZE: usize = 128;
pub const GAPS_LAX_INDEL: block_aligner::scores::Gaps = block_aligner::scores::Gaps { open: -3, extend: -1 };
pub const MAX_OL_POLISHING: usize = 75;
pub const READ_BLOCK_SIZE_FOR_COVERAGE: usize = 50_000;
pub const OVERLAP_HANG_LENGTH: usize = 750;

//At most 1/20 k-mers are snpmers. 
pub const MAX_FRACTION_OF_SNPMERS_IN_READ: f64 = 1./20.;
pub const SUPP_ALIGNMENT_SCORE_THRESHOLD: i32 = 2000;

pub const POLISHED_CONTIGS_NAME: &str = "initial_polished.fa";
pub const USE_SOLID_KMERS: bool = false;

pub const MAX_KMER_COUNT_IN_READ: usize = 500;
pub const MAX_MULTIPLICITY_KMER: usize = MAX_KMER_COUNT_IN_READ;
pub const QUALITY_SEQ_BIN: usize = 4;

pub const MINIMUM_MINIMIZER_FRACTION: f64 = 0.10;

pub const MAGIC_EXIST_STRING: &str = "exist";

pub const SAMPLES: usize = 20;
pub const BEAM_STEPS: usize = 10;
pub const SAFE_LENGTH_BACK : usize = 300_000;
pub const MAX_LENGTH_SEARCH: usize = 1_000_000;

