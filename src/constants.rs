use block_aligner;

//pub const ID_CUTOFF: f64 = 0.999;
pub const FORWARD_READ_SAFE_SEARCH_CUTOFF: usize = 10000;
//foward_read_safe_search_cutoff: usize = 10000;
pub const MAX_GAP_CHAINING: usize = 200;
pub const SAMPLING_RATE_COV: usize = 10;
pub const MINIMIZER_END_NTH_COV: usize = 20;
pub const QUANTILE_UNITIG_WEIGHT: f64 = 0.50;
pub const MID_BASE_THRESHOLD_READ: u8 = 17; // 98%
pub const MID_BASE_THRESHOLD_INITIAL: u8 = 10; // 90%
pub const MAX_BUBBLE_UNITIGS_FINAL_STAGE: usize = 5;
pub const TS_DASHES_BLANK_COLONS_DOT_BLANK: &str = "%Y-%m-%d %H:%M:%S%.3f";
pub const MIN_CHAIN_SCORE_COMPARE: i32 = 50;
pub const MIN_READ_LENGTH: usize = 1000;
pub const ENDPOINT_MAPPING_FUZZ : u32 = 200;
// seed with 42 and 31 0s
pub const RNG_SEED: [u8; 32] = [42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
pub const PSEUDOCOUNT: f64 = 3.;
pub const ID_THRESHOLD_ITERS: usize = 3;
pub const IDENTITY_THRESHOLDS: [f64; ID_THRESHOLD_ITERS] = [0.995, 0.9975, 1.0];
//pub const COV_MULTI_WEIGHTS: [f64; ID_THRESHOLD_ITERS] = [0.0, 0.0, 1.0];
pub const COV_MULTI_WEIGHTS: [f64; ID_THRESHOLD_ITERS] = [0.333, 0.333, 0.333];
pub const MIN_COV_READ: usize = 4;
pub const SUB_MATRIX: block_aligner::scores::NucMatrix = block_aligner::scores::NucMatrix::new_simple(3, -4);
pub const GAPS: block_aligner::scores::Gaps = block_aligner::scores::Gaps { open: -3, extend: -2 };
pub const MIN_BLOCK_SIZE: usize = 16;
pub const MAX_BLOCK_SIZE: usize = 64;