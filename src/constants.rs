use block_aligner;

//pub const ID_CUTOFF: f64 = 0.999;
pub const FORWARD_READ_SAFE_SEARCH_CUTOFF: usize = 10000;
//foward_read_safe_search_cutoff: usize = 10000;
pub const MAX_GAP_CHAINING: usize = 200;
pub const SAMPLING_RATE_COV: usize = 10;
pub const MINIMIZER_END_NTH_COV: usize = 30;
pub const MINIMIZER_END_NTH_OVERLAP: usize = 20;
pub const QUANTILE_UNITIG_WEIGHT: f64 = 0.50;
pub const MID_BASE_THRESHOLD_READ: u8 = 23; // 98%
pub const MID_BASE_THRESHOLD_INITIAL: u8 = 10; // 90%
pub const MAX_BUBBLE_UNITIGS_FINAL_STAGE: usize = 5;
pub const TS_DASHES_BLANK_COLONS_DOT_BLANK: &str = "%Y-%m-%d %H:%M:%S%.3f";
pub const MIN_CHAIN_SCORE_COMPARE: i32 = 150;
pub const MIN_READ_LENGTH: usize = 1000;
pub const ENDPOINT_MAPPING_FUZZ: u32 = 200;
pub const LONG_OVERLAP_LENGTH: usize = 20_000;
// seed with 42 and 31 0s
pub const RNG_SEED: [u8; 32] = [
    42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];
pub const PSEUDOCOUNT: f64 = 3.;
pub const ID_THRESHOLD_ITERS: usize = 3;
//pub const IDENTITY_THRESHOLDS: [f64; ID_THRESHOLD_ITERS] = [0.995, 0.9975, 1.0];
pub const IDENTITY_THRESHOLDS: [f64; ID_THRESHOLD_ITERS] = [0.99, 0.9975, 1.0];
//pub const COV_MULTI_WEIGHTS: [f64; ID_THRESHOLD_ITERS] = [0.0, 0.0, 1.0];
pub const COV_MULTI_WEIGHTS: [f64; ID_THRESHOLD_ITERS] = [0.333, 0.333, 0.333];
/// Minimum spacing kept between rungs of the identity ladder (see `identity_ladder`).
pub const MIN_LADDER_GAP: f64 = 0.0005;

/// The DP1/DP2/DP3 identity ladder, with the strict (top) rung taken from
/// `--snpmer-threshold-strict` instead of being pinned at 100%. `strict_frac` is a
/// fraction (0.998), not a percentage.
///
/// With the default `--snpmer-threshold-strict 100` this returns `IDENTITY_THRESHOLDS`
/// unchanged. Lowering it also pulls down any rung that would otherwise sit above the new
/// top, because `sufficient_gap` in `map_processing` assumes the ladder is increasing.
pub fn identity_ladder(strict_frac: f64) -> [f64; ID_THRESHOLD_ITERS] {
    let mut ladder = IDENTITY_THRESHOLDS;
    ladder[ID_THRESHOLD_ITERS - 1] = strict_frac;
    for i in (0..ID_THRESHOLD_ITERS - 1).rev() {
        ladder[i] = ladder[i].min(ladder[i + 1] - MIN_LADDER_GAP);
    }
    ladder
}
pub const MIN_COV_READ: usize = 5;
pub const MIN_COV_READ_SMALL: usize = 3;
pub const SPECIFICITY_THRESHOLD: usize = 20;
pub const SUB_MATRIX: block_aligner::scores::NucMatrix =
    block_aligner::scores::NucMatrix::new_simple(1, -2);
pub const GAPS: block_aligner::scores::Gaps = block_aligner::scores::Gaps {
    open: -3,
    extend: -2,
};
pub const MIN_BLOCK_SIZE: usize = 32;
pub const MAX_BLOCK_SIZE: usize = 128;
pub const GAPS_LAX_INDEL: block_aligner::scores::Gaps = block_aligner::scores::Gaps {
    open: -3,
    extend: -1,
};
pub const MAX_OL_POLISHING: usize = 75;
pub const READ_BLOCK_SIZE_FOR_COVERAGE: usize = 50_000;
pub const OVERLAP_HANG_LENGTH: usize = 350;

//At most 1/20 k-mers are snpmers.
pub const MAX_FRACTION_OF_SNPMERS_IN_READ: f64 = 1. / 20.;
pub const SUPP_ALIGNMENT_SCORE_THRESHOLD: i32 = 2000;

pub const POLISHED_CONTIGS_NAME: &str = "initial_polished.fa";
pub const CIRC_STRICT_STRING: &str = "circular-yes";
pub const CIRC_LAX_STRING: &str = "circular-possibly";
pub const USE_SOLID_KMERS: bool = false;

pub const MAX_KMER_COUNT_IN_READ: usize = 500;

// A split k-mer is kept for SNPmer detection only if it was seen at least this many times in
// total and at least this many times on each strand. Shared by the in-memory counter and the
// on-disk counter (`myloasm-kmc`).
pub const MIN_SPLIT_KMER_COUNT_TOTAL: u32 = 3;
pub const MIN_SPLIT_KMER_COUNT_PER_STRAND: u32 = 1;
pub const MAX_MULTIPLICITY_KMER: usize = MAX_KMER_COUNT_IN_READ;
pub const QUALITY_SEQ_BIN: usize = 4;

pub const MINIMUM_MINIMIZER_FRACTION: f64 = 0.10;

pub const MAGIC_EXIST_STRING: &str = "exist";

pub const SAMPLES: usize = 20;
//pub const BEAM_STEPS: usize = 10;
pub const BEAM_STEPS: usize = 10;
pub const SAFE_LENGTH_BACK: usize = 300_000;
pub const MAX_LENGTH_SEARCH: usize = 1_000_000;

pub const MAX_ALLOWABLE_SNPMER_ERROR_MISC: usize = 2;
pub const MAX_ALLOWABLE_SNPMER_ERROR_DIVIDER: usize = 200;

pub const DEDUP_SNPMERS: bool = false;

pub const HUFFMAN_SAMPLE_MAX: usize = 500_000;

pub const CLI_HEADINGS: [&str; 5] = [
    "Technology Presets",
    "Basic Algorithmic Parameters",
    "Miscellaneous Options",
    "Overlap Parameters (advanced)",
    "Graph Parameters (advanced)",
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_ladder_default_is_unchanged() {
        // The default --snpmer-threshold-strict is 100, which must leave the ladder
        // (and therefore DP1/DP2/DP3 and every downstream threshold) exactly as it was.
        assert_eq!(identity_ladder(1.0), IDENTITY_THRESHOLDS);
    }

    #[test]
    fn identity_ladder_stays_increasing() {
        for strict in [1.0, 0.998, 0.997, 0.995, 0.99] {
            let ladder = identity_ladder(strict);
            assert_eq!(ladder[ID_THRESHOLD_ITERS - 1], strict);
            assert!(ladder.windows(2).all(|w| w[0] < w[1]), "{:?}", ladder);
        }
    }
}
