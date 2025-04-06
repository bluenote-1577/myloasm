use clap::{Parser, ValueEnum};

use crate::constants::{IDENTITY_THRESHOLDS, ID_THRESHOLD_ITERS};

#[derive(Parser, Debug)]
#[command(
    name = "myloasm",
    about = "myloasm - high-resolution long-read metagenomic assembly with even low-fidelity reads.",
    version,
    author
)]

#[derive(Default, Clone)]
pub struct Cli {
    /// Input FASTQ files
    #[arg(num_args = 1..)]
    pub input_files: Vec<String>,

    /// K-mer size (must be odd)
    #[arg(short, long, default_value = "21")]
    pub kmer_size: usize,

    /// Compression ratio
    #[arg(short, long, default_value = "9")]
    pub c: usize,

    /// Number of threads to use for processing
    #[arg(short, long, default_value = "20")]
    pub threads: usize,

    /// Enable homopolymer compression
    #[arg(long, hide=true)]
    pub homopolymer_compression: bool,

    /// Output directory for results
    #[arg(short, long, default_value = "output")]
    pub output_dir: String,

    /// Verbosity level. Written to the .log file in output directory.
    #[arg(short, long, value_enum, default_value = "debug")]
    pub log_level: LogLevel,

    /// Length of tip to remove
    #[arg(long, default_value_t = 20000, help_heading = "Graph Parameters")]
    pub tip_length_cutoff: usize,

    /// Number of reads in tips to remove
    #[arg(long, default_value_t = 3, help_heading = "Graph Parameters")]
    pub tip_read_cutoff: usize,

    // No polishing
    #[arg(long, default_value_t=false, help_heading = "Overlap Parameters")]
    pub no_polish: bool,

    #[arg(long, default_value_t=false, help_heading = "Overlap Parameters")]
    pub no_snpmers: bool,

    /// Disallow reads with < this accuracy from Q-scores for the overlap step. 
    #[arg(long, default_value_t=90.)]
    pub quality_value_cutoff: f64,

    /// Snpmer identity threshold for containment
    #[arg(long, default_value_t=IDENTITY_THRESHOLDS[ID_THRESHOLD_ITERS - 1] * 100., help_heading = "Overlap Parameters")]
    pub snpmer_threshold_strict: f64,

    /// Mininum snpmer identity threshold for overlaps
    #[arg(long, default_value_t=IDENTITY_THRESHOLDS[0] * 100., help_heading = "Overlap Parameters")]
    pub snpmer_threshold_lax: f64,

    /// Error rate for snpmers for binomial test
    #[arg(long, default_value_t=0.025, help_heading = "Overlap Parameters")]
    pub snpmer_error_rate_lax: f64,

    /// Error rate for snpmers for binomial test
    #[arg(long, default_value_t=0.00, help_heading = "Overlap Parameters")]
    pub snpmer_error_rate_strict: f64,


    #[arg(long, default_value_t=30, help_heading = "Overlap Parameters")]
    pub contain_subsample_rate: usize,

    /// Cut overlaps with > (c * this) number of bases between minimizers
    #[arg(long, default_value_t=8., help_heading = "Overlap Parameters")]
    pub absolute_minimizer_cut_ratio: f64,

    /// Cut overlaps with > 5 times more bases between minimizers than the best overlap
    #[arg(long, default_value_t=5., help_heading = "Overlap Parameters")]
    pub relative_minimizer_cut_ratio: f64,

    /// Minimum overlap length for graph construction
    #[arg(long, default_value_t=500, help_heading = "Overlap Parameters")]
    pub min_ol: usize,

    /// Disables a SNPmer error overlap rescue heuristic during graph construction
    #[arg(long, help_heading = "Overlap Parameters")]
    pub disable_error_overlap_rescue: bool,

    /// Maximum bubble length to pop; keep alternates
    #[arg(long, default_value_t=500000, help_heading = "Graph Parameters")]
    pub max_bubble_threshold: usize,

    /// Small bubble length to pop; discard alternates
    #[arg(long, default_value_t=50000, help_heading = "Graph Parameters")]
    pub small_bubble_threshold: usize,

    /// Small bubble length to pop; discard alternates
    #[arg(long, default_value_t=1.0, help_heading = "Graph Parameters")]
    pub z_edge_threshold: f64,

    #[arg(long, default_value_t=200, help_heading = "Alignment Parameters")]
    pub maximal_end_fuzz: usize, 

    /// Bloom filter size in GB
    #[arg(short, long, default_value_t=10.)]
    pub bloom_filter_size: f64,

    /// Minimum number of reads in output contigs
    #[arg(long, default_value_t = 1)]
    pub min_reads_contig: usize,
    
    /// HiFi mode (--snpmer-threshold-strict 100 --snpmer-error-rate 0.001)
    #[arg(long, help_heading = "Preset Parameters")]
    pub hifi: bool,

    /// R9 (old nanopore) mode (--snpmer-error-rate 0.05)
    #[arg(long, help_heading = "Preset Parameters", hide=true)]
    pub r941: bool
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum LogLevel {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum Preset{
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

impl Default for LogLevel {
    fn default() -> Self {
        LogLevel::Debug
    }
}


impl Cli {
    pub fn log_level_filter(&self) -> log::LevelFilter {
        match self.log_level {
            LogLevel::Error => log::LevelFilter::Error,
            LogLevel::Warn => log::LevelFilter::Warn,
            LogLevel::Info => log::LevelFilter::Info,
            LogLevel::Debug => log::LevelFilter::Debug,
            LogLevel::Trace => log::LevelFilter::Trace,
        }
    }

    pub fn to_string(&self) -> String {
        format!("{:?}", self)
    }
}
