use clap::{Parser, ValueEnum};

use crate::constants::{IDENTITY_THRESHOLDS, ID_THRESHOLD_ITERS, MINIMIZER_END_NTH_OVERLAP};

#[derive(Parser, Debug)]
#[command(
    name = "myloasm",
    about = "myloasm - high-resolution long-read metagenomic assembly with even low-fidelity reads.",
    version,
    author
)]

#[derive(Default, Clone)]
pub struct Cli {
    /// Input read file(s) -- multiple files are concatenated
    #[arg(num_args = 1.., required = true, value_name = "FASTQ/FASTA (.gz)")]
    pub input_files: Vec<String>,

    
    /// Number of threads to use for processing
    #[arg(short, long, default_value = "20", help_heading = "Basic Parameters")]
    pub threads: usize,
    
    /// Compression ratio (1/c k-mers selected). 
    #[arg(short, long, default_value = "11", help_heading = "Basic Parameters")]
    pub c: usize,

    /// Relaxed compression ratio during containment; must be > c
    #[arg(long, default_value_t=44, help_heading = "Misc")]
    pub contain_subsample_rate: usize,

    /// R9 (old nanopore) mode for low (~90%) accuracy reads. Experimental. 
    #[arg(long, help_heading = "Preset Parameters")]
    pub r941: bool,

    /// Output contigs with >= this number of reads
    #[arg(long, default_value_t = 1, help_heading = "Basic Parameters")]
    pub min_reads_contig: usize,

    /// Output directory for results
    #[arg(short, long, default_value = "myloasm-out", help_heading = "Basic Parameters")]
    pub output_dir: String,

    /// Disallow reads with < % identity for graph building (estimated from base qualities). 
    #[arg(long, default_value_t=90., help_heading = "Basic Parameters")]
    pub quality_value_cutoff: f64,

    /// Minimum overlap length for graph construction
    #[arg(long, default_value_t=500, help_heading = "Basic Parameters")]
    pub min_ol: usize,

    /// HiFi mode -- less aggressive chimeric read removal
    #[arg(long, help_heading = "Preset Parameters")]
    pub hifi: bool,
    
    /// Bloom filter size in GB
    #[arg(short, long, default_value_t=10., help_heading = "Misc")]
    pub bloom_filter_size: f64,

    /// Verbosity level. Written to the .log file in output directory.
    #[arg(short, long, value_enum, default_value = "debug", help_heading = "Misc")]
    pub log_level: LogLevel,

    /// No polishing (not recommended)
    #[arg(long, default_value_t=false, help_heading = "Misc")]
    pub no_polish: bool,

    /// Disable usage of SNPmers (not recommended)
    #[arg(long, default_value_t=false, help_heading = "Misc")]
    pub no_snpmers: bool,

    /// Base length of tip to remove 
    #[arg(long, default_value_t = 20000, help_heading = "Graph Parameters (advanced)")]
    pub tip_length_cutoff: usize,

    /// Number of reads in tips to remove
    #[arg(long, default_value_t = 3, help_heading = "Graph Parameters (advanced)")]
    pub tip_read_cutoff: usize,

    /// Snpmer identity threshold for containment and strict overlaps
    #[arg(long, default_value_t=IDENTITY_THRESHOLDS[ID_THRESHOLD_ITERS - 1] * 100., help_heading = "Overlap Parameters (advanced)")]
    pub snpmer_threshold_strict: f64,

    /// Snpmer identity threshold for relaxed overlaps
    #[arg(long, default_value_t=IDENTITY_THRESHOLDS[0] * 100., help_heading = "Overlap Parameters (advanced)")]
    pub snpmer_threshold_lax: f64,

    /// Binomial test error threhsold for relaxed overlaps
    #[arg(long, default_value_t=0.025, help_heading = "Overlap Parameters (advanced)")]
    pub snpmer_error_rate_lax: f64,

    /// Binomial test error threhsold for strict overlaps
    #[arg(long, default_value_t=0.00, help_heading = "Overlap Parameters (advanced)")]
    pub snpmer_error_rate_strict: f64,

    
    /// Cut overlaps with > (c * this) number of bases between minimizers
    #[arg(long, default_value_t=8., help_heading = "Overlap Parameters (advanced)")]
    pub absolute_minimizer_cut_ratio: f64,

    /// Cut overlaps with > (this) times more bases between minimizers than the best overlap
    #[arg(long, default_value_t=5., help_heading = "Overlap Parameters (advanced)")]
    pub relative_minimizer_cut_ratio: f64,

    /// Disables a SNPmer error overlap rescue heuristic during graph construction
    #[arg(long, help_heading = "Overlap Parameters (advanced)")]
    pub disable_error_overlap_rescue: bool,
    
    /// Base bubble popping threshold 
    #[arg(long, default_value_t=50000, help_heading = "Graph Parameters (advanced)")]
    pub small_bubble_threshold: usize,

    /// Cut z-edges that are < this times smaller than the adjacent overlaps
    #[arg(long, default_value_t=1.0, help_heading = "Graph Parameters (advanced)")]
    pub z_edge_threshold: f64,

    // ------ HIDDEN ARGUMENTS -----
    
    /// K-mer size (must be odd and < 24)
    #[arg(short, long, default_value = "21", help_heading = "Advanced", hide = true)]
    pub kmer_size: usize,

    /// Soft clips with < this # of bases are allowed for alignment
    #[arg(long, default_value_t=300, help_heading = "Alignment Parameters", hide = true)]
    pub maximal_end_fuzz: usize, 

    /// Maximum bubble length to pop; keep alternates
    #[arg(long, default_value_t=500000, help_heading = "Graph Parameters", hide = true)]
    pub max_bubble_threshold: usize,

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
