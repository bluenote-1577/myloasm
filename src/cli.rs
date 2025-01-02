use clap::{Parser, ValueEnum};

#[derive(Parser, Debug)]
#[command(
    name = "myloasm",
    about = "myloasm - strain-resolved metagenomic long-read assembly with noisy sequences.",
    version,
    author
)]
pub struct Cli {
    /// Input FASTQ files
    #[arg(num_args = 1..)]
    pub input_files: Vec<String>,

    /// K-mer size (must be odd)
    #[arg(short, long, default_value = "21")]
    pub kmer_size: usize,

    /// Compression ratio
    #[arg(short, long, default_value = "10")]
    pub c: usize,

    /// Number of threads to use for processing
    #[arg(short, long, default_value = "10")]
    pub threads: usize,

    /// Enable homopolymer compression
    #[arg(long, hide=true)]
    pub homopolymer_compression: bool,

    /// Output directory for results
    #[arg(short, long, default_value = "output")]
    pub output_dir: String,

    /// Verbosity level
    #[arg(short, long, value_enum, default_value = "info")]
    pub log_level: LogLevel,

    /// Length of tip to remove
    #[arg(long, default_value_t = 20000, help_heading = "Graph Parameters")]
    pub tip_length_cutoff: usize,

    /// Number of reads in tips to remove
    #[arg(long, default_value_t = 3, help_heading = "Graph Parameters")]
    pub tip_read_cutoff: usize,

    /// Do not use snpmers; standard overlap assembly
    #[arg(long, default_value_t=false, help_heading = "Overlap Parameters")]
    pub no_snpmers: bool,

    /// Disallow reads with < this accuracy from Q-scores for the overlap step. 
    #[arg(long, default_value_t=90.)]
    pub quality_value_cutoff: f64,

    /// Snpmer identity threshold for overlaps
    #[arg(long, default_value_t=99.9, help_heading = "Overlap Parameters")]
    pub snpmer_threshold: f64,

    /// Error rate for snpmers for binomial test
    #[arg(long, default_value_t=0.025, help_heading = "Overlap Parameters")]
    pub snpmer_error_rate: f64,

    #[arg(long, default_value_t=20, help_heading = "Overlap Parameters")]
    pub contain_subsample_rate: usize,

    /// Maximum bubble length to pop; keep alternates
    #[arg(long, default_value_t=300000, help_heading = "Graph Parameters")]
    pub max_bubble_threshold: usize,

    /// Small bubble length to pop; discard alternates
    #[arg(long, default_value_t=50000, help_heading = "Graph Parameters")]
    pub small_bubble_threshold: usize,


    /// Bloom filter size in GB
    #[arg(short, long, default_value_t=3.)]
    pub bloom_filter_size: f64,

    /// Minimum number of reads in output contigs
    #[arg(long, default_value_t = 2)]
    pub min_reads_contig: usize,

    /// Don't map to reads minimap2 (TODO)
    #[arg(long, default_value_t=false)]
    pub no_minimap2: bool,

    /// HiFi mode (--snpmer-threshold 100 --snpmer-error-rate 0.001)
    #[arg(long, help_heading = "Preset Parameters", hide=true)]
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
}
