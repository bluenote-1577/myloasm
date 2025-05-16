use clap::{Parser, ValueEnum};
use clap_num::number_range;
use crate::constants::CLI_HEADINGS;

use crate::constants::{IDENTITY_THRESHOLDS, ID_THRESHOLD_ITERS};

#[derive(Parser, Debug)]
#[command(
    name = "myloasm",
    about = 
    "myloasm - high-resolution metagenomic assembly with noisy long reads.\n\nEXAMPLE (Nanopore R10): myloasm nanopore_reads.fq.gz -o output_directory -t 50\nEXAMPLE (PacBio HiFi): myloasm pacbio_reads.fq.gz -o output_directory -t 50 --hifi",
    version,
    author
)]

#[derive(Default, Clone)]
pub struct Cli {
    /// Input read file(s) -- multiple files are concatenated
    #[arg(num_args = 1.., required = true, value_name = "FASTQ/FASTA (.gz)")]
    pub input_files: Vec<String>,

    /// (DEFAULT) R10 nanopore mode for sup/hac data (> ~97% median accuracy). Specifying this flag does not do anything for now.
    #[arg(long, help_heading = CLI_HEADINGS[0])]
    pub nano_r10: bool,

    /// R9 (old nanopore) mode for low (~90%) accuracy reads. Experimental. 
    #[arg(long, help_heading = CLI_HEADINGS[0], hide = true)]
    pub nano_r9: bool,

    /// PacBio HiFi mode -- assumes less chimericism and higher accuracy
    #[arg(long, help_heading = CLI_HEADINGS[0])]
    pub hifi: bool,
        
    /// Output directory for results; created if it does not exist
    #[arg(short, long, default_value = "myloasm-out")]
    pub output_dir: String,

    /// Number of threads to use for processing
    #[arg(short, long, default_value = "20")]
    pub threads: usize,

     /// Do not dump large intermediate data to disk (intermediate data is useful for rerunning)
    #[arg(long)]
    pub clean_dir: bool,
   
    /// Compression ratio (1/c k-mers selected). Must be <= 15  
    #[arg(short, long, default_value = "11", help_heading = CLI_HEADINGS[1])]
    pub c: usize,

    /// Output contigs with >= this number of reads
    #[arg(long, default_value_t = 1, help_heading = CLI_HEADINGS[1])]
    pub min_reads_contig: usize,
    
    /// Disallow reads with < % identity for graph building (estimated from base qualities) 
    #[arg(long, default_value_t=90., help_heading = CLI_HEADINGS[1])]
    pub quality_value_cutoff: f64,

    /// Minimum overlap length for graph construction
    #[arg(long, default_value_t=500, help_heading = CLI_HEADINGS[1])]
    pub min_ol: usize,
        
    /// Bloom filter size in GB. Increase for massive datasets
    #[arg(short, long, default_value_t=10., help_heading = CLI_HEADINGS[1])]
    pub bloom_filter_size: f64,

    /// More aggressive filtering of low-abundance k-mers. May be non-deterministic
    #[arg(long, help_heading = CLI_HEADINGS[1])]
    pub aggressive_bloom: bool,

    
    /// Verbosity level. Warning: trace is very verbose
    #[arg(short, long, value_enum, default_value = "debug")]
    pub log_level: LogLevel,

    /// No polishing (not recommended)
    #[arg(long, default_value_t=false, help_heading = CLI_HEADINGS[2], hide = true)]
    pub no_polish: bool,

    /// Disable usage of SNPmers (not recommended)
    #[arg(long, default_value_t=false, help_heading = CLI_HEADINGS[2], hide = true)]
    pub no_snpmers: bool,
    
    /// Snpmer identity threshold for containment and strict overlaps
    #[arg(long, default_value_t=IDENTITY_THRESHOLDS[ID_THRESHOLD_ITERS - 1] * 100., help_heading =CLI_HEADINGS[3])]
    pub snpmer_threshold_strict: f64,

    /// Snpmer identity threshold for relaxed overlaps
    #[arg(long, default_value_t=IDENTITY_THRESHOLDS[0] * 100., help_heading =CLI_HEADINGS[3])]
    pub snpmer_threshold_lax: f64,

    /// Binomial test error parameter for relaxed overlaps
    #[arg(long, default_value_t=0.025, help_heading =CLI_HEADINGS[3])]
    pub snpmer_error_rate_lax: f64,

    /// Binomial test error parameter strict overlaps
    #[arg(long, default_value_t=0.00, help_heading =CLI_HEADINGS[3])]
    pub snpmer_error_rate_strict: f64,

    /// Relaxed compression ratio during containment; must be > c
    #[arg(long, default_value_t=44, help_heading = CLI_HEADINGS[3])]
    pub contain_subsample_rate: usize,
    
    /// Cut overlaps with > (c * this) number of bases between minimizers on average
    #[arg(long, default_value_t=8., help_heading =CLI_HEADINGS[3])]
    pub absolute_minimizer_cut_ratio: f64,

    /// Cut overlaps with > (this) times more bases between minimizers than the best overlap on average
    #[arg(long, default_value_t=5., help_heading =CLI_HEADINGS[3])]
    pub relative_minimizer_cut_ratio: f64,

    /// Disables a SNPmer error overlap rescue heuristic during graph construction
    #[arg(long, help_heading =CLI_HEADINGS[3])]
    pub disable_error_overlap_rescue: bool,
    
    /// Base bubble popping length threshold; this gets multiplied by 5-30x during progressive graph cleaning
    #[arg(long, default_value_t=50000, help_heading = CLI_HEADINGS[4])]
    pub small_bubble_threshold: usize,

    /// Cut z-edges that are < this times smaller than the adjacent overlaps
    #[arg(long, default_value_t=1.0, help_heading = CLI_HEADINGS[4])]
    pub z_edge_threshold: f64,

    /// Base length of tip to remove; this gets multiplied by 5-30x during simplification
    #[arg(long, default_value_t = 20000, help_heading = CLI_HEADINGS[4])]
    pub tip_length_cutoff: usize,

    /// Number of reads in tips to remove; this gets multiplied by 5-30x during simplification
    #[arg(long, default_value_t = 3, help_heading = CLI_HEADINGS[4])]
    pub tip_read_cutoff: usize,


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
