use clap::{Parser, ValueEnum};

#[derive(Parser, Debug)]
#[command(
    name = "phaller",
    about = "A tool for analyzing genomic sequences and building overlap graphs",
    version,
    author
)]
pub struct Cli {
    /// Input FASTQ files
    #[arg(required = true)]
    pub input_files: Vec<String>,

    /// K-mer size (must be odd)
    #[arg(short, long, default_value = "21")]
    pub kmer_size: usize,

    /// Compression ratio
    #[arg(short, long, default_value = "10")]
    pub c: usize,

    /// Number of threads to use for processing
    #[arg(short, long, default_value = "5")]
    pub threads: usize,

    /// Enable homopolymer compression
    #[arg(long)]
    pub homopolymer_compression: bool,

    /// Output directory for results
    #[arg(short, long, default_value = "assembly_out")]
    pub output_dir: String,

    /// Verbosity level
    #[arg(short, long, value_enum, default_value = "info")]
    pub log_level: LogLevel,

    /// Length of tip to remove
    #[arg(long, default_value_t = 20000)]
    pub tip_length_cutoff: usize,

    /// Number of reads in tips to remove
    #[arg(long, default_value_t = 4)]
    pub tip_read_cutoff: usize,
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
