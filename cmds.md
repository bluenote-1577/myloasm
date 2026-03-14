# Command-Line Help for `myloasm`

This document contains the help content for the `myloasm v0.5` command-line program.

**Command Overview:**

* [`myloasm`↴](#myloasm)

## `myloasm`

myloasm - high-resolution metagenomic assembly with noisy long reads. See online documentation for full options. 

EXAMPLE (Nanopore R10): myloasm nanopore_reads.fq.gz -o output_directory -t 50
EXAMPLE (PacBio HiFi): myloasm pacbio_reads.fq.gz -o output_directory -t 50 --hifi

**Usage:** `myloasm [OPTIONS] <FASTQ/FASTA (.gz)>...`

### **Arguments:**

* `<FASTQ/FASTA (.gz)>` — Input read file(s) -- multiple files are concatenated

### **Basic Algorithmic Parameters:**

* `-c`, `--c <C>` — Compression ratio (1/c k-mers selected). Must be <= 15
 [`11`] 

* `--kmc-db <KMC_DB>` — Use precomputed KMC database at this path for kmer counting. Must use -b and -k21 for KMC db creation with version v3
* `--quality-value-cutoff <QUALITY_VALUE_CUTOFF>` — Disallow reads with < % identity for graph building (estimated from base qualities)
 [`90`] 

* `--min-ol <MIN_OL>` — Minimum overlap length for graph construction
 [`500`] 

* `-b`, `--bloom-filter-size <BLOOM_FILTER_SIZE>` — Bloom filter size in GB. Increase for massive datasets if initial k-mer counting is a bottleneck (default: automatic estimation)
* `--aggressive-bloom` — More aggressive filtering of low-abundance k-mers. May be non-deterministic
* `--new-polish-trimming` — New mode: trim windows during polishing. Takes slightly longer, may incrementally improve polishing for some datasets
 [`false`] 

* `--parallel-graph-bridging` — Allow for parallel graph resolution of bridged repeats. This will make the assembly slightly worse, but may resolve a bottleneck for huge, complex metagenomes
* `--high-freq-kmer-threshold <HIGH_FREQ_KMER_THRESHOLD>` — Remove highest frequency k-mers (1 / this)
 [`100000`] 

* `--min-qual-polishing <MIN_QUAL_POLISHING>` — Disallow reads with < % identity for polishing (set to > 0 otherwise polishing may stall)
 [`75`] 

* `-k`, `--kmer-size <KMER_SIZE>` — K-mer size (must be odd and < 24)
 [`21`] 


### **Graph Parameters (advanced):**

* `--small-bubble-threshold <SMALL_BUBBLE_THRESHOLD>` — Base bubble popping length threshold; this gets multiplied by 5-30x during progressive graph cleaning
 [`50000`] 

* `--z-edge-threshold <Z_EDGE_THRESHOLD>` — Cut z-edges that are < this times smaller than the adjacent overlaps
 [`1`] 

* `--tip-length-cutoff <TIP_LENGTH_CUTOFF>` — Base length of tip to remove; this gets multiplied by 5-30x during simplification
 [`20000`] 

* `--tip-read-cutoff <TIP_READ_CUTOFF>` — Number of reads in tips to remove; this gets multiplied by 5-30x during simplification
 [`3`] 

* `--max-bubble-threshold <MAX_BUBBLE_THRESHOLD>` — Maximum bubble length to pop; keep alternates
 [`500000`] 


### **Miscellaneous Options:**

* `--no-polish` — No polishing (not recommended)
 [`false`] 

* `--no-snpmers` — Disable usage of SNPmers (not recommended)
 [`false`] 


### **Options:**

* `-o`, `--output-dir <OUTPUT_DIR>` — Output directory for results; created if it does not exist
 [`myloasm-out`] 

* `-t`, `--threads <THREADS>` — Number of threads to use for processing
 [`20`] 

* `--clean-dir` — Do not dump large intermediate data to disk (intermediate data is useful for rerunning)
* `-l`, `--log-level <LOG_LEVEL>` — Verbosity level. Warning: trace is very verbose
 [`debug`] 


  Possible values: `error`, `warn`, `info`, `debug`, `trace`

* `--markdown-help` — Print this markdown document

### **Output thresholds:**

* `--min-reads-contig <MIN_READS_CONTIG>` — Output contigs with >= this number of reads
 [`1`] 

* `--singleton-coverage-threshold <SINGLETON_COVERAGE_THRESHOLD>` — Remove singleton contigs with <= this estimated coverage depth (DP1 coverage; 99% identity coverage)
 [`3`] 

* `--secondary-coverage-threshold <SECONDARY_COVERAGE_THRESHOLD>` — Remove contigs with <= this estimated coverage depth and <= 2 reads (DP1 coverage; 99% identity coverage)
 [`1`] 

* `--absolute-coverage-threshold <ABSOLUTE_COVERAGE_THRESHOLD>` — Remove all contigs with <= this estimated coverage depth (DP1 coverage; 99% identity coverage)
* `--dereplication-ani <DEREPLICATION_ANI>` — Mark contigs with >= this average nucleotide identity (ANI) to a larger contig as alternate
 [`99`] 

* `--dereplication-length <DEREPLICATION_LENGTH>` — Mark contigs with > 90% aligned, < this length, and >= --dereplication-ani as alternate
 [`500000`] 


### **Overlap Parameters (advanced):**

* `--read-map-batch-size <READ_MAP_BATCH_SIZE>` — Batch size of indexing for read-to-read mapping and overlap stage. Higher = faster, but more memory
 [`1000000`] 

* `--snpmer-threshold-strict <SNPMER_THRESHOLD_STRICT>` — Snpmer identity threshold for containment and strict overlaps
 [`100`] 

* `--snpmer-threshold-lax <SNPMER_THRESHOLD_LAX>` — Snpmer identity threshold for relaxed overlaps
 [`99`] 

* `--snpmer-error-rate-lax <SNPMER_ERROR_RATE_LAX>` — Binomial test error parameter for relaxed overlaps
 [`0.025`] 

* `--snpmer-error-rate-strict <SNPMER_ERROR_RATE_STRICT>` — Binomial test error parameter strict overlaps
 [`0`] 

* `--contain-subsample-rate <CONTAIN_SUBSAMPLE_RATE>` — Relaxed compression ratio during containment; must be > c
 [`44`] 

* `--absolute-minimizer-cut-ratio <ABSOLUTE_MINIMIZER_CUT_RATIO>` — Cut overlaps with > (c * this) number of bases between minimizers on average
 [`8`] 

* `--relative-minimizer-cut-ratio <RELATIVE_MINIMIZER_CUT_RATIO>` — Cut overlaps with > (this) times more bases between minimizers than the best overlap on average
 [`5`] 

* `--disable-error-overlap-rescue` — Disables a SNPmer error overlap rescue heuristic during graph construction
* `--maximal-end-fuzz <MAXIMAL_END_FUZZ>` — Soft clips with < this # of bases are allowed for alignment
 [`300`] 


### **Technology Presets:**

* `--nano-r10` — (DEFAULT) R10 nanopore mode for sup/hac data (> ~97% median accuracy). Specifying this flag does not do anything for now
* `--nano-r9` — R9 (old nanopore) mode for low (~90%) accuracy reads. Experimental
* `--hifi` — PacBio HiFi mode -- assumes less chimericism and higher accuracy



<hr/>

<small><i>
    This document was generated automatically by
    <a href="https://crates.io/crates/clap-markdown"><code>clap-markdown</code></a>.
</i></small>

