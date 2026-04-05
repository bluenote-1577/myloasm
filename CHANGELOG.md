# v0.5.1 (4-5-2026) - Extra checkpoints, faster graph cleaning + bugfix

Major changes

- Fixed a bug for crashing during the polishing stage (https://github.com/bluenote-1577/myloasm/issues/20)
- Made graph cleaning much faster for massive, complex datasets (https://github.com/bluenote-1577/myloasm/issues/18)
- Added more checkpoints

Results may change slightly. 

# v0.5.0 (3-12-2026) - RAM usage + read-mapping optimizations

Major changes

- Using huffman coding to reduce RAM further. Should optimize by ~10-20%. 
- Read-to-read mapping is accelerated by 3-4x in certain cases. This should reduce stalling behavior on certain datasets with high-depth, recombinant populations. Results will change slightly.  
- Fixed a major bug where myloasm would stall during the polishing step. This was caused by base quality runs of PHRED score 33. 

Minor changes

- Logging improved slightly
- (ALPHA) There's a way to optimize memory for the k-mer counting step slightly using kmc as a pre-assembly step. Will improve this; message me if your assembly is dying early. 

# v0.4.0 (1-28-2026) - At least 30% peak RAM than before at baseline, but possibly up to ~80% reduced (i.e., even 5-fold) for large metagenomes

Major changes

- Data structure changes (liberal usage of variable-length quantity encodings) **reduced memory by ~30-40% at baseline.**
- Reworking of read mapping logic: temporary mappings do not have to be stored in memory temporarily. For large metagenomes, **this could reduce peak RAM by multiple folds.** 

Minor changes

- Reworked the logic for removing contained reads slightly; it now proceeds in batches. This gave slightly better results on preliminary data. 

BREAKING:

- The checkpoints from previous runs for v0.3.0 and below will no longer work. That is, if you want to resume a disrupted v0.3.0 run using `myloasm exist ...`,  you will have to use v0.3.0. 

# v0.3.0 (12-19-2025) - Better polishing for multi-strain communities & runtime on complex metagenomes

Major changes

- Polishing may be much improved when multiple co-existing strains are assembled for nanopore data. Otherwise, results should not change much. 
- Bloom filter size is now automatically chosen by default (but can be overridden)
- Significantly improved runtime on complex metagenomes (graph cleaning + read mapping steps)
- Procedure for dereplicating spurious contigs prior to polishing changed slightly; should improve memory and time for complex metagenomes

Minor changes: 

- Added a new filter that removes poorly assembled, low-coverage contigs with mapping issues.
- Added more options for dereplication of alternate contigs

# v0.2.0 - More conservative and cleaner contig outputs with better polishing

- allow fasta reads (https://github.com/bluenote-1577/myloasm/issues/3)
- allow reads with non ACGT; send non-ACGT --> A character (https://github.com/bluenote-1577/myloasm/issues/4)
- fixed some bugs with polishing. Should be more robust and have less ~200 bp gaps. 
- remove ALL contigs with cov <= 1 (cli parameter now)
- remove singleton contigs with cov <= 3  (cli parameter now)
- fixed specific polishing issues for small circular genomes
- more aggressive dereplication of multiplied plasmids. 
- changed contig output format to output `duplicated-yes|no|possibly`. k-mer multiplicity is now in an extra field due to decimal being unruly especially during contig processing. 
