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
