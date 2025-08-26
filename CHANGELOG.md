v0.2.0 - More conservative and cleaner contig outputs with better polishing

- allow fasta reads (https://github.com/bluenote-1577/myloasm/issues/3)
- allow reads with non ACGT; send non-ACGT --> A character (https://github.com/bluenote-1577/myloasm/issues/4)
- fixed some bugs with polishing. Should be more robust and have less ~200 bp gaps. 
- remove ALL contigs with cov <= 1 (cli parameter now)
- remove singleton contigs with cov <= 3  (cli parameter now)
- fixed specific polishing issues for small circular genomes
- more aggressive dereplication of multiplied plasmids. 
- changed contig output format to output `duplicated-yes|no|possibly`. k-mer multiplicity is now in an extra field due to decimal being unruly especially during contig processing. 
