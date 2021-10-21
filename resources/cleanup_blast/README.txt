resources/cleanup_blast directory
==================================

Description:
------------
This directory contains the TVV sequence files and respective BLAST database files for use with the 'cleanup_blast.sh' script. 

The BLAST database will serve as the complete set of reference sequences that will be used to clean up the initial read binning procedure (carried out by BBsplit). This cleanup step will yield the final assignment of viral read counts by TVV species (e.g., Tvag isolate XYZ has 123 TVV1 reads, 456 TVV2 reads, 0 TVV3 reads...).


Contents:
------------
Two databases are currently available:
  - 'minimal_tvv_db' : one reference sequence per TVV species
  - 'comprehensive_tvv_db' : all coding-complete TVV genomes from NCBI GenBank plus additional  TVV sequences derived from Tvag transcriptome datasets, as per AR Manny et al, 2021.

Files: 
------------
  - .fasta files are the human-readable viral genome sequence files
  - .nhr, .nin, .nsq files are BLAST databases built from the fasta files and are just machine-readable

