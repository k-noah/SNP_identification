# SNP_identification
Script that identifies SNPs &amp; indels from aligned fasta files.

## Getting Started

Environment requires Python3, with packages numpy, pandas and biopython.
Executable script may work on Windows, without python installation.

## Requirements

- files are in .fasta_aln format
- only 2 strains per file
- reference strain is the first sequence
- fasta description starts with the strain information, ends with the gene
  e.g. 0C2E_isolate_RP1_3623|ASV59_RS01495|nan|protein tonB|WP_023086850.1
       will output with 0C2E_isolate_RP1_3623 as the strain, WP_023086850.1 as the gene
