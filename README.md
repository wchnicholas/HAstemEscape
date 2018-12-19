## ANALYSIS FOR HK68 H3N2 HA STEM DEEP MUTATIONAL SCANNING
This study aims to search for mutations that can escape HA stem-binding broadly neutralizing antibodies (bnAbs). The repository here describes the analysis for the deep mutational scanning experiment of HA2 residues 42, 45, 46, 47, 48, 49, 52, and 111 of A/Hong Kong/1/1968 (HK68).

### INPUT FILE
* All sequencing raw reads, which can be downloaded from NIH SRA database [PRJNA510654](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA510654), should be placed in fastq/ folder. The filename for read 1 should match those described in [./doc/SampleID.tsv](./doc/SampleID.tsv). The filename for read 2 should be the same as read 1 except "R1" is replaced by "R2".
* [./doc/SampleID.tsv](./doc/SampleID.tsv): Describes the sample identity for each fastq file.

### ANALYSIS PIPELINE
