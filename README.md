## ANALYSIS FOR HK68 H3N2 HA STEM DEEP MUTATIONAL SCANNING
This study aims to search for mutations that can escape HA stem-binding broadly neutralizing antibodies (bnAbs). The repository here describes the analysis for the deep mutational scanning experiment of HA2 residues 42, 45, 46, 47, 48, 49, 52, and 111 of A/Hong Kong/1/1968 (HK68).

### INPUT FILE
* All raw sequencing reads, which can be downloaded from NIH SRA database [PRJNA510654](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA510654), should be placed in fastq/ folder. The filename for read 1 should match those described in [./doc/SampleID.tsv](./doc/SampleID.tsv). The filename for read 2 should be the same as read 1 except "R1" is replaced by "R2".
* [./doc/SampleID.tsv](./doc/SampleID.tsv): Describes the sample identity for each fastq file.
* [./doc/WTCodon.tsv](./doc/WTCodon.tsv): Describes the wild type nucleotide sequence for the residues of interest
* [./doc/ResiBarcode.tsv](./doc/ResiBarcode.tsv): Describes the nucleotide sequence for HA2 residues 43, 44, 50, and 112, which serves as an internal barcode. Silent mutations were introduced into HA2 residues 43, 44, 50, and 112, to index which residues were being randomized. 

### ANALYSIS PIPELINE
1. [./HK68\_Stem\_read\_to\_count.py](./HK68_Stem_read_to_count.py): Convert raw reads to variant counts
  - Input files: 
    - Raw sequencing reads in fastq/ folder
    - [./doc/SampleID.tsv](./doc/SampleID.tsv)
    - [./doc/WTCodon.tsv](./doc/WTCodon.tsv)
    - [./doc/ResiBarcode.tsv](./doc/ResiBarcode.tsv)
