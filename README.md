## ANALYSIS FOR HK68 H3N2 HA STEM DEEP MUTATIONAL SCANNING
This study aims to search for mutations that can escape HA stem-binding broadly neutralizing antibodies (bnAbs). The repository here describes the analysis for the deep mutational scanning experiment of HA2 residues 42, 45, 46, 47, 48, 49, 52, and 111 of A/Hong Kong/1/1968 (HK68).

### INPUT FILE
* All raw sequencing reads, which can be downloaded from NIH SRA database [PRJNA510654](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA510654), should be placed in fastq/ folder. The filename for read 1 should match those described in [./doc/SampleID.tsv](./doc/SampleID.tsv). The filename for read 2 should be the same as read 1 except "R1" is replaced by "R2".
* [Deep mutational scanning data from Bloom lab](https://github.com/jbloomlab/HA\_stalkbnAb\_MAP) are in the [./Bloom\_data/](./Bloom\_data/) folder
* [./doc/SampleID.tsv](./doc/SampleID.tsv): Describes the sample identity for each fastq file.
* [./doc/WTCodon.tsv](./doc/WTCodon.tsv): Describes the wild type nucleotide sequence for the residues of interest
* [./doc/ResiBarcode.tsv](./doc/ResiBarcode.tsv): Describes the nucleotide sequence for HA2 residues 43, 44, 50, and 112, which serves as an internal barcode. Silent mutations were introduced into HA2 residues 43, 44, 50, and 112, to index which residues were being randomized. 
* [./doc/HK68\_WTheatmap.tsv](./doc/HK68\_WTheatmap.tsv): Describes the location of the boxed (wild type residues) when plotting heatmap HK68.
* [./doc/WSN\_WTheatmap.tsv](./doc/WSN\_WTheatmap.tsv): Describes the location of the boxed (wild type residues) when plotting heatmap for WSN.
* [./doc/Perth09\_WTheatmap.tsv](./doc/Perth09\_WTheatmap.tsv): Describes the location of the boxed (wild type residues) when plotting heatmap for Perth09.

### ANALYSIS PIPELINE
1. [./script/HK68\_Stem\_read\_to\_count.py](./script/HK68_Stem_read_to_count.py): Converts raw reads to variant counts.
    - Input files: 
      - Raw sequencing reads in fastq/ folder
      - [./doc/SampleID.tsv](./doc/SampleID.tsv)
      - [./doc/WTCodon.tsv](./doc/WTCodon.tsv)
      - [./doc/ResiBarcode.tsv](./doc/ResiBarcode.tsv)
    - Output files: count/count\_\*.tsv
2. [./script/HK68\_Stem\_count\_to\_fitness.py](./script/HK68_Stem_count_to_fitness.py): Computes relative fitness based on read counts.
    - Input files: count/count\_\*.tsv
    - Output files:
      - [./result/Fitness\_S.tsv](./result/Fitness\_S.tsv)
      - [./result/Fitness\_D.tsv](./result/Fitness\_D.tsv)
3. [./script/HK68\_Stem\_compute\_resistance.py](./script/HK68\_Stem\_compute\_resistance.py): Computes relative resistance based on relative fitness.
    - Input files:
      - [./result/Fitness\_S.tsv](./result/Fitness\_S.tsv)
      - [./result/Fitness\_D.tsv](./result/Fitness\_D.tsv)
    - Output files:
      - [./result/Resist\_S.tsv](./result/Resist\_S.tsv)
      - [./result/Resist\_D.tsv](./result/Resist\_D.tsv)
4. [./script/HK68\_format\_heatmap.py](./script/HK68\_format\_heatmap.py): Generates the input file for double mutant heatmap plotting.
    - Input file: [./result/Resist\_D.tsv](./result/Resist\_D.tsv)
    - Output file: [./result/Heatmap\_D.tsv](./result/Heatmap\_D.tsv)
5. [./script/Perth09\_count2fit.py](./script/Perth09\_count2fit.py): Computes relative fitness from Perth09 count data.
    - Input file: Bloom\_data/Perth09\/\*.csv
    - Output file: [./result/Fitness\_Perth09.tsv](./result/Fitness\_Perth09.tsv)
6. [./script/WSN\_count2fit.py](./script/WSN\_count2fit.py): Computes relative fitness from WSN count data.
    - Input file: Bloom\_data/WSN/\*.csv
    - Output file: [./result/Fitness\_WSN.tsv](./result/Fitness\_WSN.tsv)
7. [./script/Compare\_strains.py](./script/Compare\_strains.py): Compile single mutant fitness data from HK68, Perth09, and WSN into a single file
    - Input files:
      - [./result/Resist\_S.tsv](./result/Resist\_S.tsv)
      - [./result/Fitness\_Perth09.tsv](./result/Fitness\_Perth09.tsv)
      - [./result/Fitness\_WSN.tsv](./result/Fitness\_WSN.tsv)
    - Output file: [./result/Fit\_compare.tsv](./result/Fit\_compare.tsv)

### PLOTTING
1. [./script/Plot\_rep\_cor.R](./script/Plot\_rep\_cor.R): Plots correlation between replicates 
    - Input files: 
      - [./result/Fitness\_S.tsv](./result/Fitness\_S.tsv)
      - [./result/Fitness\_D.tsv](./result/Fitness\_D.tsv)
    - Output file: [./graph/QC\_replicates\_cor.png](./graph/QC\_replicates\_cor.png)
2. [./script/Plot\_sample\_cor.R](./script/Plot\_sample\_cor.R): Compare relative resistance profiles across different antibody selections
   - Input files: 
      - [./result/Resist\_S.tsv](./result/Resist\_S.tsv)
      - [./result/Resist\_D.tsv](./result/Resist\_D.tsv)
   - Output file: 
      - [./graph/QC\_antibodies\_cor\_S.png](./graph/QC\_antibodies\_cor\_S.png)
      - [./graph/QC\_antibodies\_cor\_D.png](./graph/QC\_antibodies\_cor\_D.png)
3. [./script/Plot\_HK68\_single\_fit\_heatmap.R](./script/Plot\_HK68\_single\_fit\_heatmap.R): Plots the heatmap of single mutant relative fitness.
    - Input files:
      - [./doc/HK68\_WTheatmap.tsv](./doc/HK68\_WTheatmap.tsv)
      - [./result/Resist\_S.tsv](./result/Resist\_S.tsv)
    - Output file: [./graph/heatmap\_fit\_single.png](./graph/heatmap\_fit\_single.png)
4. [./script/Plot\_HK68\_resist\_heatmap.R](./script/Plot\_HK68\_resist\_heatmap.R): Plots the heatmap for single mutant relative resistance and also for the single mutant relative resistance with a background mutation.
    - Input files:
      - [./doc/HK68\_WTheatmap.tsv](./doc/HK68\_WTheatmap.tsv)
      - [./result/Resist\_S.tsv](./result/Resist\_S.tsv)
      - [./result/Heatmap\_D.tsv](./result/Heatmap\_D.tsv)
    - Output files: graph/heatmap\_\*esp\_single\*.png
5. [./script/Plot\_double\_heatmap.R](./script/Plot\_double\_heatmap.R): Plots the fitness heatmap and resistance heatmap for double mutant.
    - Input file: [./result/Heatmap\_D.tsv](./result/Heatmap\_D.tsv)
    - Output files: graph/heatmap\_\*\_double.png
6. [./script/Plot\_compare\_strains.R](./script/Plot\_compare\_strains.R): Generates scatter plots for comparing the relative fitness data of HK68, Perth09, and WSN.
    - Input file: [./result/Fit\_compare.tsv](./result/Fit\_compare.tsv)
    - Output files: graph/scatter\_fit\_\*.png
7. [./script/Plot\_resi45\_fit.R](./script/Plot\_resi45\_fit.R): Plots the fitness of different variants at HA2 residue 45 for HK68, Perth09, and WSN. 
    - Input file: [./result/Fit\_compare.tsv](./result/Fit\_compare.tsv)
    - Output file: [./graph/resi45\_fit.png](./graph/resi45\_fit.png)
8. [./script/Plot\_resist\_heatmap\_Perth09.R](./script/Plot\_resist\_heatmap\_Perth09.R):
    - Input files:
      - [./doc/Perth09\_WTheatmap.tsv](./doc/Perth09\_WTheatmap.tsv)
      - [./result/Fitness\_Perth09.tsv](./result/Fitness\_Perth09.tsv)
    - Output files: [./graph/heatmap\_FI6v3esp\_single\_Perth09.png](./graph/heatmap\_FI6v3esp\_single\_Perth09.png)
9. [./script/Plot\_resist\_heatmap\_WSN.R](./script/Plot\_resist\_heatmap\_WSN.R): 
    - Input files:
      - [./doc/WSN\_WTheatmap.tsv](./doc/WSN\_WTheatmap.tsv)
      - [./result/Fitness\_WSN.tsv](./result/Fitness\_WSN.tsv)
    - Output files:
      - [./graph/heatmap\_FI6v3esp\_single\_WSN.png](./graph/heatmap\_FI6v3esp\_single\_WSN.png)
      - [./graph/heatmap\_CR9114esp\_single\_WSN.png](./graph/heatmap\_CR9114esp\_single\_WSN.png)
