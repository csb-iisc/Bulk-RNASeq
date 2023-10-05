# Bulk-RNASeq
A repository for the analysis script for bulk RNA Sequencing data

## Requirements (with tested versions)
- curl == 7.81
- fastqc == 0.11.9
- fastp == 0.20.1
- STAR == 2.7.3a
- python == 3.10.6
    - multiqc == 1.14
    - htseq == 2.0.2

## How to use
1. Clone this repository/copy script from `~/Scripts/Bulk-RNASeq/`
2. Edit the number of threads and organism in the script (and the path if you aren't using the sequencing cluster)
3. If you have the fastq files, else follow step 4
    - Place them in the fastq folder
    - Ensure that the paired end naming convention similar to the following, you might have to edit the script otherwise: 
        - Sample1_1.fastq.gz
        - Sample1_2.fastq.gz
4. If you're downloading a public dataset: 
    - Go to [ENA](https://www.ebi.ac.uk/ena/browser/) and search for the BioProject accession number
    - Download report in tsv format and place it in the same folder as the script
5. Run the script
```bash 
./Bulk-RNASeq.sh
```
6. Wait
7. Profit
