# Bulk-RNASeq
A repository for the script to process raw fastq to counts/TPM for bulk RNA Sequencing data.

This branch contains script to generate genome index. The script has been adopted from the [cellranger reference build pipeline](https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps).

## Requirements (with tested versions)
- curl == 8.5.0
- fastqc == 0.12.1
- fastp == 0.23.4
- STAR == 2.7.11b
- python == 3.12.3
    - multiqc == 1.18
    - htseq == 2.0.5

Install with
```bash
sudo apt install curl fastqc faspt rna-star python3-htseq multiqc

```

**Warning:** The STAR on Ubuntu 24.04 repo has a bug that doesn't parallelize for genomeGenerate, other functions do not seem to be affected. Download and use the [GitHub release](https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip) for this part.


## Genome versions on Sequencing Cluster
- hg38: GENCODE v44/Ensembl 110
- mm10: GENCODE v33/Ensembl 110

## Generating Genome Index
1. Clone this repository and checkout the `GenomeGen` branch.
2. Edit the input parameters:
    - `genome_path`: Path to save the genome and index files in. 
    - `thr`: number of threads to parallelize over.
    - `org`: Organism name (hg38 or mm10, other organism might need further modifications of script).
    - `genome_release`: Genome version as per Ensembl.
    - `gencode_release`: Annotation version as per GENCODE. Refer to [release history](https://www.gencodegenes.org/human/releases.html).

3. Run the script
```bash 
./STAR_genomeGenerate.sh
```
4. The script downlods the Genome and saves it in the `genome_path` as `<org>-<genome_release>` if it is not present already.
5. The script processes the genome and stores the index files in `<org>_STARIndexed` for further use.
