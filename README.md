# Bulk-RNASeq
A repository for the script to process raw fastq to counts/TPM for bulk RNA Sequencing data.

## Requirements (with tested versions)
- curl == 8.5.0
- fastqc == 0.12.1
- fastp == 0.23.4
- STAR == 2.7.11b
- samtools == 1.19.2
- python == 3.12.3
    - multiqc == 1.18
    - htseq == 2.0.5
    - pandas == 2.1.4
    - numpy == 1.26.4

## Input files
- If you have the fastq files, place them in the fastq folder.
- If you're downloading a public dataset: 
    - Go to [ENA](https://www.ebi.ac.uk/ena/browser/) and search for the BioProject accession number
    - Ensure that you have selected only the RNA seq files
    - Download report in tsv format
- For paired end sequencing, the script supports the following naming conventions. You might have to rename the files or edit the script otherwise.
    - SRA: `SampleName_1.fastq.gz` and `SampleName_2.fastq.gz`
    - Illumina: `SampleName_S1_L001_R1_001.fastq.gz` and `SampleName_S1_L001_R2_001.fastq.gz`

## Usage
```bash 
RNASeq.py [-h] [--thr THR] [--filereport FILEREPORT] [--num_downloads NUM_DOWNLOADS] [--preserve_fastq | --no-preserve_fastq] [--preserve_bam | --no-preserve_bam] [--fastq_dir FASTQ_DIR] [--fastqc_dir FASTQC_DIR] [--fastp_dir FASTP_DIR] [--fastqp_dir FASTQP_DIR] [--bam_dir BAM_DIR] [--counts_dir COUNTS_DIR] [--summary_dir SUMMARY_DIR] [--org ORG] [--genome GENOME] [--genome_dir GENOME_DIR] [--gtf_file GTF_FILE] [--genome_length GENOME_LENGTH]
```

**Examples:**
- Download fastq files from public data set and processes. Deletes the fastq and bam files. Organism: Human
```bash
./RNASeq.py --no-preserve_fastq --no-preserve_bam --filereport filereport_read_run_<BioProjectNumber>_tsv.txt --org hg38
```

- Process the fastq files in `./fastq` directory. Does not delete the fastq and bam files. Organism: Mouse
```bash
./RNASeq.py --preserve_fastq --preserve_bam --org mm10
```

**Important Options**
- `--filereport FILEREPORT`
    File path for the ENA file report in TSV format (if downloading).

- `--thr THR`
    Number of threads to use for processing. Default: `maximum cores - 2`. 

- `--org ORG`
    Organism to use. Available options: hg38 (human), mm10 (mouse). Default: `hg38`. 

- `--preserve_fastq, --no-preserve_fastq`
    Controls whether fastq files are preserved after trimming and alignment. Default: `no-preserve`. 

- `--preserve_bam, --no-preserve_bam`
    Controls whether bam files are preserved after counting. Default: `preserve`.

## Options

### Misc
- `-h, --help`
    Shows the help message with all available options. 

- `--thr THR`
    Number of threads to use for processing. Default: `maximum cores - 2`. 

### Download
- `--filereport FILEREPORT`
    File path for the ENA file report in TSV format (if downloading).

- `--num_downloads NUM_DOWNLOADS`
    Number of simultaneous fastq files to download. Increasing this might result in IP ban from ENA. Default: `10`. 

### Files

- `--preserve_fastq, --no-preserve_fastq`
    Controls whether fastq files are preserved after trimming and alignment. Default: `no-preserve`.  

- `--preserve_bam, --no-preserve_bam`
     Controls whether bam files are preserved after counting. Default: `preserve`.

### Working directory

- `--fastq_dir FASTQ_DIR`
    Path to the directory containing raw fastq files or directory to download fastq files to. Default: `./fastq`. 

- `--fastqc_dir FASTQC_DIR`
    Directory to store FastQC reports. Default: `./fastqc`. 

- `--fastp_dir FASTP_DIR`
    Directory for Fastp output files. Default: `./fastp`. 

- `--fastqp_dir FASTQP_DIR`
    Path to fastq files after trimming. Default: `./fastqp`. 

- `--bam_dir BAM_DIR`
    Directory for storing bam alignment files. Default: `./bam`. 

- `--counts_dir COUNTS_DIR`
    Output directory for raw counts and tpm. Default: `./counts`. 

- `--summary_dir SUMMARY_DIR`
    Directory to save summary files. Default: current `./summary`.

### Genome

- `--org ORG`
    Organism to use. Available options: hg38 (human), mm10 (mouse). Default: `hg38`. 

- `--genome GENOME`
    Path to the genome reference file. Default: `/mnt/nvme_sequencing/RNASeq`.

- `--genome_dir GENOME_DIR`
    Explicitly provide the directory containing genome files. Overrides automatic selection from genome and organism settings. 

- `--gtf_file GTF_FILE`
    Explicitly specify a GTF file for annotation. Overrides automatic selection from genome and organism settings.


## Generating Genome Index
Checkout the `GenomeGen` branch and follow instructions there.
