#!/usr/bin/env python
import os
import sys
import time
import glob
import subprocess
import argparse
import pandas as pd
import numpy as np

def main():
    parser = argparse.ArgumentParser(prog='RNASeq.py', description='RNA-Seq pipeline')
    parser.add_argument('--thr', type=int, default=None, help='Number of threads to use')
    parser.add_argument('--filereport', type=str, help='File report tsv file from ENA')
    parser.add_argument('--num_downloads', type=int, default=10, help='Number of simultaneous fastq files to download')
    parser.add_argument('--preserve_fastq', action=argparse.BooleanOptionalAction, default=False, help='Preserve fastq files after alignment')
    parser.add_argument('--preserve_bam', action=argparse.BooleanOptionalAction, default=True, help='Preserve bam files after counting')
    parser.add_argument('--fastq_dir', type=str, default='./fastq', help='Path to the fastq files')
    parser.add_argument('--fastqc_dir', type=str, default='./fastqc', help='Path to the fastqc files')
    parser.add_argument('--fastp_dir', type=str, default='./fastp', help='Path to the fastp files')
    parser.add_argument('--fastqp_dir', type=str, default='./fastqp', help='Path to the fastq files after trimming')
    parser.add_argument('--bam_dir', type=str, default='./bam', help='Path to the bam files')
    parser.add_argument('--counts_dir', type=str, default='./counts', help='Path to the counts files')
    parser.add_argument('--summary_dir', type=str, default='./summary', help='Path to the summary files')
    parser.add_argument('--org', type=str, default='hg38', help='Organism to use. Available options: hg38, mm10')
    parser.add_argument('--genome', type=str, default='/mnt/nvme_sequencing/RNASeq', help='Path to the genome')
    parser.add_argument('--genome_dir', type=str, default=None, help='Explicitly provide genome dir. Overrides autoselected from genome/org')
    parser.add_argument('--gtf_file',type=str, default=None, help='Explicitly provide gtf file. Overrides autoselected from genome/org')
    parser.add_argument('--genome_length',type=str, default=None, help='Explicitly provide gene length file. Overrides autoselected from genome/org')
    args = parser.parse_args()
    thr = args.thr if args.thr else (os.cpu_count()-2)
    check_programs()
    mkdirs(args.fastq_dir, args.fastqc_dir, args.fastp_dir, args.fastqp_dir, args.bam_dir, args.summary_dir)
    if args.filereport:
        download_fastq(args.filereport, args.num_downloads, args.fastq_dir)
    fastqc(args.fastq_dir, args.fastqc_dir, thr)
    df = file_pair(args.fastq_dir)
    trimming(df, args.fastq_dir, args.fastp_dir, args.fastqp_dir, thr, args.preserve_fastq)
    fastq_dir = args.fastqp_dir
    fastqc(fastq_dir, args.fastp_dir, thr)
    align(df, fastq_dir, args.bam_dir, thr, args.genome_dir, args.org, args.genome, args.preserve_fastq)
    sort_bam(args.bam_dir, thr)
    bam_to_counts(args.bam_dir, args.counts_dir, thr, args.gtf_file, args.org, args.genome, args.preserve_bam)
    counts_to_tpm(args.counts_dir, args.genome_length, args.org, args.genome)
    multiqc(args.summary_dir, args.fastqc_dir, args.fastp_dir, args.bam_dir, args.counts_dir)

def check_programs():
    # Check if the programs are installed
    programs = ['curl', 'fastqc','multiqc', 'fastp', 'STAR', 'samtools', 'htseq-count']
    for program in programs:
        try:
            subprocess.run([program, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except FileNotFoundError:
            print(f'{program} is not installed.')
            sys.exit(1) 

def mkdirs(*path_list):
    # Make directories if they don't exist
    for dir_path in path_list:
        os.makedirs(dir_path, exist_ok=True)

curl_commands = [
    'curl',
    '--retry', '999',
    '--retry-max-time', '0',
    '--retry-all-errors',
    '-C', '-', # Resume download
    '-Z', # Parallel mode
    '--remote-name-all', # Save files with the original name
    ]

def download_fastq(filereport, ndl, fastq_dir):
    print(f'Downloading fastq files in {filereport}')
    df = pd.read_csv(filereport, sep='\t')
    if 'fastq_ftp' not in df.columns:
        print('Error: fastq_ftp column not found in the file report')
        sys.exit(1)
    fastq_files = df['fastq_ftp'].str.split(';', expand=True).stack().tolist()
    new_commands = curl_commands + ['--parallel-max', str(ndl)] + fastq_files
    subprocess.run(new_commands, cwd=fastq_dir)

def fastqc(fastq_dir, fastqc_dir, thr):
    print(f'Running Quality Check with FastQC on {fastq_dir}')
    print(time.strftime('%A %d %B %Y %I:%M:%S %p %Z'))
    fastq_files = sorted(glob.glob(f"{fastq_dir}/*.fastq.gz"))
    fastqc_command = [
        'fastqc', 
        '-t', str(thr), 
        '-o', fastqc_dir, 
        ] + fastq_files
    subprocess.run(fastqc_command)

def file_pair(fastq_dir):
    fastq_files = sorted(glob.glob('*.fastq.gz', root_dir=fastq_dir))
    df = pd.DataFrame(fastq_files, columns=['fname'])
    # Extract sample name and read number using regex
    df[['sname','read']] = df['fname'].str.extract(r'^(.*?)(?:(?:_S\d+_L\d+_R|_)(\d+)(?:_\d+)?)?\.fastq\.gz$')
    # Fill missing sample names with full filename (single-end case)
    df['sname'] = df['sname'].fillna(df['fname'].str.replace('.fastq.gz$', '', regex=True))
    df['read'] = df['read'].fillna('1')
    # Pivot table to get Read 1 and Read 2 columns
    df = df.pivot(index='sname', columns='read', values='fname')
    # Set as paired if both columns exist
    df['paired'] = df.notna().all(axis=1)
    # If file in 2 but not 1 move to 1
    df['1'] = df['1'].fillna(df['2'])
    # Drop the 2 column if not paired
    df.loc[~df['paired'], '2'] = None
    return df

def trimming(df, fastq_dir, fastp_dir, fastqp_dir, thr, preserve_fastq):
    print(f"Trimming with fastp started on {fastq_dir}")
    print(time.strftime('%A %d %B %Y %I:%M:%S %p %Z'))
    fastp_command = [
        'fastp',
        '-w', str(thr),
        '-r',
        '-l', '36',
        ]
    trimming_single(df[~df['paired']], fastp_command, fastq_dir, fastp_dir, fastqp_dir, preserve_fastq)
    trimming_paired(df[df['paired']], fastp_command, fastq_dir, fastp_dir, fastqp_dir, preserve_fastq)

def trimming_single(df, fastp_command, fastq_dir, fastp_dir, fastqp_dir, preserve_fastq):
    for sample_name, fastq_files in df.iterrows():
        new_commands = fastp_command + [
            '-i', f'{fastq_dir}/{fastq_files["1"]}',
            '-o', f'{fastqp_dir}/{fastq_files["1"]}',
            '-j', f'{fastp_dir}/{sample_name}_fastp.json',
            '-h', f'{fastp_dir}/{sample_name}_fastp.html',
            ]
        subprocess.run(new_commands)
        if not preserve_fastq:
            os.remove(f'{fastq_dir}/{fastq_files["1"]}')

def trimming_paired(df, fastp_command, fastq_dir, fastp_dir, fastqp_dir, preserve_fastq):
    for sample_name, fastq_files in df.iterrows():
        new_commands = fastp_command + [
            '-c',
            '--detect_adapter_for_pe',
            '-i', f'{fastq_dir}/{fastq_files["1"]}',
            '-I', f'{fastq_dir}/{fastq_files["2"]}',
            '-o', f'{fastqp_dir}/{fastq_files["1"]}',
            '-O', f'{fastqp_dir}/{fastq_files["2"]}',
            '-j', f'{fastp_dir}/{sample_name}_fastp.json',
            '-h', f'{fastp_dir}/{sample_name}_fastp.html',
            ]
        subprocess.run(new_commands)
        if not preserve_fastq:
            os.remove(f'{fastq_dir}/{fastq_files["1"]}')
            os.remove(f'{fastq_dir}/{fastq_files["2"]}')

def load_genome(genome_dir, bam_dir, remove=False):
    load = 'Remove' if remove else 'LoadAndExit'
    subprocess.run([
        'STAR',
        '--genomeLoad', load,
        '--genomeDir', genome_dir,
        '--outFileNamePrefix', bam_dir,
        ])

def align(df, fastq_dir, bam_dir, thr, genome_dir, org, genome, preserve_fastq):
    print(f"Mapping with STAR started on {fastq_dir}")
    print(time.strftime('%A %d %B %Y %I:%M:%S %p %Z'))
    genome_dir = genome_dir if genome_dir else genome+'/'+org+'_STARIndexed/'
    # Check if genome_dir is present
    if not os.path.exists(genome_dir):
        print(f'Error: Genome directory {genome_dir} does not exist.')
        sys.exit(1)
    load_genome(genome_dir, f'{bam_dir}/', remove=False)
    STAR_commands = [
        'STAR',
        '--runThreadN', str(thr),
        '--genomeDir', genome_dir,
        '--genomeLoad', 'LoadAndKeep',
        '--readFilesCommand', 'zcat',
        '--outSAMattributes', 'All',
        '--outSAMtype', 'BAM', 'Unsorted',
        '--quantMode', 'GeneCounts',
    ]
    align_single(df[~df['paired']], STAR_commands, bam_dir, fastq_dir, preserve_fastq)
    align_paired(df[df['paired']], STAR_commands, bam_dir, fastq_dir, preserve_fastq)
    load_genome(genome_dir, f'{bam_dir}/', remove=True)

def align_single(df, STAR_commands, bam_dir, fastq_dir, preserve_fastq):
    for sample_name, fastq_files in df.iterrows():
        new_commands = STAR_commands + [
            '--outFileNamePrefix', f'{bam_dir}/{sample_name}_',
            '--readFilesIn', f'{fastq_dir}/{fastq_files["1"]}'
            ]
        subprocess.run(new_commands)
        if not preserve_fastq:
            os.remove(f'{fastq_dir}/{fastq_files["1"]}')

def align_paired(df, STAR_commands, bam_dir, fastq_dir, preserve_fastq):
    for sample_name, fastq_files in df.iterrows():
        new_commands = STAR_commands + [
            '--outFileNamePrefix', f'{bam_dir}/{sample_name}_',
            '--readFilesIn', f'{fastq_dir}/{fastq_files["1"]}', f'{fastq_dir}/{fastq_files["2"]}'
            ]
        subprocess.run(new_commands)
        if not preserve_fastq:
            os.remove(f'{fastq_dir}/{fastq_files["1"]}')
            os.remove(f'{fastq_dir}/{fastq_files["2"]}')

def sort_bam(bam_dir, thr):
    print(f'Sorting with samtools started on {bam_dir}')
    print(time.strftime('%A %d %B %Y %I:%M:%S %p %Z'))
    bam_files = sorted(glob.glob(f'{bam_dir}/*_Aligned.out.bam'))
    sam_commands = [
        'samtools', 'sort',
        '-@', str(thr),
    ]
    for bam_file in bam_files:
        out_bam = bam_file.replace('_Aligned.out.bam', '_sorted.bam')
        new_commands = sam_commands + ['-n', f'{bam_file}', '-o', f'{out_bam}']
        subprocess.run(new_commands)
        os.remove(bam_file)

def bam_to_counts(bam_dir, counts_dir, thr, gtf_file, org, genome, preserve_bam):
    print(f'Getting Raw Counts with HTSeq from {bam_dir}')
    print(time.strftime('%A %d %B %Y %I:%M:%S %p %Z'))
    gtf_file = gtf_file if gtf_file else genome+'/'+org+'_STARIndexed/annotation.gtf'
    bam_files = sorted(glob.glob(f'{bam_dir}/*_sorted.bam'))
    stranded = strandedness(bam_files[0])
    print('Stranded-ness:', stranded)
    htseq_command = [
        'htseq-count',
        '-n', str(thr),
        '-s', stranded,
        '--with-header',
        '--additional-attr=gene_name',
        '-c', f'{counts_dir}/RawCounts_htseq.tsv',
        ] + bam_files + [gtf_file]
    subprocess.run(htseq_command)
    if not preserve_bam:
        os.remove(bam_files)

def strandedness(bam_file):
    reads_file = bam_file.replace('_sorted.bam', '_ReadsPerGene.out.tab')
    df = pd.read_csv(reads_file, sep='\t', header=None, index_col=0)
    strdict = {1: 'no', 2: 'yes', 3: 'reverse'}
    stridx = df.iloc[4:,:].sum(axis=0).idxmax()
    return strdict[stridx]

def multiqc(summary_dir, *dirs):
    print('Creating summary report with MultiQC')
    print(time.strftime('%A %d %B %Y %I:%M:%S %p %Z'))
    multiqc_command = [
        'multiqc',
        '-o', summary_dir,
        ] + list(dirs)
    subprocess.run(multiqc_command)

def counts_to_tpm(count_dir, genome_length, org, genome):
    genome_length = genome_length if genome_length else genome+'/'+org+'_STARIndexed/GeneLength_kb.tsv'
    df = pd.read_csv(f'{count_dir}/RawCounts_htseq.tsv', sep='\t', index_col=[0,1])
    df.drop(df.tail(5).index, inplace=True)
    gAnnot = pd.read_csv(genome_length, delimiter='\t', names=['GeneID','Len'], index_col=0, header=None)
    # Normalize to TPM
    if (len(df)!=len(gAnnot)):
        raise('Gene IDs not proper')
    df = df.divide(gAnnot['Len'], axis='rows', level=0)
    df = df.divide(df.sum(), axis='columns') * 1e6
    df = np.log2(df+1)
    # Save TPM counts
    df.to_csv(f'{count_dir}/TPMCounts.tsv',sep='\t')
    
if __name__ == '__main__':
    main()
