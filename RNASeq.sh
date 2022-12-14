#!/bin/bash

############################################################################################
############################################################################################
##   A complete automated RNA-Seq pipeline                                                ##
##   Author : Susmita Mandal                                                              ##
##   Date: 13/10/2020                                                                     ##
##   Make three folders fastq, BAM and htseq                                              ##
##   Download the fastq files in the fastq folder and run this script                     ##
##   When prompted Organism: type hg38 or mm10 and ENTER                                  ##
##   When prompted Cluster threads: type the no of cores available to you. eg. 10,20,40   ##
############################################################################################
############################################################################################

# Input parameters
path='/mnt/csb-seq/'
STAR='_STARIndexed'
org='mm10'
thr=14
ndl=10

# exit when any command fails
set -e

# Setting file descriptor size to 999999 & file size to unlimited
ulimit -n 999999
ulimit -f unlimited

# Make directories if they don't exist
for dir in {fastq,fastqc,fastqp,fastp,BAM,htseq}; do
	[ -d $dir ] || mkdir $dir
done

cd fastq

# Download fastq files if given the tsv file
if [ -f ../filereport_read_run_*_tsv.txt ]; then
	echo "TSV file found, attempting to download"
	awk -F '\t' '
	NR==1 {
    		for (i=1; i<=NF; i++) {
        		f[$i] = i
    		}
	}\
	{ print $(f["fastq_ftp"])}
	' ../filereport_read_run_*_tsv.txt | tail -n +2 -q | awk -F ';' '
		{for(i=1;i<=NF;i++) print $i}
	' > ../download_links.txt
	curl --retry 999 --retry-max-time 0 --retry-all-errors -C - -Z --parallel-max $ndl --remote-name-all $(cat ../download_links.txt)
fi


# Check Quality of fastq using fastqc
echo "Quality Check b4 Pre-Processing"
date

fastqc \
	-t $thr\
	-o ../fastqc\
	*.fastq.gz

multiqc	-o ../fastqc/ ../fastqc/

echo "Pre-processing"
date

if [ `ls *_1.fastq.gz 2>/dev/null | wc -l` -gt 0 ] ; then
    for i in $(ls *.fastq.gz | sed -e 's/_1.fastq.gz//' -e 's/_2.fastq.gz//' | sort -u); do
		fastp -w $thr -r -l 36 -c --detect_adapter_for_pe -i ${i}_1.fastq.gz -I ${i}_2.fastq.gz -o ../fastqp/${i}_1.fastq.gz -O ../fastqp/${i}_2.fastq.gz -j ../fastp/${i}_fastp.json -h ../fastp/${i}_fastp.html
		
		rm ${i}_1.fastq.gz ${i}_2.fastq.gz
	done
else
	for i in $(ls *.fastq.gz | sed -e 's/.fastq.gz//' | sort -u);do
		fastp -w $thr -r - l 36 -i ${i}.fastq.gz -o ../fastqp/${i}.fastq.gz -j ../fastp/${i}_fastp.json -h ../fastp/${i}_fastp.html

		rm ${i}.fastq.gz
	done
fi

cd ../fastqp

echo "Quality Check after Pre-Processing"
date

fastqc \
	-t $thr\
	-o ../fastp\
	*.fastq.gz


multiqc	-o ../fastp/ ../fastp/

# Mapping fastq using STAR-aligner
echo "Mapping Started"
date

STAR --genomeLoad LoadAndExit --genomeDir $path$org$STAR

## For paired-end reads
if [ `ls *_1.fastq.gz 2>/dev/null | wc -l` -gt 0 ] ; then

        for i in $(ls *.fastq.gz | sed -e 's/_1.fastq.gz//' -e 's/_2.fastq.gz//' | sort -u); do

        STAR \
	--runThreadN $thr\
	--genomeDir $path$org$STAR\
	--genomeLoad LoadAndKeep\
	--readFilesCommand zcat\
	--readFilesIn ${i}_1.fastq.gz ${i}_2.fastq.gz\
	--outSAMattributes All\
	--outSAMtype BAM Unsorted\
	--quantMode GeneCounts\
	--outFileNamePrefix ../BAM/${i}_
	
	rm ${i}_1.fastq.gz ${i}_2.fastq.gz
        
done

else
## For normal reads

       for i in $(ls *.fastq.gz | sed -e 's/.fastq.gz//' | sort -u);do

	STAR \
	--runThreadN $thr\
	--genomeDir $path$org$STAR\
	--genomeLoad LoadAndKeep\
	--readFilesCommand zcat\
	--readFilesIn ${i}.fastq.gz\
	--outSAMattributes All\
	--outSAMtype BAM Unsorted\
	--quantMode GeneCounts\
	--outFileNamePrefix ../BAM/${i}_

	rm ${i}.fastq.gz

done

fi

STAR --genomeLoad Remove --genomeDir $path$org$STAR

date
echo "Mapping Done"

cd ../BAM/

##Sorting the BAM file by name (htseq default: name)

date
echo "Sorting Started"

for i in $(ls *.bam | sed -e 's/_Aligned.out.bam//' | sort -u)

do

	samtools sort\
	-@ thr\
	-n\
	${i}_Aligned.out.bam\
	-o ${i}.bam

	rm -rf ${i}_Aligned.out.bam

done


echo "Sorting Done"
date

gtf=$(if [[ $org == "hg38" ]]; then
        printf "/Homo_sapiens.GRCh38.99.gtf"
else
        printf "/Mus_musculus.GRCm38.99.gtf"
fi)

DIR=$(pwd | sed -e 's/BAM//' | awk '{print $0"/htseq/"}')

OUTPUT=$(for i in $(ls *ReadsPerGene.out.tab | head -n1 | sed -e 's/_ReadsPerGene.out.tab//')

do
	grep "ENS" ${i}_ReadsPerGene.out.tab | cut -f2-4 | numsum -c | awk '{if ($1>$2 && $1>$3) Stranded = "no"; else if ($2>$1 && $2>$3) Stranded = "yes"; else Stranded = "reverse"; print Stranded;}'
done)

echo "Stranded-ness : ${OUTPUT}"

GTF=$path$org$gtf

echo "HTSeq Started"

date
parallel -j $thr 'htseq-count -s {1} {2} {3} > {4}{2}_htseq.txt' :::  $OUTPUT ::: *.bam ::: $GTF ::: $DIR
date

echo "HTSeq Done"

# Preparing htseq files for further use

cd ../htseq/

GenLen=$path$org/$org'_GeneLength_kb.tsv'
../process-count.py $GenLen
rm *.bam_htseq.txt

echo "RNASeq pipeline complete"
