#!/usr/bin/env bash

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

path='/path/to/indexed/genome/'
STAR='_STARIndexed'
echo "Organism: $1" 
read org
echo "Cluster Threads: $2" 
read thr

## Mapping fastq using STAR-aligner

ulimit -n 999999
ulimit -f unlimited

echo "Mapping Started"
date

if [ `ls *_1.fastq.gz 2>/dev/null | wc -l` -gt 0 ] ; then

        for i in $(ls *.fastq.gz | sed -e 's/_1.fastq.gz//' -e 's/_2.fastq.gz//' | sort -u); do

        STAR \
	--runThreadN $thr\
	--genomeDir $path$org$STAR\
	--readFilesCommand zcat\
	--readFilesIn ${i}_1.fastq.gz ${i}_2.fastq.gz\
	--outSAMattributes All\
	--outSAMtype BAM Unsorted\
	--quantMode GeneCounts\
	--outFileNamePrefix ../BAM/${i}_
        
done

else

       for i in $(ls *.fastq.gz | sed -e 's/.fastq.gz//' | sort -u)

do

	STAR \
	--runThreadN 112\
	--genomeDir $path$org$STAR\
	--readFilesCommand zcat\
	--readFilesIn ${i}.fastq.gz\
	--outSAMattributes All\
	--outSAMtype BAM Unsorted\
	--quantMode GeneCounts\
	--outFileNamePrefix ../BAM/${i}_

done

fi

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
parallel -j $thr 'htseq-count --additional-attr=gene_name -s {1} -f bam {2} {3} > {4}{2}_htseq.txt' :::  $OUTPUT ::: *.bam ::: $GTF ::: $DIR
date

echo "HTSeq Done"

# Preparing htseq files for further use

cd ../htseq/

for i in $(ls *.txt | sed -e 's/.bam_htseq.txt//' | sort -u)

do
	head -n -5 ${i}.bam_htseq.txt > ${i}.htseq.txt
	sed -i 1i"GeneID\tGeneName\t${i%.htseq.txt}" ${i}.htseq.txt
done

tmp=$(ls *.htseq.txt | head -n1 | sed -e 's/.htseq.txt//')
cut -f1,2 $tmp.htseq.txt > tmp1.txt

paste *.htseq.txt | awk '{ for (i=3;i<=NF;i+=3) {printf "%s ",$i} ;print ""}' | sed 's/ /\t/g' > tmp2.txt

dir=$(dirname $PWD | sed 's:.*/::')
paste tmp1.txt tmp2.txt > $dir.htseq.txt

sed -e s/GeneID//g $dir.htseq.txt | sed -e s/GeneName//g > ${dir}_RawCounts.tsv

echo "RNASeq pipeline complete"
