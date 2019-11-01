#!/bin/bash

# program paths
trimmomatic='java -jar /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/Trimmomatic-0.38/trimmomatic-0.38.jar'
fastqc='/home/proj/biosoft/software/FastQC/fastqc'
STAR='/home/proj/biosoft/software/STAR_2.3.0e.Linux_x86_64/STAR'
bamstats='java -jar /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/BAMStats-1.25/BAMStats-GUI-1.25.jar'


# input values
base=$1
ref_genome=$2
library=$3
bwa_index=$4
SNP_file=$5

	
# trimming
$trimmomatic PE -basein $base\_S*_R1_001.fastq.gz -threads 10 -phred33 -baseout $base.trimmed.fastq.gz LEADING:25 TRAILING:25 SLIDINGWINDOW:10:25 MINLEN:75
$fastqc

# mapping (STAR aligner)
$STAR --genomeDir /home/proj/biosoft/GENOMIC/HUMAN/STAR_INDEX --readFilesIn my_reads.fastq --readFilesCommand zcat #--runThreadN 6
