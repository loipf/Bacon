#!/bin/bash


# program paths
trimmomatic='java -jar /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/Trimmomatic-0.38/trimmomatic-0.38.jar'

fastq1=$1
fastq2=$2

outputFolder="$(dirname "$fastq1")"

fileName1="$(basename $fastq1 '.fastq.gz' )"
#fileName1="${fileName1%.*}"


fileName2="$(basename $fastq2 '.fastq.gz' )"
#fileName2="${fileName2%.*}"

$trimmomatic PE -threads 3 -phred33 $fastq1 $fastq2 -baseout $outputFolder/$fileName1'_trimmed_paired.fastq.gz' $outputFolder/$fileName1'_trimmed_unpaired.fastq.gz' $outputFolder/$fileName2'_trimmed_paired.fastq.gz' $outputFolder/$fileName2'_trimmed_unpaired.fastq.gz' LEADING:25 TRAILING:25 SLIDINGWINDOW:10:25 MINLEN:50
rm -f $outputFolder/$fileName1'_trimmed_unpaired.fastq.gz'
rm -f $outputFolder/$fileName2'_trimmed_unpaired.fastq.gz'







