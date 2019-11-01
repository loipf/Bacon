#!/bin/bash


outputFolder="/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/rawdata/RNASeq/mirtrap"

# qsub all jobs in grid

for folder in $outputFolder/*
do
	for dir in $folder/*fw.fastq.gz
	do
		fileName="$(basename $dir '_fw.fastq.gz' )"
		dir2="$(dirname $dir)"
		dir2=$dir2"/"$fileName"_rw.fastq.gz"
		
		bash /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/mapping/applyTrim.sh $dir $dir2

		#qsub -N map  -l vf=20480M -o $outputFolder/info_$fileName.out -e $outputFolder/error_$fileName.out /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/mapping/applyTrim.sh $dir $dir2
		
	done
done



















