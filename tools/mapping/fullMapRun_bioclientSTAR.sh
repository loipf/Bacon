#!/bin/bash


# fastqFile=$1
# outputFolder="$(dirname "$fastqFile")"

# call with 
# screen -S starMap -dm bash -c 'script -c "bash ~/Desktop/neap/finishedScripts/mapping/fullMapRun_bioclientSTAR.sh" scriptFileSTAR.txt'

# iterate over all, CHANGE REFERENCES
# path to all fastq files
outputFolder="/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/rawdata/RNASeq/mirtrap/H103"

# qsub all jobs in grid
for dir in $outputFolder/*.fastq.gz
do
	fileName="$(basename $dir '.fastq.gz')"
	#echo $dir
	# HANDLE PAIRED READS !!!!
	bash /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/mapping/trimAndMap_bioclientSTAR.sh human $dir

done





