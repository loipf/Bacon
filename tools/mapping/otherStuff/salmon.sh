#!/bin/bash
#qsub -N salmon -o /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/RNASeq/GSE60217/salmon/salmon_info.out -e /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/RNASeq/GSE60217/salmon/salmon_error.out /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/salmon.sh

salmon="/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/salmon-0.10.2/src/salmon"
index_path="/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/human_grch38_92/salmon_transcriptome"
path="/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/RNASeq/GSE60217"

for file in $path/*.fastq
do
	base=$(basename $file | tail -c 17 | head -c 10)
	fileName="$path/salmon/$base"

	# Salmon
	mkdir $fileName
	$salmon quant -l A -i $index_path/hsapiens_index -r $file -p 8 --seqBias --gcBias -g $index_path/transcript_mapping.txt -o $fileName
done

