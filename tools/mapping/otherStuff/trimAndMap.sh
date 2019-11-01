#!/bin/bash

# program paths
trimmomatic='java -jar /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/Trimmomatic-0.38/trimmomatic-0.38.jar'
fastqc='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/FastQC/fastqc'
STAR='/home/proj/biosoft/software/STAR_2.3.0e.Linux_x86_64/STAR'
NextGenMap="/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm"
HISAT='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/hisat2-2.1.0/hisat2'
samtools='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/samtools-1.8/samtools'
bamstats='java -jar /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/BAMStats-1.25/BAMStats-GUI-1.25.jar'
FeatureCount='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/subread-1.4.6-p5-source/bin/featureCounts'
stringtie='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/stringtie-1.3.4/stringtie'
stringtieToCM='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/stringtieOut_to_countM.py'
DGEAnalyzer='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/DGEAnalyzer/DGEAnalyzeR.R'

# input values
ref_genome=$1
HISATindex=$2 # '/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/human_grch38_92/hisat_indices/'
STARindex=$2 # '/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/human_grch38_92/star_indices_overhang100/' # 100 is fine in most cases
gff=$3 #/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/RNASeq/GSE87554/GRCh38_latest_genomic.gff
gtf=$4 #/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/human_grch38_92/Homo_sapiens.GRCh38.92.gtf

#input
#output
#minlen

	
# trimming

#loop
	#select SE/PE	
	$trimmomatic SE -phred33 $input $output LEADING:25 TRAILING:25 SLIDINGWINDOW:10:25 MINLEN:$minlen
	$trimmomatic PE -phred33 $input -baseout $output LEADING:25 TRAILING:25 SLIDINGWINDOW:10:25 MINLEN:$minlen

	$fastqc
#end loop


#select mapper
	# STAR aligner, old
	$STAR --genomeDir /home/proj/biosoft/GENOMIC/HUMAN/STAR_INDEX --readFilesIn $inputFastq --readFilesCommand zcat --runThreadN 6

	# nextgenmap SE
	$NextGenMap -q $trimOut -r $ref_genome -o $outputSam -t 4
	# nextgenmap PE
	$NextGenMap-utils interleave -1 $trimOut1 -2 $trimOut2 -o $trimCombined
	$NextGenMap -q $trimCombined -p -r $ref_genome -o $outputSam -t 4
	
	# HISAT2 SE  p: #threads, q: input fastq, dta: tailord trans assembler, unpaired -U
	$HISAT -p 5 -q --dta --new-summary -x $HISATindex -U $inputFastq -S $ouputSam
	# HISAT2 PE  p: #threads, q: input fastq, dta: tailord trans assembler, paired -1 -2
	$HISAT -p 5 -q --dta --new-summary -x $HISATindex -1 $inputFastq1 -2 $inputFastq1 -S $ouputSam
	# make .sam to .bam and sort, does not work on bioclient
	$samtools sort -@ 4 -o $ouputBam $inputSam

	# STAR SE, only works on bioclient, otherwise too slow, outputName=$outputBam'Aligned.sortedByCoord.out.bam'
	$STAR --runThreadN 5 --genomeDir $STARindex --readFilesIn $inputFastq --outFileNamePrefix $outputBam --outSAMtype BAM SortedByCoordinate
	# STAR PE, only works on bioclient, otherwise too slow, outputName=$outputBam'Aligned.sortedByCoord.out.bam'
	$STAR --runThreadN 5 --genomeDir $STARindex --readFilesIn $inputFastq1 $inputFastq2 --outFileNamePrefix $outputBam --outSAMtype BAM SortedByCoordinate


	rm $sam
#end select


#select count tool
	#FeatureCount, -g gene only necessary for nxtgenmap
	$FeatureCount -a $gff -o $countMatrix -T 5 -g gene $inputSam 
	
	# stringTie, e: estimates abundance of ref transcript; bam file has to sorted
	$stringtie $inputBam -p 5 -G $gtf -e -o $outputGtf
	# converts all .gtf files with stringtie in name to countmatrix 	
	$stringtieToCM -i $outputGtfDirectory -g mapper_stringtie.txt


#end select


#DGE analysis R
Rscript $DGEAnalyzer







