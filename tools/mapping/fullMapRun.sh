#!/bin/bash


# fastqFile=$1
# outputFolder="$(dirname "$fastqFile")"

Rpath='/home/proj/biosoft/software/R-3.4.3/bin/Rscript'
ensembl2gene='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/ensembl_to_geneName.R'
salmon='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/salmon-0.11.0'
uniqueGeneNames='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/uniqueGeneNamesMatrix.R'
FeatureCountToCM='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/subreadOut_to_countM.R'
stringtieToCM='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/stringtieOut_to_countM.py'


# iterate over all, CHANGE REFERENCES
# path to all fastq files
outputFolder="/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/rawdata/RNASeq/mirtrap"


for folder in $outputFolder/H*    # CHANGE H TO M AND ALSO IN LINE 33
do
	# qsub all jobs in grid
	for dir in $folder/*fw.fastq.gz  # $outputFolder/*fastq.gz
	do
		# fileName="$(basename $dir '.fastq.gz')"  #SE

		fileName="$(basename $dir '_fw.fastq.gz' )"	# PE
		dir2="$(dirname $dir)"
		dir2=$dir2"/"$fileName"_rw.fastq.gz"
		
		#echo $dir   # 10240 10 gb
		#bash /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/trimAndMap_full.sh human $dir
		qsub -N map  -l vf=4000M,h_rt=4:00:00 -o $outputFolder/info_$fileName.out -e $outputFolder/error_$fileName.out /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/mapping/trimAndMap_full.sh human $dir $dir2
		#bash /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/mapping/trimAndMap_full.sh human $dir $dir2
	done
done





######################
### call after all .fastq files are run through


# converts feature count output to countmatrix
# $Rpath $FeatureCountToCM $outputFolder/star_subread/ $outputFolder/star_subread_matrix.txt
# $Rpath $FeatureCountToCM $outputFolder/hisat_subread/ $outputFolder/hisat_subread_matrix.txt
# $Rpath $FeatureCountToCM $outputFolder/nextgenmap_subread/ $outputFolder/nextgenmap_subread_matrix.txt

# converts all .gtf files with stringtie in name to countmatrix 	
# python $stringtieToCM -i $outputFolder/nextgenmap_stringtie/ -g $outputFolder/nextgenmap_stringtie_matrix.txt
# python $stringtieToCM -i $outputFolder/star_stringtie/ -g $outputFolder/star_stringtie_matrix.txt
# python $stringtieToCM -i $outputFolder/hisat_stringtie/ -g $outputFolder/hisat_stringtie_matrix.txt

# create gene count matrix from kallisto transcript output 
# $Rpath $kallisto/kallisto_trans2gene.R $outputFolder/kallisto

# create count matrix from salmon files
# bash $salmon/salmon_multiple2cm.sh $outputFolder



# transform count matrices with ensembl ids to gene names
# echo "transforming transcripts to genes"
# $Rpath $ensembl2gene $outputFolder/hisat_subread_matrix.txt
# $Rpath $ensembl2gene $outputFolder/star_subread_matrix.txt
# $Rpath $ensembl2gene $outputFolder/nextgenmap_subread_matrix.txt
# $Rpath $ensembl2gene $outputFolder/hisat_stringtie_matrix.txt
# $Rpath $ensembl2gene $outputFolder/star_stringtie_matrix.txt
# $Rpath $ensembl2gene $outputFolder/nextgenmap_stringtie_matrix.txt
# $Rpath $ensembl2gene $outputFolder/salmon_matrix.txt
# $Rpath $ensembl2gene $outputFolder/kallisto_matrix.txt

# keep only unique gene names
# $Rpath $uniqueGeneNames $outputFolder




