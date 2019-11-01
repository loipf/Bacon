#!/bin/bash


#fastqFile=$1
#outputFolder="$(dirname "$fastqFile")"

Rpath='/home/proj/biosoft/software/R-3.4.3/bin/Rscript'
ensembl2gene='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/ensembl_to_geneName.R'
salmon='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/salmon-0.11.0'
uniqueGeneNames='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/uniqueGeneNamesMatrix.R'



# iterate over all, CHANGE REFERENCES
# path to all fastq files
outputFolder="/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/rawdata/RNASeq/GSE60217_2"


for dir in $outputFolder/*.fastq.gz
do
	#fileName="$(basename $dir '.fastq.gz')"
	#echo $fileName	
	bash /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/trimAndMap_full.sh human $dir
	

done





######################
### call after all .fastq files are run through

# create gene count matrix from kallisto transcript output -> error !!
$Rpath $kallisto/kallisto_trans2gene.R $outputFolder/kallisto


# create count matrix from salmon files
bash $salmon/salmon_multiple2cm.sh $outputFolder


# transform count matrices with ensembl ids to gene names
echo "transforming transcripts to genes"
$Rpath $ensembl2gene $outputFolder/hisat_subread_matrix.txt
$Rpath $ensembl2gene $outputFolder/star_subread_matrix.txt
$Rpath $ensembl2gene $outputFolder/nextgenmap_subread_matrix.txt
$Rpath $ensembl2gene $outputFolder/hisat_stringtie_matrix.txt
$Rpath $ensembl2gene $outputFolder/star_stringtie_matrix.txt
$Rpath $ensembl2gene $outputFolder/nextgenmap_stringtie_matrix.txt
$Rpath $ensembl2gene $outputFolder/salmon_matrix.txt
$Rpath $ensembl2gene $outputFolder/kallisto_matrix.txt


$Rpath $uniqueGeneNames $outputFolder


