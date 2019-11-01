Rpath='/home/proj/biosoft/software/R-3.4.3/bin/Rscript'
ensembl2gene='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/ensembl_to_geneName.R'
salmon='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/salmon-0.11.0'


# iterate over all, CHANGE REFERENCES
# path to all fastq files
outputFolder="/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/rawdata/RNASeq/GSE87554"

counter=1

for dir in $outputFolder/*.fastq.gz
do
	fileName="$(basename $dir '.fastq.gz')"
	echo $dir
	echo $fileName
	#bash /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/trimAndMap_full.sh human $dir > $outputFolder/$fileName'_log'$counter.txt
	#counter=$(($counter + 1))

done

