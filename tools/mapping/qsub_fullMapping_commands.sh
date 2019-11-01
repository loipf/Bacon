#!/bin/bash



qsub -N fullMap  -l vf=20480M -o /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/rawdata/RNASeq/GSE60217/fullMap_info.out -e /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/rawdata/RNASeq/GSE60217/fullMap_error.out /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/fullMapRun2.sh




#qsub -N hisatMap -P long_proj -o /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/rawdata/RNASeq/GSE87554/fullMap_info.out -e /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/rawdata/RNASeq/GSE87554/fullMap_error.out /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/trimAndMap_full.sh

#qsub -N fullMapping -P long_proj -o /home/l/loipfingers/Desktop/2testFolder/qsub_mapInfo.out -e /home/l/loipfingers/Desktop/2testFolder/qsub_mapError.out bash /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/trimAndMap_full.sh human /home/l/loipfingers/Desktop/2testFolder/SRR4340399.fastq


#fastqFile=$1
#outputFolder="$(dirname "$fastqFile")"

#specify space and time
#does not work 
#qsub -N fullMapping -o $outputFolder/$fastqFile'_mapInfo.out' -e $outputFolder/$fastqFile'_mapError.out' /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/trimAndMap_full.sh human $fastqFile


#qsub -N fullMapping -o /home/l/loipfingers/Desktop/2testFolder/qsub_mapInfo.out -e /home/l/loipfingers/Desktop/2testFolder/qsub_mapError.out /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/trimAndMap_full.sh


#qsub -N hisatMap -o /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/RNASeq/GSE60217/hisat2_stringtie/fullMap_info.out -e /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/RNASeq/GSE60217/hisat2_stringtie/fullMap_error.out /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/trimAndMap_full.sh





#qsub -N hisatMap /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/trimAndMap_full.sh







