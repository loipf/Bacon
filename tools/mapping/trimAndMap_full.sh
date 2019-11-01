#!/bin/bash

#checklist:
# single: trim, hisat, samtools, star, kallisto
# stringtie with hisat, star nextgenmap
# featureCount test R script !!
# paired: 



# program paths
trimmomatic='java -jar /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/Trimmomatic-0.38/trimmomatic-0.38.jar'
fastqc='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/FastQC/fastqc'
#STAR='/home/proj/biosoft/software/STAR_2.3.0e.Linux_x86_64/STAR'
STAR='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/STAR-2.6.0c/bin/Linux_x86_64_static/STAR'
NextGenMap='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm'
HISAT='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/hisat2-2.1.0/hisat2'
samtools='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/samtools-1.8/samtools' # doesnt run on bioclient
#samtools='/home/proj/biosoft/software/samtools-1.2/samtools' # doesnt run on grid
bamstats='java -jar /home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/BAMStats-1.25/BAMStats-GUI-1.25.jar'
FeatureCount='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/subread-1.4.6-p5-source/bin/featureCounts'
FeatureCountToCM='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/subreadOut_to_countM.R'
stringtie='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/stringtie-1.3.4/stringtie'
stringtieToCM='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/stringtieOut_to_countM.py'
subread='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/subread-1.4.6-p5-source/bin/featureCounts'
kallisto='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/kallisto_v0.44.0'
salmon='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/tools/salmon-0.11.0'
DGEAnalyzer='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/DGEAnalyzer/DGEAnalyzeR.R'
Rpath='/home/proj/biosoft/software/R-3.4.3/bin/Rscript'
ensembl2gene='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/ensembl_to_geneName.R'


# input values
organism=$1
fastq1=$2
#organism='human'
#fastq1='/home/l/loipfingers/Desktop/2testFolder/SRR4340399.fastq'
fastq2=$3
trimMinLen=50
threadNo=3

outputFolder="$(dirname "$fastq1")"
fileName1="$(basename $fastq1)"
fileName1="${fileName1%.*}"
singleEnd=true


# check if paired end
if [ $# -ge 3 ]
then
	fileName1="$(basename $fastq1 '.fastq.gz' )"
	fileName2="$(basename $fastq2 '.fastq.gz' )"
	singleEnd=false
fi


# load according index files
if [ $organism = 'human' ]
then
	ref_genome='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/human_grch38_92/genome.fa'
	HISATindex='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/human_grch38_92/hisat_indices/genome'
	STARindex='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/human_grch38_92/star_indices_overhang100/' # 100 is fine in most cases
	kallistoIndex='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/human_grch38_92/kallisto_indices/kallistoIndex.idx'
	salmonIndex='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/human_grch38_92/salmon_transcriptome'
	#gff= #/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/RNASeq/GSE87554/GRCh38_latest_genomic.gff
	gtf='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/human_grch38_92/Homo_sapiens.GRCh38.92.gtf'

elif [ $organism = 'mouse' ]
then
	ref_genome='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/mouse_grcm38_92/genome.fa'
	HISATindex='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/mouse_grcm38_92/hisat_indices/genome'
	STARindex='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/mouse_grcm38_92/star_indices_overhang100/' # 100 is fine in most cases
	kallistoIndex='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/mouse_grcm38_92/kallisto_indices/kallistoIndex.idx'
	salmonIndex='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/mouse_grcm38_92/salmon_transcriptome'
	#gff #/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/RNASeq/GSE87554/GRCh38_latest_genomic.gff
	gtf='/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/mouse_grcm38_92/Mus_musculus.GRCm38.92.gtf'
fi


#:'  # comment start

# trimming
if [ "$singleEnd" = true ]
then
	echo "trim"	
	$trimmomatic SE -threads $threadNo -phred33 $fastq1 $outputFolder/$fileName1'_trimmed.fastq.gz' LEADING:25 TRAILING:25 SLIDINGWINDOW:10:25 MINLEN:$trimMinLen
else
	if [ ! -e $fileName1'_trimmed_paired.fastq.gz' ]	
		echo "trim"
		$trimmomatic PE -threads $threadNo -phred33 $fastq1 $fastq2 $outputFolder/$fileName1'_trimmed_paired.fastq.gz' $outputFolder/$fileName1'_trimmed_unpaired.fastq.gz' $outputFolder/$fileName2'_trimmed_paired.fastq.gz' $outputFolder/$fileName2'_trimmed_unpaired.fastq.gz' LEADING:25 TRAILING:25 SLIDINGWINDOW:10:25 MINLEN:$trimMinLen
		rm -f $outputFolder/$fileName1'_trimmed_unpaired.fastq.gz'
		rm -f $outputFolder/$fileName2'_trimmed_unpaired.fastq.gz'
	fi
fi

# $fastqc



#select mapper
if [ "$singleEnd" = true ]
then
	echo "map"
	# nextgenmap SE --bam 
	$NextGenMap -t $threadNo -q $outputFolder/$fileName1'_trimmed.fastq.gz' -r $ref_genome -o /home/proj/biocluster/tmp/$fileName1'_nextgenmap.sam'

	# HISAT2 SE  p: #threads, q: input fastq, dta: tailord trans assembler, unpaired -U
	$HISAT -p $threadNo -q --dta --new-summary -x $HISATindex -U $outputFolder/$fileName1'_trimmed.fastq.gz' -S /home/proj/biocluster/tmp/$fileName1'_hisat.sam'

	# STAR SE, only works on bioclient, otherwise too slow, outputName=$outputBam'Aligned.sortedByCoord.out.bam'
	# $STAR --runThreadN $threadNo --genomeDir $STARindex --readFilesIn $outputFolder/$fileName1'_trimmed.fastq.gz' --readFilesCommand zcat --outFileNamePrefix $outputFolder/$fileName1
else
	# nextgenmap PE
	$NextGenMap -t $threadNo -1 $outputFolder/$fileName1'_trimmed_paired.fastq.gz' -2 $outputFolder/$fileName2'_trimmed_paired.fastq.gz' -r $ref_genome -o /home/proj/biocluster/tmp/$fileName1'_nextgenmap.sam'

	# HISAT2 PE  p: #threads, q: input fastq, dta: tailord trans assembler, paired -1 -2
	$HISAT -p $threadNo -q --dta --new-summary -x $HISATindex -1 $outputFolder/$fileName1'_trimmed_paired.fastq.gz' -2 $outputFolder/$fileName2'_trimmed_paired.fastq.gz' -S /home/proj/biocluster/tmp/$fileName1'_hisat.sam'
		
	# STAR PE, only works on bioclient, otherwise too slow, outputName=$outputBam'Aligned.sortedByCoord.out.bam'
	# $STAR --runThreadN $threadNo --genomeDir $STARindex --readFilesIn $outputFolder/$fileName1'_trimmed_paired.fastq.gz' $outputFolder/$fileName2'_trimmed_paired.fastq.gz' --readFilesCommand zcat --outFileNamePrefix $fileName1
	

fi


# post-pre processing

	echo "post-process"
	# nxtgenmap: make .sam to .bam and sort
	# $samtools view -u /home/proj/biocluster/tmp/$fileName1'_nextgenmap.sam' | samtools sort -@ $threadNo -o $outputFolder/$fileName1'_nextgenmap.bam'

	$samtools sort -@ $threadNo -o $outputFolder/$fileName1'_nextgenmap.bam' /home/proj/biocluster/tmp/$fileName1'_nextgenmap.sam'

	rm -f /home/proj/biocluster/tmp/$fileName1'_nextgenmap.sam'


	# from hisat2: make .sam to .bam and sort
	# $samtools view -u /home/proj/biocluster/tmp/$fileName1'_hisat.sam' | samtools sort -@ $threadNo -o $outputFolder/$fileName1'_hisat.bam'

	$samtools sort -@ $threadNo -o $outputFolder/$fileName1'_hisat.bam' /home/proj/biocluster/tmp/$fileName1'_hisat.sam'

	rm -f /home/proj/biocluster/tmp/$fileName1'_hisat.sam'

	# star output to bam and rename, rm unnecessary files from star
	# $samtools view -u $outputFolder/$fileName1'Aligned.out.sam' | samtools sort -@ $threadNo -o $outputFolder/$fileName1'_star.bam'
	#rm -f $outputFolder/$fileName1'Aligned.out.sam'
	#rm -f $outputFolder/$fileName1'SJ.out.tab'
	#rm -f $outputFolder/$fileName1'Log.out'
	#rm -f $outputFolder/$fileName1'Log.progress.out'
	#rm -rf $outputFolder/$fileName1'_STARtmp'

	

# create folders for gtf output if they dont exist
	mkdir -p $outputFolder/nextgenmap_subread
	mkdir -p $outputFolder/nextgenmap_stringtie
	#mkdir -p $outputFolder/star_subread
	#mkdir -p $outputFolder/star_stringtie
	mkdir -p $outputFolder/hisat_subread
	mkdir -p $outputFolder/hisat_stringtie
	mkdir -p $outputFolder/kallisto
	mkdir -p $outputFolder/salmon





#select count tool
	echo "count tool"
	#FeatureCount, -g gene only necessary for nxtgenmap !!!!!!!!!!!!!!!!!!!!!
	### $FeatureCount -a $gtf -o $countMatrix -T 5 -g gene $inputSam    -g gene !! -t exon ?
	# $FeatureCount -a $gtf -o $outputFolder/star_subread/$fileName1'_star_subread.txt' -T $threadNo $outputFolder/$fileName1'_star.bam'
	$FeatureCount -a $gtf -o $outputFolder/hisat_subread/$fileName1'_hisat_subread.txt' -T $threadNo $outputFolder/$fileName1'_hisat.bam'
	$FeatureCount -a $gtf -o $outputFolder/nextgenmap_subread/$fileName1'_nextgenmap_subread.txt' -T $threadNo $outputFolder/$fileName1'_nextgenmap.bam' 

#' # comment end

	
	# stringTie, e: estimates abundance of ref transcript; bam file has to sorted
	# $stringtie $outputFolder/$fileName1'_star.bam' -p $threadNo -G $gtf -e -o $outputFolder/star_stringtie/$fileName1'_star_stringtie.gtf'
	$stringtie $outputFolder/$fileName1'_hisat.bam' -p $threadNo -G $gtf -e -o $outputFolder/hisat_stringtie/$fileName1'_hisat_stringtie.gtf'
	$stringtie $outputFolder/$fileName1'_nextgenmap.bam' -p $threadNo -G $gtf -e -o $outputFolder/nextgenmap_stringtie/$fileName1'_nextgenmap_stringtie.gtf'




# alternative mapper

mkdir -p $outputFolder/kallisto/$fileName1
if [ "$singleEnd" = true ]
then
	echo "kallisto"
	# kallisto SE -b 100 bootstrap, --single need -l avg read length, -s standard deviation, -t threads
	gunzip < $outputFolder/$fileName1'_trimmed.fastq.gz' > /home/proj/biocluster/tmp/$fileName1'_trimmed.fastq'
	singleEndPar=$(bash $kallisto/kallisto_SE_fastqMeanLen.sh /home/proj/biocluster/tmp/$fileName1'_trimmed.fastq') # get .fastq mean length and stdDev for kallisto SE run 
	rm /home/proj/biocluster/tmp/$fileName1'_trimmed.fastq'
	$kallisto/kallisto quant -i $kallistoIndex -o $outputFolder/kallisto/$fileName1 -b 100 --single $singleEndPar -t $threadNo $outputFolder/$fileName1'_trimmed.fastq.gz'
	#rm $outputFolder/kallisto/$fileName1/abundance.h5  # r package h5 input doesnt work
	
	$salmon/bin/salmon quant -l A -i $salmonIndex/hsapiens_index -r $outputFolder/$fileName1'_trimmed.fastq.gz' -p $threadNo --seqBias --gcBias -g $salmonIndex/transcript_mapping.txt -o $outputFolder/salmon/$fileName1

else
	# kallisto PE -b 100 bootstrap, --single need -l avg read length, -s standard deviation, -t threads
	$kallisto/kallisto quant -i $kallistoIndex -o $outputFolder/kallisto/$fileName1 -b 100 -t $threadNo $outputFolder/$fileName1'_trimmed_paired.fastq.gz' $outputFolder/$fileName2'_trimmed_paired.fastq.gz'
	#rm $outputFolder/kallisto/$fileName1/abundance.h5  # r package h5 input doesnt work

	$salmon/bin/salmon quant -l A -i $salmonIndex/hsapiens_index -1 $outputFolder/$fileName1'_trimmed_paired.fastq.gz' -2 $outputFolder/$fileName2'_trimmed_paired.fastq.gz' -p $threadNo --seqBias --gcBias -g $salmonIndex/transcript_mapping.txt -o $outputFolder/salmon/$fileName1

fi

#'


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






# create gene count matrix from kallisto transcript output -> error !!
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




