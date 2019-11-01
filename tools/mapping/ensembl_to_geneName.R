
# converts count matrix with ensembl ids to gene names

args = commandArgs(trailingOnly=TRUE)
inputFile <- args[1]

input <- read.table(inputFile, header=TRUE, sep="\t")
input <- data.frame(lapply(input, as.character), stringsAsFactors=FALSE)
geneTable <- read.table('/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/genomes/human_grch38_92/human_ensembl_tran_gene_name.txt',  header=TRUE, sep="\t")
geneTable <- data.frame(lapply(geneTable, as.character), stringsAsFactors=FALSE)


# get gene names of ensembl table
getName <- function(x) {
  row = geneTable[geneTable$ens_gene==x,]
  if(nrow(row)==0) {
	return(x)
  }
  else{
    return(row$ext_gene[1])
  }
}

geneNames <- sapply(input$GENE, function(x) getName(x) )
input$GENE <- geneNames

write.table(input,inputFile, quote = FALSE, sep="\t", row.names=FALSE)



