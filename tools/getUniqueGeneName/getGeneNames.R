
library(data.table)

# input gene list and printed out table
# in return: inputId, humanId, mouseId
# if organism specified -> does only output one list

path="/home/stefan/bioinformatics/neap/getUniqueGeneName/"
path = "/Users/stefanloipfinger/Documents/GitHub/Bacon/tools/getUniqueGeneName/"




# input gene list
args = commandArgs(trailingOnly=TRUE)
geneList <- args[1]
outputFile <- args[2]
organism <- args[3] # human|mouse
#geneList <- c("PAX1,TZSD,ZZ3, UI,P53,ALMS1,A2m,AAR2,ABCC10,A3GALT2,ABCC6")

geneList <- strsplit(geneList,',')[[1]]

tableHuman = fread(paste0(path, "human_gene_names.txt"))
tableMouse = fread(paste0(path, "mouse_gene_names.txt"))


# extract HGNC gene name 
getGeneName <- function(g, org){
  g <- toupper(g)
  if(org=="human"){ df <- tableHuman}  else { df <- tableMouse}
  s = df[toupper(df$Symbol)==g]$Symbol
  if(identical(s, character(0))) {
    l <- grepl(paste0("^",g,"(\\||$)|\\|",g,"(\\||$)"),toupper(df$Synonyms), perl=TRUE)   # only get specifc synonym
    l <- df[l]
    if(nrow(l)==0){
      return('')
    } else {
      s <- l$Symbol
    }
  }
  return(s)
}


if(!is.na(organism)) {
  l <- sapply(geneList, getGeneName, organism)
  print(paste(l, collapse=","))
} else {
      out = data.frame(input=character(), humanGene=character(), mouseGene=character())
      humanGene = sapply(geneList, getGeneName, "human")
      mouseGene = sapply(geneList, getGeneName, "mouse")    
      
      out = as.data.frame(cbind(humanGene,mouseGene))
      out$input=rownames(out)
	rownames(out) <- c()
      out = out[,c(3,1,2)]
	if(outputFile=='') {  print(out) }
	else { write.table(out,outputFile, sep="\t", row.names=FALSE, quote=FALSE)
	}
}















