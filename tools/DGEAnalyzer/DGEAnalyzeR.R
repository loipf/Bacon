#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
dge_file = args[1]
targets_file = args[2]
design = args[3]
control = args[4]
condotion = args[5]

library(limma)
library(edgeR)
library(DESeq2)
library(KEGGREST)
library(biomaRt)
library(LSD)
library(pheatmap)
library(hexbin)
library(plot3D)

# load data of an RNAseq experiment and store it to an internal expression data object
# 
# dge_file: path to the DGE Matrix File ("DGE_Matrix.txt" by default)
# targets_file: path to targets file containing information about the experiment (e.g. group descriptors; "Targets.txt" by default)
# rowMin: minimal number of counts within one row (gene); the default value is 2
# 
# value: a list with three elements: "Targets" contains experimental information, rowMin the minimal number of counts within one row,
# "Matrix" the raw count matrix
loadData <- function(dge_file="DGE_Matrix.txt",targets_file="Targets.txt",rowMin=2){
	dge = read.csv(dge_file,header=T,sep='\t')
	# replace any '+' in the targets file by '.' ('+' in a name cannot be handled in R)
	targets = as.data.frame(apply(read.csv(targets_file,header=T,sep='\t'),2,function(x){x=gsub("\\+","\\.",x)}),stringsAsFactors=F)
	# if there are any fields with label 'Control', change these to 'control' (the default value for reference samples)
	for(i in 1:nrow(targets)){
	  for(j in 1:ncol(targets)){
	    if(targets[i,j] == "Control") targets[i,j] = "control"
	  }
	}
	colnames(targets)[1] = "Sample"
	rownames(targets) = targets$Sample
	
	rownames(dge) = dge$GENE
	dge = dge[,colnames(dge)!='GENE']
	dge = dge[,targets$Sample]
	dge = dge[apply(dge,1,sum) >= rowMin,colnames(dge) %in% targets$Sample]
	
	return(list(Targets=targets,rowMin=rowMin,Matrix=dge))
}


# transforms the raw count matrix of an expression data object by different methods (log, rlog, VST and CPM)
# 
# data: internal expression data object
# design: name of the column in the targets file that contains group information
# 
# value: internal expression data object with additional elements for each transformation and the design column
transformDGE <- function(data,design){
	data$log = log2(data$Matrix+1)
	
	dds <- DESeqDataSetFromMatrix(countData=data$Matrix,colData=data$Targets,as.formula(paste0("design=~",design)))
	data$rlog = assay(rlog(dds,blind=F))
	data$vst = assay(vst(dds,blind=F))
	data$cpm = cpm(data$Matrix,log=T)
	data$design = design
	
	return(data)
}


# rename samples by using informative columns from the targets information; if the columns are not sufficient to provide unique
# names, an enumeration for identical names is added; this can be done before or after transformation
# 
# data: internal expression data object
# fields: vector of column names to be considered for sample renaming
# 
# value: internal expression data object with changed sample names
createNames = function(data,fields){
	ids = c()
	for(i in rownames(data$Target)){
		ids = c(ids,paste(data$Targets[i,fields],collapse='_'))
	}
	if(anyDuplicated(ids)){
	  warning('Given fields are not sufficient to provide unique sample names. Enumeration is being added.')
	  for(i in ids[!duplicated(ids)]){
  	  nr=1
  	  for(j in which(ids==i)){
  	    ids[j] = paste0(ids[j],"_",nr)
  	    nr = nr+1
  	  }
  	  nr=1
  	}
	}
	names(ids) = data$Targets$Sample
	
	rownames(data$Targets) = ids
	cn=ids[colnames(data$Matrix)]
	for(i in names(data)[!(names(data)%in%c('Targets','rowMin','design'))]){
		colnames(data[[i]]) = cn
	}
	
	return(data)
}



# test for differentially expressed genes in RNASeq data using the DESeq2 package; currently, this only works with single factor experiments
# 
# data: internal expression data object
# ref: String giving the identifier of the group to be treated as reference (default is "control")
# condition: String giving the identifier of the group to test differential expression for
# 
# value: DESeqDataSet object with results of differential expression analysis
DifferentialExpressionDESeq = function(data,ref="control",condition){
  design = data$design
  targets = data$Targets[data$Targets[[design]]==ref | data$Targets[[design]]==condition,]
  mat = data$Matrix[,rownames(targets)]
  mat = mat[apply(mat,1,sum) >= data$rowMin,]
  dds <- DESeqDataSetFromMatrix(countData=mat,colData=targets,as.formula(paste0("design=~",design)))
  dds[[design]] <- relevel(dds[[design]],ref=ref)
  dds <- DESeq(dds)
  dds <- lfcShrink(dds, coef=paste0(design,"_",condition,"_vs_",ref))
  
  return(dds)
}


DifferentialExpressionLimma = function(data){
  #dge <- DGEList(counts=data$Matrix)
  design = model.matrix(~1+as.factor(data$Targets[[data$design]]))
  #keep <- filterByExpr(dge,design)
  #dge <- dge[keep,,keep.lib.sizes=F]
  #dge <- calcNormFactors(dge)
  
  i <- voom(data$Matrix,design)
  
  return(eBayes(lmFit(i,design)))
}


volcanoplotter = function(diff_mat,cutoff=2){
  diff_mat = diff_mat[!is.na(diff_mat$padj),]
  color = c("black","blue","red")
  col.code = rep(1,nrow(diff_mat))
  for(i in 1:nrow(diff_mat)){
    if(diff_mat$log2FoldChange[i]<=-cutoff & diff_mat$padj[i]<=0.05) col.code[i]=2
    if(diff_mat$log2FoldChange[i]>=cutoff & diff_mat$padj[i]<=0.05) col.code[i]=3
  }
  
  plot(diff_mat$log2FoldChange,-log10(diff_mat$padj),pch=20,xlab="log2 foldchange",ylab="-log10(p-value)",col=color[col.code])
  abline(h=-log10(0.05),lty=2)
  abline(v=c(-cutoff,cutoff),lty=2,col=c("blue","red"))
}

#########################################################################################################################################
# program code

# DESeq2
dat = loadData(dge_file,targets_file)
dat = transformDGE(dat,design)
dat = createNames(dat,design)
diff = testDifferentialExpression(dat,control,condition)

write.table(diff,"diffGenes.tab",quote=F,sep="\t")


# limma
fit = DifferentialExpressionLimma(dat)

tiff("VolcanoPlot.tif",600,400)
  volcanoplotter(data.frame(log2FoldChange=fit$coefficients[,2],padj=p.adjust(fit$p.value,"fdr")))
dev.off()
write.table(topTable(fit),"topGenes.txt",quote=F,row.names=F,sep='\t')