library(limma)
library(DESeq2)
library(Biobase)
library(GEOquery)

# load series and platform data from GEO
gset <- getGEO("GSE9874", GSEMatrix =TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

dge = exprs(gset)
genenames = featureData(gset)@data[,"Gene Symbol"]
names(genenames) = featureData(gset)@data[,"ID"]
rownames(dge) = genenames[rownames(dge)]
t = phenoData(gset)@data$title
t2 = as.character(t)
t2[grep("with atherosclerosis",t)] = "disease"
t2[grep("without atherosclerosis",t)] = "control"
targets = data.frame(Barcode=rownames(phenoData(gset)@data),Genotype=t2)

dat = list(Targets=targets,rowMin=2,Matrix=dge,log=log2(dge),design="Genotype")
design = model.matrix(~1+as.factor(t2))

fit = lmFit(dat$log,design)
fit = eBayes(fit)
tiff("VolcanoPlot.tif",600,400)
  volcanoplot(fit,coef=2)
dev.off()
write.table(topTable(fit),"topGenes.txt",quote=F,row.names=F,sep='\t')

control = dat$log[,targets$Barcode[targets$Genotype=="control"]]
disease = dat$log[,targets$Barcode[targets$Genotype=="disease"]]
p = c()
for(i in 1:nrow(dge)){
  p[i] = t.test(control[i,],disease[i,],paired=T)$p.value
}
p.adj = p.adjust(p,method="fdr")

plot(apply(disease,1,mean)-apply(control,1,mean),-10*log10(p.adj),pch=20,xlab="Log2 Fold Change",ylab="-log10(P-value)",main="Differential Gene Expression")

