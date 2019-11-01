#############################################################################################################################
# 
# Enrichment AnalyzeR v1.0
# Toolbox for the analysis of geneset enrichment from RNAseq derived count matrices
# (c) Maximilian Zwiebel, 2018
# 
#############################################################################################################################

library(KEGGREST)
library(biomaRt)
library(DESeq2)
library(globaltest)

# dataframe containing GO evidence codes and their meaning; should not be changed 
evcode = data.frame(Shortcut=c("EXP","IDA","IPI","IMP","IGI","IEP","ISS","ISO","ISA","ISM","IGC","IBA","IBD","IKR","IRD","RCA","TAS","NAS","IC","ND","IEA"),Description=c("Inferred from Experiment","Inferred from Direct Assay","Inferred from Physical Interaction","Inferred from Mutant Phenotype","Inferred from Genetic Interaction","Inferred from Expression Pattern","Inferred from Sequence or structural Similarity","Inferred from Sequence Orthology","Inferred from Sequence Alignment","Inferred from Sequence Model","Inferred from Genomic Context","Inferred from Biological aspect of Ancestor","Inferred from Biological aspect of Descendant","Inferred from Key Residues","Inferred from Rapid Divergence","Inferred from Reviewed Computational Analysis","Traceable Author Statement","Non-traceable Author Statement","Inferred by Curator","No biological Data available","Inferred from Electronic Annotation"))

# path to the executable GSEA .jar file
gsea_path = "C:/Program Files (x86)/GSEA/gsea-3.0.jar"
# path to the folder containing the most current gmt files of MSigDB
msig_path = "C:/Users/Maximilian/Desktop/RNAseq/MSigDB"
# path to the file containing the most current geneset annotations from Ensembl
ensembl_path = "C:/Users/Maximilian/Desktop/RNAseq/GO/Ensembl_all.txt"
# prefixes of the gmt files and names of the corresponding datasets; should only be changed if new sets are added
msig_ids = c("h.all.v","c2.cp.b","c3.mir.","c3.tft.")
msig_sets = c("Hallmark","Biocarta","miRNA_Targets","TranscriptionFactorTargets")


# update geneset databases and descriptors
# 
# organism: either human (default) or mouse; in case of updating MSigDB, both organisms are updated simultaniously
# db: one of either GO, Reactome, KEGG, MSigDB or all (default)
# 
# value: depending on the selected database, the global variables db_organism and db_organism_description are updated;
# the databases are organized as lists with the geneset identifiers as names and the corresponding gene names as list values,
# stored as vectors; for GO, the list contains three sublists for the different ontologies (GO_MF, GO_BP and GO_CC); for MSigDB
# the list contains several sublists, named according to the fields in msig_sets; the corresponding descriptions are stored
# in dataframes with geneset identifiers mapped to their descriptors if available; otherwithe, the dataframe only contains
# NAs; in order to maintain compatibility with the other functions, these datastructures should be preserved
updateDatabases = function(organism=c("human","mouse"),db=c("all","GO","Reactome","KEGG","MSigDB")){
  organism = match.arg(organism)
  db = match.arg(db)
  
  # GO, currently loaded from Ensembl, has to be updated
  if(db=="all" | db=="GO"){
    #if(organism == "human") o = "hsapiens" else o = "mmusculus"
    ens_data = read.csv(ensembl_path,header=T,sep="\t",stringsAsFactors=F)
    colnames(ens_data) = c("external_gene_name","go_id","name_1006","go_linkage_type","namespace_1003")
    ens_data = ens_data[ens_data$go_id!="",]
    
    GO = list(GO_MF=list(),GO_BP=list(),GO_CC=list())
    for(i in unique(ens_data$go_id)){
      term = unique(ens_data[ens_data$go_id == i,c("external_gene_name","go_linkage_type","namespace_1003")])
      if(term$namespace_1003[1]=="molecular_function") GO[[1]][[i]] = term[,1:2]
      else if(term$namespace_1003[1]=="biological_process") GO[[2]][[i]] = term[,1:2]
      else GO[[3]][[i]] = term[,1:2]
    }
    
    go_descr = unique(ens_data[,c("go_id","name_1006")])
    colnames(go_descr) = c("ID","Description")
    
    if(organism == "human"){
      assign("GO_human",GO,envir=.GlobalEnv)
      assign("GO_human_description",go_descr,envir=.GlobalEnv)
    } else{
      assign("GO_mouse",GO,envir=.GlobalEnv)
      assign("GO_mouse_description",go_descr,envir=.GlobalEnv)
    }
  }
  
  # Reactome; currently not implemented
  if(db=="all" | db=="Reactome"){
    
  }
  
  # KEGG; fully automatic download of the KEGG database using the KEGGREST package
  if(db=="all" | db=="KEGG"){
    if(organism == "human") id = "hsa" else id = "mmu"
    kegg = keggLink("pathway",id) # get mapping of KEGG pathway identifiers to genes
    kegg = substr(kegg,9,13)
    
    pathways = keggList("pathway") # get descriptors for pathway identifiers
    names(pathways) = substr(names(pathways),9,13)
    pathways = pathways[names(pathways)%in%kegg]
    genes = keggList(id) # get all genes from specified organism as listed in KEGG
    genes = strsplit(genes,",|;",perl=T)
    names(kegg) = sapply(names(kegg),function(x){x=genes[[x]][1]}) # replace KEGG gene identifiers by common names
    
    KEGG=list()
    for(i in unique(kegg)){ # parse geneset mapping table and construct databse list structure
      KEGG[[i]] = names(kegg)[kegg==i]
    }
    kegg_descr = data.frame(ID=names(pathways),Description=pathways,stringsAsFactors=F)
    no_descr = setdiff(names(KEGG),kegg_descr[,1])
    for(i in no_descr){
      kegg_descr = rbind(kegg_descr,data.frame(ID=i,Description=NA))
    }
    
    if(organism == "human"){
      assign("KEGG_human",KEGG,envir=.GlobalEnv)
      assign("KEGG_human_description",kegg_descr,envir=.GlobalEnv)
    } else{
      assign("KEGG_mouse",KEGG,envir=.GlobalEnv)
      assign("KEGG_mouse_description",kegg_descr,envir=.GlobalEnv)
    }
  }
  
  # MSigDB; updates both the human and mouse MSigDB using the files provided in msig_path; the mouse database is updated
  # by mapping to orthologues from Ensembl
  if(db=="all" | db=="MSigDB"){
    MSIG = list()
    for(i in list.files(msig_path)){
      msig = list()
      set = msig_sets[which(msig_ids==substr(i,1,7))]
      if(length(set)==0) next
      
      l=max(count.fields(paste0(msig_path,"/",i),sep = "\t"))
      tab = read.table(paste0(msig_path,"/",i),header=F,sep="\t",col.names=paste0("V",seq_len(l)),fill=T)
      for(j in 1:nrow(tab)){
        x = tab[j,3:ncol(tab)]
        msig[[as.character(tab[j,1])]] = x[x!=""]
      }
      MSIG[[set]] = msig
    }
    msig_descr = data.frame(NA,NA)
    assign("MSigDB_human",MSIG,envir=.GlobalEnv)
    assign("MSigDB_human_description",msig_descr,envir=.GlobalEnv)
    
    # orthology based geneset assignment for mouse data; only one2one matches are used
    ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
    ortho = getBM(attributes=c("external_gene_name","mmusculus_homolog_associated_gene_name","mmusculus_homolog_orthology_type"),filter="external_gene_name",values=unique(unlist(unlist(MSIG))),mart=ensembl)
    ortho = ortho[ortho$mmusculus_homolog_orthology_type=="ortholog_one2one",1:2]
    MSIG = lapply(MSIG,function(x){lapply(x,function(y){y[y%in%ortho[,1]]})})
    MSIG = lapply(MSIG,function(x){lapply(x,function(y){sapply(y,function(z){z=ortho[ortho[,1]==z,2]},USE.NAMES=F)})})
    assign("MSigDB_mouse",MSIG,envir=.GlobalEnv)
    assign("MSigDB_mouse_description",msig_descr,envir=.GlobalEnv)
  }
}


# display the evidence code dataframe on the console
showEvidenceCodes = function(){
  print(evcode)
}


# test for differentially expressed genes in RNASeq data using the DESeq2 package; currently, this only works with single factor experiments
# 
# data: internal expression data object
# ref: String giving the identifier of the group to be treated as reference (default is "control")
# condition: String giving the identifier of the group to test differential expression for
# 
# value: DESeqDataSet object with results of differential expression analysis
testDifferentialExpression = function(data,ref="control",condition){
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


expressionHeatmap = function(data,diff_mat,top=50){
  diff_mat = diff_mat[!is.na(diff_mat$padj) & diff_mat$padj<=0.05,]
  diff_mat = diff_mat[order(abs(diff_mat$log2FoldChange)),][(nrow(diff_mat)-top+1):nrow(diff_mat),]
  
  pheatmap(data$rlog[rownames(diff_mat),])
}

# auxiliary function to return a list of genesets from a certain database; should only be called internally from enrichment
# analysis functions to get the desired genesets
# 
# database: name of the database to use (one of either GO, Reactome, KEGG or MSigDB)
# organism: either human or mouse
# evidence: vector of GO evidence codes to restrict database entries onto; with NULL (the default), all evidence codes are used
# ontology: GO ontology to use (one of either cellular component, biological process or molecular function)
# 
# value: list of genesets from the selected database
getDB = function(database,organism,evidence,ontology){
  if(database=="GO"){
    if(organism=="human") db = GO_human else db = GO_mouse
    
    if(ontology=="molecular_function") db = db[["GO_MF"]]
    else if(ontology=="biological_process") db = db[["GO_BP"]]
    else db = db[["GO_CC"]]
    
    if(is.null(evidence)){
      db = lapply(db,function(x){unique(x[,1])})
    } else{
      db = lapply(db,function(x){unique(x[x$go_linkage_type %in% evidence,1])})
    }
  } else if(database=="Reactome"){
    if(organism=="human") db = Reactome_human else db = Reactome_mouse
  } else if(database=="KEGG"){
    if(organism=="human") db = KEGG_human else db = KEGG_mouse
  } else{
    if(organism=="human") db = MSigDB_human else db = MSigDB_mouse
  }
  
  return(db)
}


# auxiliary function to return the genesets descriptions belonging to a certain database; should only be called internally
# from enrichment analysis functions to get the necessary descriptors
# 
# database: name of the database to use (one of either GO, Reactome, KEGG or MSigDB)
# organism: either human or mouse
# 
# value: dataframe containing geneset descriptors for the selected database
getDescription = function(database,organism){
  if(database=="GO"){
    if(organism=="human") descr = GO_human_description else descr = GO_mouse_description
  } else if(database=="Reactome"){
    if(organism=="human") descr = Reactome_human_description else descr = Reactome_mouse_description
  } else if(database=="KEGG"){
    if(organism=="human") descr = KEGG_human_description else descr = KEGG_mouse_description
  } else{
    if(organism=="human") descr = MSigDB_human_description else descr = MSigDB_mouse_description
  }
  
  return(descr)
}


# analysis of geneset enrichment using Fisher's Exact Test
# 
# diff_mat: DESeqDataSet object returned by testDifferentialExpression
# max_gs_size: maximal size of a geneset to be considered in the analysis (300 by default)
# database: name of the database to use (one of either GO (default), Reactome, KEGG or MSigDB)
# evidence: vector of GO evidence codes to restrict database entries onto; with NULL (the default), all evidence codes are used
# ontology: GO ontology to use (one of either cellular component, biological process or molecular function (default))
# organism: either human or mouse
# 
# value: a dataframe with metadata information and BH-adjusted p-Values for all genesets in the selected database
fisherEnrichment = function(diff_mat,max_gs_size=300,database=c("GO","Reactome","KEGG","MSigDB"),evidence=NULL,ontology=c("molecular_function","biological_process","cellular_compartment"),organism=c("human","mouse")){
  # get variables, database and geneset descriptors
  database = match.arg(database)
  organism = match.arg(organism)
  ontology = match.arg(ontology)
  
  db = getDB(database,organism,evidence,ontology)
  descr = getDescription(database,organism)
  
  # mutually restrict expression data and database to the same genes
  genes = intersect(unique(unlist(db)),rownames(diff_mat))
  diff_mat = diff_mat[genes,]
  db = lapply(db,function(x){x[x%in%genes]})
  n = length(genes) # the size of the unsiverse is defined as the total number of measured genes with at least one corresponding annotation
  
  # seperately test geneset enrichment for upregulated and downregulated genes
  # upregulation
  diff_up = diff_mat[which(diff_mat$log2FoldChange>0 & diff_mat$padj<=0.1),]
  de_up = nrow(diff_up)
  up = fisher(db,diff_up,de_up,n,"up",max_gs_size,descr)
  
  # downregulation
  diff_down = diff_mat[which(diff_mat$log2FoldChange<0 & diff_mat$padj<=0.1),]
  de_down = nrow(diff_down)
  down = fisher(db,diff_down,de_down,n,"down",max_gs_size,descr)
  
  # combine the results of both analyses
  df = rbind(up,down)
  df = df[order(df$FDR),]
  
  return(df)
}


# auxiliary function to perform Fisher's Exact Test for geneset enrichment for either up- or downregulated genes; is only called
# internally by fisherEnrichment()
# 
# db: name of the database to use (one of either GO (default), Reactome, KEGG or MSigDB)
# diff_mat: DESeqDataSet object returned by testDifferentialExpression
# de: number of differentially expressed genes
# n: size of the universe
# direction: either "up" (overrepresented genesets) or "down" (underrepresented genesets)
# max_gs_size: maximal size of a geneset to be considered in the analysis (300 by default)
# descr: dataframe containing geneset descriptors belonging to db (as specified before)
# 
# value: a dataframe with metadata information and BH-adjusted p-Values for all up-/downregulated genesets in the selected database
fisher = function(db,diff_mat,de,n,direction=c("up","down"),max_gs_size,descr){
  df = data.frame(Term=character(),Description=character(),Overlap=integer(),Termsize=integer(),Universe=integer(),OR=double(),ConfIntLower=double(),ConfIntUpper=double(),pValue=double(),Direction=character(),Genes=character(),stringsAsFactors=F)
  for(i in names(db)){
    geneset = db[[i]]
    gs = length(geneset)
    if(gs<=max_gs_size & gs>=1){
      de_gs = length(intersect(rownames(diff_mat),geneset))
      if(is.na(de_gs)) de_gs=0
      f = fisher.test(matrix(c(de_gs,de-de_gs,gs-de_gs,n-de-gs+de_gs),2,2,byrow=T))
      df = rbind(df,data.frame(Term=i,Description=descr[descr[,1]==i,2],Overlap=de_gs,Termsize=gs,Universe=n,OR=f$estimate,ConfIntLower=f$conf.int[1],ConfIntUpper=f$conf.int[2],pValue=f$p.value,Direction=direction,Genes=paste(geneset,collapse=',')))
    }
  }
  
  rownames(df) = df$Term
  df$FDR = p.adjust(df$pValue,method="fdr")
  df=df[order(df$FDR),c(1:9,12,10:11)]
  
  return(df)
}


# a wrapper function to perform GSEA using the corresponding java tool
# 
# data: internal expression data object
# condition: String giving the identifier of the "treatment" group
# ref: String giving the identifier of the reference group (default is "control")
# max_gs_size: maximal size of a geneset to be considered in the analysis (300 by default)
# database: name of the database to use (one of either GO (default), Reactome, KEGG or MSigDB)
# evidence: vector of GO evidence codes to restrict database entries onto; with NULL (the default), all evidence codes are used
# ontology: GO ontology to use (one of either cellular component, biological process or molecular function (default))
# organism: either human or mouse
# preranked: logical indicating if preranked GSEA should be performed (FALSE by default)
# diff_mat: DESeqDataSet object returned by testDifferentialExpression (has to be provided if preranked os TRUE)
# 
# value: a dataframe with metadata information and BH-adjusted p-Values for all genesets in the selected database
gseaEnrichment = function(data,condition,ref="control",max_gs_size=300,database=c("GO","Reactome","KEGG","MSigDB"),evidence=NULL,ontology=c("molecular_function","biological_process","cellular_compartment"),organism=c("human","mouse"),preranked=F,diff_mat=NULL){
  # get variables
  if(preranked==T & is.null(diff_mat)) stop("A differential expression matrix has to be provided in preranked mode.")
  database = match.arg(database)
  organism = match.arg(organism)
  ontology = match.arg(ontology)
  design = data$design
  
  # get database information and geneset description (subset of GO terms by evidence code, if specified)
  db = getDB(database,organism,evidence,ontology)
  db = lapply(db,function(x){x[x%in%rownames(data$rlog)]}) # restrict database to genes for which expression values are available
  descr = getDescription(database,organism)
  
  # produce expression matrix for GSEA analysis (according to the gct format specification)
  gsea_df = cbind(Name=rownames(data$rlog),Description=NA,data$rlog)
  write("#1.2","expression.gct")
  write(paste(dim(data$rlog),collapse="\t"),"expression.gct",append=T)
  write.table(gsea_df,"expression.gct",append=T,row.names=F,col.names=T,quote=F,sep="\t")
  
  # produce geneset file for GSEA analysis (according to the gmt format specification)
  for(i in names(db)){
    geneset = db[[i]]
    if(length(geneset)<=max_gs_size & length(geneset)>=1){
      write(paste(i,"NA",paste(geneset,collapse="\t"),sep="\t"),"genesets.gmt",append=T)
    }
  }
  
  # produce phenotype information file for GSEA analysis (according to the cls format specification)
  write(paste(ncol(data$rlog),length(unique(data$Targets[,design])),1,sep=" "),"phenotypes.cls")
  write(paste("#",paste(unique(data$Targets[,design]),collapse=" "),sep=" "),"phenotypes.cls",append=T)
  write(paste(data$Targets[,design],collapse=" "),"phenotypes.cls",append=T)
  
  # perform GSEA analysis using the GSEA java tool and save results in a new folder (located in the wd), labeled with timestamp
  path=getwd()
  if(preranked==F){
    shell(paste0("java -cp \"",gsea_path,"\" -Xmx4G xtools.gsea.Gsea -res ",path,"/expression.gct -cls phenotypes.cls#",condition,"_versus_",ref," -gmx ",path,"/genesets.gmt -collapse false -mode Max_probe -norm None -nperm 1000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis -metric tTest -sort real -order descending -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 15 -zip_report false -out ",path," -gui false"))
    
    name = paste0("GSEA_",gsub(" |:","_",substr(date(),5,19),perl=T))
    shell(paste0("move my_analysis.Gsea.* ",name))
    inpath_up = paste0(name,"/",grep(paste0("gsea_report_for_",condition,".*.xls"),list.files(paste0(path,"/",name)),value=T))
    inpath_down = paste0(name,"/",grep(paste0("gsea_report_for_",ref,".*.xls"),list.files(paste0(path,"/",name)),value=T))
  } else{ # preranked GSEA
    # produce preranked geneset list (according to the rnk format specification)
    # ranking score = foldchange * log10(1/p-value)
    prList=data.frame(Gene=rownames(diff_mat),Rankscore=diff_mat$log2FoldChange*log10(1/diff_mat$pvalue))
    prList=prList[!is.na(prList[,2]),]
    write.table(prList,"preranked.rnk",col.names=F,row.names=F,quote=F,sep="\t")
    shell(paste0("java -cp \"",gsea_path,"\" -Xmx4G xtools.gsea.GseaPreranked -gmx ",path,"/genesets.gmt -norm None -nperm 1000 -rnk ",path,"/preranked.rnk -scoring_scheme weighted -rpt_label my_analysis -create_svgs false -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false -out ",path," -gui false"))
    
    name = paste0("GSEA_preranked_",gsub(" |:","_",substr(date(),5,19),perl=T))
    shell(paste0("move my_analysis.GseaPreranked* ",name))
    inpath_up = paste0(name,"/",grep(paste0("gsea_report_for_na_pos.*.xls"),list.files(paste0(path,"/",name)),value=T))
    inpath_down = paste0(name,"/",grep(paste0("gsea_report_for_na_neg.*.xls"),list.files(paste0(path,"/",name)),value=T))
  }
  
  # read GSEA output and produce dataframe
  infile_up <- read.csv(inpath_up,header=T,sep="\t",colClasses=c("character","character","character","integer","double","double","double","double","double","integer","character",NA),stringsAsFactors=F)
  infile_up = cbind(infile_up,Direction="up")
  infile_down <- read.csv(inpath_down,header=T,sep="\t",colClasses=c("character","character","character","integer","double","double","double","double","double","integer","character",NA),stringsAsFactors=F)
  infile_down = cbind(infile_down,Direction="down")
  infile = rbind(infile_up,infile_down)
  infile$description = apply(infile,1,function(x){return(as.character(descr[descr[,1]==x[1],2]))})
  infile = infile[,c(1,14,4,5,6,7,8,9,10,13)]
  colnames(infile) = c("Term","Description","Termsize","ES","NES","pValue","FDR","FWER","RankAtMax","Direction")
  rownames(infile) = infile$Term
  
  shell("del *.gct *.gmt *.cls *.rnk")
  return(infile[order(infile$FDR),])
}


# analysis of geneset enrichment using the Globaltest method described by Goeman et al. (2005)
# 
# data: internal expression data object
# condition: String giving the identifier of the "treatment" group
# ref: String giving the identifier of the reference group (default is "control")
# max_gs_size: maximal size of a geneset to be considered in the analysis (300 by default)
# database: name of the database to use (one of either GO (default), Reactome, KEGG or MSigDB)
# evidence: vector of GO evidence codes to restrict database entries onto; with NULL (the default), all evidence codes are used
# ontology: GO ontology to use (one of either cellular component, biological process or molecular function (default))
# organism: either human or mouse
# 
# value: a dataframe with metadata information and BH-adjusted p-Values for all genesets in the selected database
globalTest = function(data,condition,ref="control",max_gs_size=300,database=c("GO","Reactome","KEGG","MSigDB"),evidence=NULL,ontology=c("molecular_function","biological_process","cellular_compartment"),organism=c("human","mouse")){
  database = match.arg(database)
  organism = match.arg(organism)
  ontology = match.arg(ontology)
  design = data$design
  
  # select database (subset of GO terms by evidence code, if specified)
  db = getDB(database,organism,evidence,ontology)
  db = lapply(db,function(x){x[x%in%rownames(data$rlog)]})
  descr = getDescription(database,organism)
  
  tab = data$Targets[data$Targets[,design] %in% c(ref,condition),]
  outcome = as.numeric(as.factor(tab[,design]))-1
  mat = t(data$rlog[,rownames(tab)]) # for globaltest, the count matrix has to be be transposed
  
  db = db[lengths(db)<=max_gs_size & lengths(db)>=1]
  res = result(p.adjust(gt(outcome,mat,model="logistic",subsets=db,directional=T,standardize=T),"BH")) # all genesets are tested with the same function
  df = data.frame(Term=names(db),Description=sapply(names(db),function(x){return(descr[descr[,1]==x,2])}),Termsize=lengths(db),Statistic=res$Statistic,Expected=res$Expected,Std.dev=res$Std.dev,pValue=res[,1],FDR=res$BH,Genes=unlist(lapply(db,function(x){return(paste(x,collapse=","))})))
  
  rownames(df) = df$Term
  df=df[order(df$FDR),]
  return(df)
}


