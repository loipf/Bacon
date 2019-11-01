# give in geneList (+ mouse strain) and analyse for human connected genes

# better and faster to pre-create table
#args = commandArgs(trailingOnly=TRUE)
#geneList = args[1]
#mouseStrain = args[2]
#strictMouseStrain = args[3]

# eventuell genelist per file eingeben?
geneList = c('OPA3','A1CF','A4GALT','Gja8','Aorls1') # create for whole genelist later
mouseStrain = ''  # oder keinen
strictMouseStrain = FALSE


library(data.table)
#library(biomaRt)


# load in tables
path = '/home/l/loipfingers/Desktop/self_neap/HMDC'
setwd(path)

# handle special loading
humanMouseTable <- fread("HMD_HumanMousePheno.txt")
colnames(humanMouseTable) <- c(colnames(humanMouseTable)[2:length(humanMouseTable)],"del")
humanMouseTable[[length(humanMouseTable)]] <- NULL

diseaseTable <- fread("MGI_Geno_DiseaseDO.txt")
nonDiseaseTable <- fread("MGI_Geno_NotDiseaseDO.txt")
mgiDiseaseTable <- fread("MGI_diseaseConnection.txt")

sexTable <- fread("others/MGI_Pheno_Sex.txt")
mpNameTable <- fread("others/MP_name.txt", header = FALSE)
strainTable <- fread("MGI_Strain.txt")
strainTable2 <- fread("others/MGI_MousePhenotypicAllele.txt", quote = "")
strainTable3 <- fread("others/MGI_MouseQTLAllele.txt")



# get sex specifity and according pubmed ID, left out
getSexSpecific <- function(mgi,mp) {
  s = sexTable[sexTable$genotype_ID==mgi & sexTable$MP_ID==mp]
  ma = s[s$sex=="M"]$`sex-specific_YesNo`
  fe = s[s$sex=="F"]$`sex-specific_YesNo`
  if(length(ma)==0) {
    ma='-'
  }
  if(length(fe)==0){
    fe='-'
  }
  outS = paste0(ma,fe,collapse = '')
  outP = paste0(s$PubMed_MGI, collapse=',')
  print(s)
  return(c(outS,outP))
}

# get name of mouse phenotypes MP
getMPname <- function(mp) {
  singleMps = strsplit(mp," ") [[1]]
  ph = mpNameTable[mpNameTable$V1 %in% singleMps,]$V2
  return(paste0(ph, collapse = ","))
}

# get genetic background name
getStrainName <- function(mgi) {
  singleMgis = strsplit(mgi," ") [[1]]
  ph1 = strainTable[strainTable$MGI_ID %in% singleMgis,]$strainName
  ph2 = strainTable2[strainTable2$MGI_ID %in% singleMgis | strainTable2$MGI_ID2 %in% singleMgis,]$allele_symbol
  ph3 = strainTable3[strainTable3$MGI_ID %in% singleMgis | strainTable3$MGI_ID2 %in% singleMgis,]$allele_name
  ph = unique(c(ph1,ph2,ph3))
  return(paste0(ph, collapse = ","))
}

getDiseaseCon <- function(homoGeneId) {
  ph = mgiDiseaseTable[mgiDiseaseTable$HomoloGene_ID == homoGeneId,]
  disOverlap = ph[ph$DO_disease_ID == ph[duplicated(ph$DO_disease_ID)]$DO_disease_ID]
  t = disOverlap$NCBI_Taxon_ID
  if (9606 %in% t & 10090 %in% t) {
    if(length(t)>3) {
      print(paste0("problem: multipleDOs ",homoGeneId))
    }
    return(disOverlap)
  }
}


geneList <- humanMouseTable$human_marker # create table for all marker genes found

### iterate over gene list
# add mgi http://www.informatics.jax.org/allele/genoview/MGI:3706990
# paste0(names(mainTable),collapse=", ")
mainTable <- data.frame(geneID_human=character(), geneID_mouse=character(),  humanMouseConnection=character(), diseaseName=character(), DO_ID=character(), OMIM_ID=character(), mouseStrain=character(), mouseMGI_ID=character(), mouseMP_ID=character(), MPnames=character())
for( ge in geneList) {
  row1 = humanMouseTable[humanMouseTable$human_marker==ge | humanMouseTable$mouse_marker == ge]
  if(nrow(row1)==0) {
      outputRow = data.frame(ge,"","","","","","","","","")  # change according #col
      names(outputRow) <- names(mainTable)
  }
  else {
    mgi <- getStrainName(row1$MGI_maker_ID) 
    mp <- getMPname(row1$MP_ID)
    dis <- getDiseaseCon(row1$HomoloGene_ID)
    
    if(!is.null(dis)) {
      outputRow = data.frame(row1$human_marker, row1$mouse_marker, 'y', dis$DO_disease_name[1], dis$DO_disease_ID[1], dis$OMIM_IDs[1], mgi,row1$MGI_maker_ID, row1$MP_ID, mp )
      names(outputRow) <- names(mainTable)
    }
    else {
      outputRow = data.frame(row1$human_marker, row1$mouse_marker, 'n',"","","",mgi,row1$MGI_maker_ID,row1$MP_ID,mp) # and many more
      names(outputRow) <- names(mainTable)
      }
  }

  mainTable <- rbind(mainTable, outputRow)
}

mainTable <- data.frame(lapply(mainTable, as.character), stringsAsFactors=FALSE)


# get only certain mice strain
delRows = c()
if (mouseStrain !='') {
  if (strictMouseStrain==TRUE) {
    for(r in 1:nrow(mainTable)) {
      s = strsplit(mainTable[r,]$mouseStrain,",")[[1]]
      if(!mouseStrain %in% s) {
        delRows <- c(delRows, r)
      }
    }
  }
  else {
    for(r in 1:nrow(mainTable)) {
      s = mainTable[r,]$mouseStrain
      if(!grepl(mouseStrain,s)) {
        delRows <- c(delRows, r)
      }
    }
  }
  
  mainTable <- mainTable[-delRows,]
}


write.table(mainTable, "outputMainTable.txt", quote=FALSE, sep = "\t")
#print(mainTable)







