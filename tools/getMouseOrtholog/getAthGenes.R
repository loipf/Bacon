
library(data.table)

# user input, geneList per file input
args = commandArgs(trailingOnly=TRUE)
geneList <- args[1]
mouseStrain <- args[2]
strictMouseStrain <- args[3]

geneList <- c("PAX1,TZSD,ZZ3, UI,P53,ALMS1,A2m,AAR2,ABCC10,A3GALT2,ABCC6")
mouseStrain <- 'all'
mouseStrain <- "MGI:2685279"
strictMouseStrain <- FALSE


#path = '/home/l/loipfingers/Desktop/self_neap/HMDC'
path = '/home/stefan/bioinformatics/neap/getMouseOrtholog'
pathGetGeneNames = '/home/stefan/bioinformatics/neap/getUniqueGeneName/'

# handle input
geneList = gsub("\\s+", "", geneList, perl=TRUE)
mouseStrain = gsub("\\s+", "", mouseStrain, perl=TRUE)

# get HDNC gene names from synonymes
system(paste0("Rscript ",pathGetGeneNames,"getGeneNames.R"," ",geneList," ",path,"/synonymTable.txt"))
hdnc <- read.table(paste0(path,"/synonymTable.txt"), sep="\t", header=TRUE,stringsAsFactors = FALSE)

geneList <- strsplit(geneList,",")[[1]]
#print(geneList)
#print(mouseStrain)
#print(strictMouseStrain)


# load table
geneTable <- read.csv(paste0(path,"/output_athSpecifc_final.csv"), sep=",", stringsAsFactors = F)
#geneTable$V1 <- NULL 
#geneTable[is.na(geneTable)] <- ""

#create output table by iterating over gene names
searchTable=data.frame(matrix(ncol = (length(geneTable)+1), nrow = 0))
colnames(searchTable) <- c("input",colnames(geneTable))

for( ge in geneList) {
  #print(ge)
  entryHdnc <- hdnc[hdnc$input==ge,]
  row1 = geneTable[geneTable$geneID_human== entryHdnc$humanGene[1] | geneTable$geneID_mouse == entryHdnc$mouseGene[1],]
  if(nrow(row1)==0) {
    outputRow = data.frame(ge,entryHdnc$humanGene[1], entryHdnc$mouseGene[1],"","","","","","","","","")  # gene not found
    names(outputRow) <- names(searchTable)
  }
  else {
    outputRow=data.frame(ge,row1)
    #print(outputRow)
    names(outputRow) <- names(searchTable)
    
  }
  
  searchTable <- rbind(searchTable, outputRow)
}

searchTable <- data.frame(lapply(searchTable, as.character), stringsAsFactors=FALSE)

mouseStrain <- "MGI:2685279"
# get only certain mice strain
mgiStrain = strsplit(mouseStrain,',')[[1]]   # filter by mgi
mgiTable = searchTable[which(searchTable$mouseMGI_ID %in% mgiStrain),]

mouseStrain = gsub("MGI:\\d+", "", mouseStrain, perl=TRUE)  # remove mgi
mouseStrain = strsplit(mouseStrain,',')[[1]]   # filter by mgi
mouseStrain = mouseStrain[mouseStrain!='']
strainTable = copy(searchTable)
delRows = c()
if (!identical(mouseStrain, character(0)) | !identical(mouseStrain,"all")) {   # filter by name
  if (strictMouseStrain=="TRUE") {
    for(r in 1:nrow(strainTable)) {
      s = strsplit(strainTable[r,]$mouseStrain,";")[[1]]
      deleleRow = TRUE
      for(p in mouseStrain) {
        if(p %in% s)  {
          deleleRow = FALSE
        }
      }
      if(deleleRow) {
        delRows <- c(delRows,r)
      }
    }
  }
  else {
    for(r in 1:nrow(strainTable)) {
      s = strainTable[r,]$mouseStrain
      deleleRow = TRUE
      for(p in mouseStrain) {
        if(grepl(p,s) ) {
          deleleRow = FALSE
        }
      }
      if(deleleRow) {
        delRows <- c(delRows,r)
      }
    }
  }
    strainTable <- strainTable[-delRows,]
}


finalTable <- searchTable
if(!identical(mouseStrain,"all")) {
  finalTable <- unique(rbind(strainTable,mgiTable))
}


# return tables
#finalTable   # full search Table
athTableMala <- finalTable[finalTable$athSpecificMalaCard!='',] # only ath specifc table malacard
athTableOwn <- finalTable[finalTable$geneID_human %in% c("A1BG","A1CF","A2M","A3GALT2","A4GALT","ABCC6"),]  # REPLACE WITH OWN GENESET




write.table(finalTable, paste0(path,"/fullAthGenes.csv"), sep=",", row.names=FALSE, quote=FALSE)

if(nrow(athTableMala)!=0) {
	write.table(athTableMala, paste0(path,"/malaAthGenes.csv"), sep=",", row.names=FALSE, quote=FALSE)
}

if(nrow(athTableOwn)!=0) {
	write.table(athTableOwn, paste0(path,"/ownAthGenes.csv"), sep=",", row.names=FALSE, quote=FALSE)
}




