
library(data.table)

# user input, geneList per file input
#geneFile
#geneList 
#mouseStrain
#strictMouseStrain

geneList = c(geneTable$geneID_human[1:40],'Aorls1')
geneList = gsub("\\s+", "", geneList, perl=TRUE)
mouseStrain = 'Aadac, MGI:2152878, MGI:1917115,A3galt2' # list of mouse strains or MGI
strictMouseStrain = TRUE

handleGeneFile <- function(path){
  f = readChar(path, file.info(path)$size)
  return(f)
}
handleGeneFile("geneListTest.txt")

# handle input
geneList = gsub("\\s+", "", geneList, perl=TRUE)
mouseStrain = gsub("\\s+", "", mouseStrain, perl=TRUE)


# merge both tables
path = '/home/l/loipfingers/Desktop/self_neap/HMDC'
setwd(path)

# load table
geneTable <- fread("output_athSpecifc_final.txt")
geneTable$V1 <- NULL 

#create output table by iterating over gene names
searchTable=data.frame(matrix(ncol = length(geneTable), nrow = 0))
colnames(searchTable) <- colnames(geneTable)

for( ge in geneList) {
  row1 = geneTable[geneTable$geneID_human==ge | geneTable$geneID_mouse == ge]
  if(nrow(row1)==0) {
    outputRow = data.frame(ge,"","","","","","","","","","")  # gene not found
    names(outputRow) <- names(searchTable)
  }
  else {
    outputRow=row1
  }
  searchTable <- rbind(searchTable, outputRow)
}

searchTable <- data.frame(lapply(searchTable, as.character), stringsAsFactors=FALSE)


# get only certain mice strain
mgiStrain = strsplit(mouseStrain,',')[[1]]   # filter by mgi
mgiTable = searchTable[which(searchTable$mouseMGI_ID %in% mgiStrain),]

mouseStrain = gsub("MGI:\\d+", "", mouseStrain, perl=TRUE)  # remove mgi
mouseStrain = strsplit(mouseStrain,',')[[1]]   # filter by mgi
mouseStrain = mouseStrain[mouseStrain!='']
strainTable = copy(searchTable)
delRows = c()
if (mouseStrain !='') {   # filter by name
  if (strictMouseStrain==TRUE) {
    for(r in 1:nrow(strainTable)) {
      s = strsplit(strainTable[r,]$mouseStrain,",")[[1]]
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
if(mouseStrain !='') {
  finalTable <- unique(rbind(strainTable,mgiTable))
}


# return tables
finalTable   # full search Table
athTable <- searchTable[searchTable$athSpecificMalaCard!=''] # only ath specifc table





