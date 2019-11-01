

library(data.table)



# merge both tables
path = '/home/l/loipfingers/Desktop/self_neap/HMDC'
setwd(path)

geneTable = fread("output_disConnection.txt")
geneTable$V1 <- NULL
omimTable = fread("athIDs.txt")

# return name of ids
checkAthDO <- function(disId) {
  r = omimTable[omimTable$DO_ID==disId]
  if(disId=='' | nrow(r)==0) {
    return('')
  }
  else{
    return(paste0(r$MalaCard_Name,' [',r$MalaCard_ID,']'))
  }
}

# return name of ids
checkAthOMIM <- function(disId) {
  if(disId=='') {
    return('')
  }
  disId = gsub('OMIM:','',disId)
  disId = strsplit(disId,';')[[1]]
  
  r = omimTable[omimTable$OMIM_ID %in% disId]
  if(nrow(r)==0) {
    return('')
  }
  else{
      r = unique(r)
      mala = paste0(r$MalaCard_Name,' [',r$MalaCard_ID,']')
      mala = paste(mala, collapse=';')
      return(mala)
  }
}

mergeCols <- function(col1, col2) {
  c = gsub(col1,'',col2)
  if(c!='') {
    return(paste0(col1,c))
  }
  else {
    return('')
  }
}


geneTable$athSpecificMalaCard_DO <- sapply(geneTable$DO_ID, FUN=function(x) checkAthDO(x))
geneTable$athSpecificMalaCard_OMIM <- sapply(geneTable$OMIM_ID, FUN=function(x) checkAthOMIM(x))
geneTable$athSpecificMalaCard <- ifelse(geneTable$athSpecificMalaCard_DO!=geneTable$athSpecificMalaCard_OMIM, paste0(geneTable$athSpecificMalaCard_DO,';', geneTable$athSpecificMalaCard_OMIM), geneTable$athSpecificMalaCard_DO)
geneTable$athSpecificMalaCard_DO <- NULL
geneTable$athSpecificMalaCard_OMIM <- NULL

write.table(geneTable, 'output_athSpecifc_final.txt', sep='\t', quote=FALSE)


