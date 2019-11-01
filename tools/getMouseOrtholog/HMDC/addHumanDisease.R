library(data.table)


g <- fread('/home/stefan/bioinformatics/neap/getMouseOrtholog/HMDC/MGI_diseaseConnection.txt')
old <- fread('/home/stefan/bioinformatics/neap/getMouseOrtholog/HMDC/output_athSpecifc_final_OLD.txt')
old$V1 <-NULL



old$diseaseName2 <- paste(old$diseaseName," [",old$DO_ID,old$OMIM_ID,"]")

old <- within(old, {
  diseaseConnection <- ifelse(diseaseName != "", paste0(old$diseaseName," [",old$DO_ID,";",old$OMIM_ID,"]"), "")
  #a2_t <- ifelse(diseaseName != "" & DO_ID!="", paste(old$diseaseName," [",old$DO_ID,";",old$OMIM_ID,"]"))
})

getDis <- function(ge,org) {
  u <- g[g$Common_Organism_Name==org & g$Symbol==ge]
  u <- unique(u)
  out <- ""
  if(nrow(u)==0) { return("")}
  else {
    out <- paste0(u$DO_disease_name," [",u$DO_disease_ID,";",gsub("\\|",";",u$OMIM_IDs),"]")
    return(paste(out,collapse=";")) }
}
  
# getDis("ZZ3","human")
# getDis("GDF5","human")


old$humanDisCon <- sapply(old$geneID_human, getDis, "human")
old$mouseDisCon <- sapply(old$geneID_mouse, getDis, "mouse, laboratory")

out <- old[,c("geneID_human","geneID_mouse","mouseHumanConnection","diseaseConnection","humanDisCon","mouseDisCon","mouseMGI_ID","mouseStrain","mouseMP_ID","MPnames","athSpecificMalaCard")]

out$mouseStrain <- gsub(",",";",out$mouseStrain)
out$mouseMGI_ID <- gsub(" ",";",out$mouseMGI_ID)
out$mouseMP_ID <- gsub(" ",";",out$mouseMP_ID)
out$MPnames <- gsub(",",";",out$MPnames)

out$mouseStrain <- gsub('"',"",out$mouseStrain)
out$mouseMGI_ID <- gsub('"',"",out$mouseMGI_ID)
out$mouseMP_ID <- gsub('"',"",out$mouseMP_ID)
out$MPnames <- gsub('"',"",out$MPnames)
out$diseaseConnection <- gsub('"',"",out$diseaseConnection)
out$humanDisCon <- gsub('"',"",out$humanDisCon)
out$mouseDisCon <- gsub('"',"",out$mouseDisCon)


write.table(out,"/home/stefan/bioinformatics/neap/getMouseOrtholog/HMDC/output_athSpecifc_final.csv", row.names = F,quote = T,sep=",")

# remove all ' signs in datatable afterwards !!!!!








