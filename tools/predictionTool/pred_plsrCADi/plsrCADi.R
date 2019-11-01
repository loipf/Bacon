

library(pls)

# 0=not possible, 1=possible, 2=new model needed
cadi_ownPossible <- 0
cadi_paperPossible <- 0
cadi_currOwnModel <- NULL
cadi_currPaperModel <- NULL


scriptPath <- ("/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/predictionTool/")
scriptPath <- ("/Users/stefanloipfinger/Documents/GitHub/Bacon/tools/predictionTool/")




# initliaze all models
plsr_initialize <- function(){
  plsr_own <- readRDS(paste0(scriptPath,"pred_plsrCADi/plsr_ownAnalysis.rds"))
  plsr_paper <- readRDS(paste0(scriptPath,"pred_plsrCADi/plsr_paperAnalysis.rds"))
  
  geneList_own <<- as.character(read.table(paste0(scriptPath, "geneLists/genes_plsrCADi_own.txt"))$V1)
  geneList_paper <<- as.character(read.table(paste0(scriptPath, "geneLists/genes_plsrCADi_paper.txt"))$V1)
}


plsr_compareGeneLists <- function(predDataGenes) {
  inters_own <<- intersect(predDataGenes, geneList_own)
  inters_paper <<- intersect(predDataGenes, geneList_paper)
  
  compareOut <- paste0("found ", length(inters_own),"/",length(geneList_own)," of own and ", length(inters_paper),"/",length(geneList_paper)," of paper genes")
  if( length(inters_own)<(0.2*length(geneList_own)) ) {
    cadi_ownPossible <<- 0
    compareOut <- paste0(compareOut, " \n -> prediction for own analysis not possible")
  } else {
    if( length(inters_own)==length(geneList_own) ){
      cadi_ownPossible <<- 1
    }
    else {
      cadi_ownPossible <<- 2
      compareOut <- paste0(compareOut, " \n -> made up  new own training model with less genes -> LOWER ACCURACY !!!")
    }
  }
  
  if( length(inters_paper)<(0.2*length(geneList_paper)) ) {
    cadi_paperPossible <<- 0
    compareOut <- paste0(compareOut, " \n -> prediction for own analysis not possible")
  } else {
    if( length(inters_paper)==length(geneList_paper) ){
      cadi_paperPossible <<- 1
    }
    else {
      cadi_paperPossible <<- 2
      compareOut <- paste0(compareOut, " \n -> made up new paper training model with less genes -> LOWER ACCURACY !!!")
    }
  }
  plsr_calculateNewModels()
  return(compareOut)
}


plsr_calculateNewModels <- function() {
  if(cadi_ownPossible==1) {
    cadi_currOwnModel <<- plsr_own
  }
    else{ if(cadi_ownPossible==2) {
      m <- read.table(paste0(scriptPath,"pred_plsrCADi/exprsMatrix_own.txt"), sep = "\t")
      trainD <- m[colnames(m) %in% c(inters_own,"class")]
      cadi_currOwnModel <<- plsr(class~ ., ncomp =10, scale = TRUE, data=trainD)
        }
    }
  if(cadi_paperPossible==1) {
    cadi_currPaperModel <<- plsr_paper
  }
    else{ if(cadi_paperPossible==2) {
      m <- read.table(paste0(scriptPath,"pred_plsrCADi/exprsMatrix_paper.txt"), sep = "\t")
      trainD <- m[colnames(m) %in% c(inters_paper,"class")]
      cadi_currPaperModel <<- plsr(class~ ., ncomp =10, scale = TRUE, data=trainD)
    }
    }

}

plsr_plot_perf <- function(mod, predValues){
  #currOwnModel 
  predplot(mod, ncomp = 2, asp = 1, line = TRUE, xlim=c(0,100), ylim = c(0,100), main="prediction of CADi")
  #plot(currOwnModel, ncomp =2, xlim=c(-10,100), ylim = c(-10,100), xlab="measured", ylab = "predicted", main="performance CADi of own geneset", line=TRUE)
  #abline(lm(currOwnModel$Yscores ~ currOwnModel$Xscores), col="red")
  
  if(length(predValues) == 2) {
    points( predValues[,2], predValues[,1],pch=20, col="red")
  }

	#print(mod)
	#print(predValues)
  #points(currOwnModel$Ymeans, pch=20, col="red")
  #points(currOwnModel$Yscores,currOwnModel$fitted.values, pch=20, col="red")  ###
  #points(currOwnModel$fitted.values, pch=20, col="red")
  #points(currOwnModel$model, pch=20, col="green")
}

plsr_plot_rmsep <- function(mod) {
  plot(RMSEP(mod), legendpos = "topright", ylim=c(0,35))
}



plsr_predTable_own <- function(mod, predM) {
  testData <- predM[colnames(predM) %in% inters_own]
  predict(mod, ncomp = 2, newdata = testData)
}

plsr_predTable_paper <- function(mod, predM) {
  testData <- predM[colnames(predM) %in% inters_paper]
  predict(mod, ncomp = 2, newdata = testData)
}







