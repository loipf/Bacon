library(shiny)
library(GEOquery)
library(DT)


### TO DO ###
# specify temp dir getGEO files get saved in
# line 224 change !!
# export files in folder not in my home folder !!! (change this if you want to use it, also EXAMPLE)
# server port specification
# outputFILE start with GENE then list without row names

# execute before running app
# setwd('/home/stefan/bioinformatics/neap/getGEOdata/')
# setwd('D:/studium/8_semester/NEAP/getGEOdata/')
# shiny::runApp('/home/stefan/bioinformatics/neap/predictionTool', port = 7817) # dont know


path <- "/home/stefan/bioinformatics/neap/"
saveDirPath <- paste0(path,"datasets/")

#path <- "/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/"
#saveDirPath <- paste0(path,"finishedDatasets/Microarray_processedData/")

# global variables
showedLoaded <- FALSE
fileLoaded <- FALSE
gsetFile <- NULL 
gset <- NULL 
importedFile <- FALSE



server <- function(input, output) {
  options(shiny.maxRequestSize=1000*1024^2)  # upload up to 1GB of data

  ###### load data #####
  
  # downloads GEO file
  loadGEOfile <- eventReactive(input$buttonGetGEO, {
  #loadGEOfile <- function() {
    print("pressed")
    print(showedLoaded)
    # check if already loaded
    if (dir.exists(paste0(saveDirPath,input$GEOacc)) & (showedLoaded==FALSE)) {
      output$GEOfileLoaded <- renderUI({ HTML(paste0("dataset ",input$GEOacc," already loaded in!<br/>load again if you want to reload it")) })
    }
    else { showedLoaded <<-TRUE }
    
    if(showedLoaded==TRUE) {
      fileLoaded <<- FALSE
  
      tryCatch( {
        withProgress(message = "download GEO file", value = 0, {
        incProgress(1/10) # loading bar
        gsetFile <<- getGEO(input$GEOacc,GSEMatrix =TRUE, AnnotGPL=TRUE) # save somewhere or delete later!!!
        incProgress(9/10) # loading bar
        fileLoaded <<- TRUE
        incProgress(10/10) # loading bar
        })
        return(TRUE)
      },
        error = function(e) {
          fileLoaded <<- FALSE 
        return(FALSE)}
      )}
    return(FALSE)
    print(showedLoaded)
   } )

  # change download text
  observeEvent(input$buttonGetGEO, {
    loadBoolean <- loadGEOfile()
    if (loadBoolean) {
      outText <- HTML(paste0("... sucessfully downloaded, choose platform"))
      output$GEOfileLoaded <- renderUI({ outText })
       }
      else {
        if(!loadBoolean & showedLoaded ==FALSE) { showedLoaded <<- TRUE  }
        else {  
          outText<- "error loading file"
          output$GEOfileLoaded <- renderUI({ outText })}
      }
    
  })
  

  
  # show platforms
  loadGEOplatforms <- function(){
      platformList <- c()
      for (i in 1:length(gsetFile)) {
        platformList <- c(platformList, annotation(gsetFile[[i]]))
      }
      names(platformList)<-sapply(platformList,paste)
      for(i in 1:length(platformList)) { # bad way but too late to think about apply
        platformList[i] <- i
      }
     return(platformList)
  }

  # if file loaded show platform + button
  observeEvent(input$buttonGetGEO, {
    if(fileLoaded) {
    output$GEOplatform <- renderUI({ selectInput("GEOplatform","platform",choices = loadGEOplatforms())  })
    output$buttonGetPlatform <- renderUI({  actionButton("buttonGetPlatform", "select") })
    }
  })
  
  
  # geo title loading with perl script
  loadGEOtitle <- eventReactive(input$buttonGetGEO, {
    print("works")
    t <- "you really should install Perl on your PC !!"
    try( {  
      cmd <- paste("perl", "getGEOname.pl",input$GEOacc)
      t <- system(cmd, intern=TRUE)
    } )

      if(t != 'no match') {
        urlT <- a(t, href=paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",input$GEOacc),target="_blank")
        return(urlT)
      }
      else {
        return(t)
      }
  })

  output$GEOtitle <- renderUI(loadGEOtitle())
  
  
  
  
  ###### choose appropriate columns ######
  output$GEOdataInfo <- eventReactive(input$buttonGetPlatform, {
    gset <<- gsetFile[[as.integer(input$GEOplatform)]]
    re <- paste0("expression set data has ",dim(gset)[2]," samples with ",dim(gset)[1]," features")
    return(re)
  } )

  # print matrices 
  getGsetPhenoTable <- eventReactive(input$buttonGetPlatform, {
    g <- gset[1:50,] # only print first 50
    g <- pData(phenoData(g))
    phenoTableCol <<- colnames(g)
    g <- sapply(g, substring, 1, 40)
  } )
  
  getGsetFeatureTable <- eventReactive(input$buttonGetPlatform, {
    g <- gset[1:50,] # only print first 50
    g <- pData(featureData(g))
    featTableCol <<- colnames(g)
    g <- sapply(g, substring, 1, 40)
  }  )
  
  # create input selector for GEO files
  printArgumentsGEO <- function(phe,fea) {
    output$chooseSampleId <- renderUI({ selectInput("chooseSampleId","choose sample ID column", choices = c(Choose="",phe), selectize = TRUE)})
    output$chooseGroup <- renderUI({ selectInput("chooseGroup","choose ath/non ath column", choices = c(Choose="",phe), selectize = TRUE)})
    output$chooseGeneSymbol <- renderUI({ selectInput("chooseGeneSymbol","choose gene symbol column", choices = c(Choose="",fea), selectize = TRUE)}) 
    output$chooseOtherColumns <- renderUI({ checkboxGroupInput("chooseOtherColumns", "other important columns", phe[1:6]) })
     }
  
  
  # show second and third wellPanel after input, also create preview tables
  output$platformLoaded <- eventReactive(input$buttonGetPlatform, {
    
    output$showPhenoTable <- DT::renderDataTable(
      DT::datatable(getGsetPhenoTable(), options = list(searching = FALSE,pageLength = 5))
    )
    
    output$showFeatureTable <- DT::renderDataTable(
      DT::datatable(getGsetFeatureTable(), options = list(searching = FALSE, pageLength = 5))
    )

    printArgumentsGEO(colnames(getGsetPhenoTable()),colnames(getGsetFeatureTable()))
    return("y") # to silly for boolean
  })
  outputOptions(output,"platformLoaded", suspendWhenHidden = FALSE)
  
  
  # remove columns in other choices because thats cool
  observeEvent(input$chooseSampleId, {
    output$chooseOtherColumns <- renderUI({ checkboxGroupInput("chooseOtherColumns", "other important columns", phenoTableCol[!phenoTableCol %in% c(input$chooseSampleId, input$chooseGroup)]) })
    #output$chooseGroup <- renderUI({ selectInput("chooseGroup","choose ath/non ath column", choices = phenoTableCol[!phenoTableCol ==input$chooseSampleId])})
  })
  
  observeEvent(input$chooseGroup, {
    output$chooseOtherColumns <- renderUI({ checkboxGroupInput("chooseOtherColumns", "other important columns", phenoTableCol[!phenoTableCol %in% c(input$chooseSampleId, input$chooseGroup)]) })
    #output$chooseSampleId <- renderUI({ selectInput("chooseSampleId","choose sample ID column", choices = phenoTableCol[!phenoTableCol ==input$chooseGroup], selected = input$SampleId)})
  })
  
  
  
  
  
  ###### handle data in right format ######
  importGEOFile <- function() {
    importedFile <<- FALSE
    idProp <- input$chooseSampleId
    groupProp <- input$chooseGroup
    geneProp <- input$chooseGeneSymbol
    otherProp <- input$chooseOtherColumns
    
    if(groupProp!="" & geneProp!="" & idProp!="") {
    tryCatch( {
      withProgress(message = "import GEO file in database", value = 0, {
      # export targets
      incProgress(1/10) # loading bar
      #create directory
      ifelse(!dir.exists(file.path(saveDirPath,input$GEOacc)), dir.create(file.path(saveDirPath,input$GEOacc)), FALSE)
        
      phenoSet <- pData(phenoData(gset))[c(idProp,groupProp,otherProp)]
      colnames(phenoSet) <- c("sampleId","group",colnames(phenoSet)[-c(1,2)]) # change names here
      write.table(phenoSet, paste0(saveDirPath,input$GEOacc,"/",input$GEOacc,"_targets.txt"), sep="\t",quote=FALSE,row.names = FALSE) # change export location !!!!
      incProgress(4/10) # loading bar
      
      # export expression matrix
      fd <- featureData(gset)
      e <- exprs(gset)
      geneList <- sapply(rownames(e), function(x) gsub("///.*","",fd[x][[geneProp]]) ) # takes a while
      e$GENE <- geneList
      e <- e[,c("GENE",colnames(e))] # not tested yet but should work
      incProgress(8/10) # loading bar
      write.table(e,paste0(saveDirPath,input$GEOacc,"/",input$GEOacc,"_matrix.txt"), sep="\t",quote=FALSE, row.names = FALSE)
      incProgress(10/10) # loading bar
      importedFile <<- TRUE
      })
      return(TRUE)
    },
    error = function(e) { 
      importedFile <<- FALSE
      return(FALSE)}
    )
      
    }else { return(FALSE) }
  }
  
  
  observeEvent(input$buttonPlatformFinished, {
      importGEOFile()
      output$platformFinsishedStatus <- renderUI(if(importedFile) { return("... sucessfully imported file")} else {"error importing file"})

      })
  

  
  
  ###### load own expression set and targets ######
  
  # example load
  output$exampleTable_exprs <- DT::renderDataTable(DT::datatable(
    read.table(paste0(path,"/getGEOdata/example_matrix.txt"), sep="\t", row.names = NULL, header=TRUE),   # !!!!!!CHANGE !!!!!!!!!!
    options = list(searching = FALSE,pageLength = 11,lengthChange = FALSE, autoWidth = TRUE,
                   columnDefs = list(list(width = '100px', targets = "_all")))
  ))
  
  output$exampleTable_targets <- DT::renderDataTable(DT::datatable(
    read.table(paste0(path,"/getGEOdata/example_targets.txt"), sep="\t", row.names = NULL,header = TRUE),
    options = list(searching = FALSE,pageLength = 11,lengthChange = FALSE, autoWidth = FALSE,
                   columnDefs = list(list(width = '100px', targets = "_all")))
  ))
  
  
  importExprsFile <- function() {
    req(input$matrixFileUp)
    req(input$targetsFileUp)
    
    tryCatch(
      { withProgress(message = "import files in database", value = 0, {
          fileM <- read.csv(input$matrixFileUp$datapath, sep = input$sepFileUp)
          incProgress(3/10) # loading bar
          fileT <- read.csv(input$targetsFileUp$datapath, sep = input$sepFileUp)  # maybe handle double id rows
          incProgress(6/10) # loading bar
          
          # create directory
          folderName <- basename(input$matrixFileUp$datapath)
          folderName <- gsub("_*matrix","",folderName)
          ifelse(!dir.exists(file.path(saveDirPath,folderName)), dir.create(file.path(saveDirPath,folderName)), FALSE)
          
          write.table(fileM,paste0(saveDirPath,folderName,"/",folderName,"_matrix.txt"), sep="\t",quote=FALSE)
          incProgress(8/10) # loading bar
          write.table(fileT,paste0(saveDirPath,folderName,"/",folderName,"_targets.txt"), sep="\t",quote=FALSE)
          incProgress(10/10) # loading bar
       })
        return(TRUE)
      },
      error = function(e) { 
        return(FALSE) }
    
    )
    
    
    return(TRUE)
  }
  
  observeEvent(input$buttonImportFileUp, {
    t <- importExprsFile()
    output$expressionLoadedStatus <- renderUI( if(t) { return("... sucessfully imported file")} else {"error importing file"})
    
  })
  
  
}




