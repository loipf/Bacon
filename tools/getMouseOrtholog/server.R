library(shiny)
library(DT)
library(plotly)



### TO DO ###
# implement own differential expressed genes
# specify Rscript path down there
# ath and mala table not automatically reseted

# execute before running app
# setwd('/home/stefan/bioinformatics/neap/getMouseOrtholog/')
# shiny::runApp('/home/stefan/bioinformatics/neap/getMouseOrtholog', port = 7817) # dont know


path = "/home/stefan/bioinformatics/neap/getMouseOrtholog/"
#path = "/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/getMouseOrtholog/"
path = "/Users/stefanloipfinger/Documents/GitHub/Bacon/tools/getMouseOrtholog/"



# global variables
fileLoaded <- FALSE
geneList <- ""
mouseStrainList <- 'all'


server <- function(input, output) {
  #options(shiny.maxRequestSize=1000*1024^2)  # upload up to 1GB of data


  
  # find genes
  loadTables <- eventReactive(input$buttonGetTables, {
    #reset everything
    fileLoaded <<- FALSE
    geneList <<- ''
    mouseStrainList <<- 'all'
    
    # check if GEOacc is already in database !!! 
    tryCatch( {
      withProgress(message = "search for orthologs", value = 0, {
      incProgress(1/10) # loading bar
        
        #delete old ones
        if (file.exists(paste0(path,"fullAthGenes.csv"))) file.remove(paste0(path,"fullAthGenes.csv"))
        if (file.exists(paste0(path,"ownAthGenes.csv"))) file.remove(paste0(path,"ownAthGenes.csv"))
        if (file.exists(paste0(path,"malaAthGenes.csv"))) file.remove(paste0(path,"malaAthGenes.csv"))
  
      geneList <<- input$geneText
      incProgress(5/10) # loading bar
      if (!is.null(input$geneFile)) {
          geneListFile <- readChar(input$geneFile$datapath, file.info(input$geneFile$datapath)$size)
          if(geneList!="") {  geneList <<- paste(c(geneList,geneListFile), collapse=",") }
          else { geneList <<- geneListFile}
      }
      fileLoaded <<- TRUE
      geneList <<- gsub("\\s+", "", geneList, perl=TRUE)
      incProgress(10/10) # loading bar
      })
      if(geneList=="") {
        fileLoaded <<- FALSE
        return(FALSE)
      } else {
        return(TRUE)
      }
      
    },
      error = function(e) {
        fileLoaded <<- FALSE 
      return(FALSE)}
    )
  })
  
  output$textFindGenes <- renderUI({
    if (loadTables()) {
      return("... sucessfully loaded")
    } 
    else {
      return("error loading list/file")
    }
  })
  


  # if file loaded show tables
  output$geneListLoaded <- eventReactive(input$buttonGetTables, {
    if(fileLoaded) {
      
      # run search and filter
      if(input$strainInput!="") {
        mouseStrainList <<- input$strainInput
      }
      
      system(paste0("Rscript ",path,"getAthGenes.R"," ",geneList," ", mouseStrainList," ", input$exactMouseName))
      #commandArgs <- function(...) geneList, mouseStrainList, input$exactMouseName
      #source('file_to_source.R')
      
      
      allTable <- read.table(paste0(path,"fullAthGenes.csv"), sep = ",", header=TRUE)
      #output$downloadAll <- renderUI({ downloadButton("downloadAll","download file") })
      output$tableAll <- DT::renderDataTable(
        DT::datatable(allTable,extensions='Buttons', filter = 'bottom',options = list(searching = TRUE,pageLength = 10, dom = "Blfrtip", buttons = list("copy", list( extend = "collection" , buttons = c("csv", "excel", "pdf"), text = "Download" ) )
        ) ))
      #output$downloadAll <- downloadHandler(filename = "fullAthGenes.tsv", content = function(file) { write.table(allTable, file, sep="\t", row.names=FALSE, quote=FALSE) } )
      
      ownTable <- NULL
      ownGenesFound <- 0
      if (file.exists(paste0(path,"ownAthGenes.csv"))) {
        #output$downloadOwn <- renderUI({ downloadButton("downloadOwn","download file") })
        ownTable <- read.table(paste0(path,"ownAthGenes.csv"), sep = ",", header=TRUE)
        output$tableOwnAth <- DT::renderDataTable(
          DT::datatable(ownTable, extensions='Buttons',filter = 'bottom', options = list(searching = TRUE,pageLength = 10, dom = "Blfrtip", buttons = list("copy", list( extend = "collection" , buttons = c("csv", "excel", "pdf"), text = "Download" ) )
          ) ))
        #output$downloadOwn <- downloadHandler(filename = "ownAthGenes.tsv", content = function(file) { write.table(ownTable, file, sep="\t", row.names=FALSE, quote=FALSE) } )
        ownGenesFound <- nrow(ownTable)
      }
      
      malaTable <- NULL
      malacardGenesFound <- 0
      if (file.exists(paste0(path,"malaAthGenes.csv"))) {
        #output$downloadMala <- renderUI({ downloadButton("downloadMala","download file") })
        malaTable <- read.table(paste0(path,"malaAthGenes.csv"), sep = ",", header=TRUE)
        output$tableMalaCardAth <- DT::renderDataTable(
          DT::datatable(malaTable, extensions='Buttons', filter = 'bottom', options = list(searching = TRUE,pageLength = 10, dom = "Blfrtip", buttons = list("copy", list( extend = "collection" , buttons = c("csv", "excel", "pdf"), text = "Download" ) )
        ) ))
        #output$downloadMala <- downloadHandler(filename = "malaAthGenes.tsv", content = function(file) { write.table(malaTable, file, sep="\t", row.names=FALSE, quote=FALSE) } )
        malacardGenesFound <- nrow(malaTable)
      }
      
      output$textOutputStatistic <- renderText({ 
        paste0("searched genes: ", nrow(allTable),"\northologs found: ", max(sum(allTable$geneID_mouse!=""),0,na.rm = TRUE),
               "\natherosclerosis connection own: ",ownGenesFound,"\natherosclerosis connection MalaCard: ", malacardGenesFound) })
      
      
      createPlots(allTable,ownTable,malaTable)
      
      return("y") # to show panel
      
      
      }
  })
  
  
  createPlots <- function(tAll, tOwn,tMala) {
    #bar plot
    tAll[] <- lapply(tAll, as.character)
    tAll[is.na(tAll)] <- ""
    
    tAllpercNotFound <- as.numeric(max(sum(tAll$mouseHumanConnection==""),0,na.rm = TRUE)/nrow(tAll))
    tAllTextNotFound <- as.character(tAll[tAll$mouseHumanConnection=="",]$input)
    tAllTextNotFound <- paste0(length(tAllTextNotFound)," genes:",paste(tAllTextNotFound,collapse=","))
    
    # f*** the bug, hardcoded
    tAllpercNotFound <- if(length(tAllTextNotFound)==0) 1 else tAllpercNotFound
    
    tAllpercFound <- as.numeric(1-tAllpercNotFound)
    barDf <- data.frame(kindof=c("all","all"),found=c("found","not_found"), perc=c(tAllpercFound,tAllpercNotFound), txt=c("",tAllTextNotFound))

    # not implemented
    # if(!is.null(tOwn)){
    #   tOwnpercFound <- (1-max(sum(tOwn$geneID_mouse==""),0,na.rm = TRUE)/nrow(tOwn))
    #   tOwnpercNotFound <- (max(sum(tOwn$geneID_mouse==""),0,na.rm = TRUE)/nrow(tOwn))
    #   tOwnTextNotFound <- as.character(tOwn[tOwn$geneID_mouse=="",]$geneID_human)
    #   tOwnTextNotFound <- paste0(length(tOwnTextNotFound)," genes:",paste(tOwnTextNotFound,collapse=","))
    #   barDfown <- data.frame(kindof=c("own","own"), found=c("found","not_found"),perc=c(tOwnpercFound,tOwnpercNotFound),txt=c("",tOwnTextNotFound))
    #   barDf <- rbind(barDf,barDfown)
    # }
    # 
    # if(!is.null(tMala)){
    #   tMalapercFound <- (1-max(sum(tMala$geneID_mouse==""),0,na.rm = TRUE)/nrow(tMala))
    #   tMalapercNotFound <- (max(sum(tMala$geneID_mouse==""),0,na.rm = TRUE)/nrow(tMala))
    #   tMalaTextNotFound <- as.character(tMala[tMala$geneID_mouse=="",]$geneID_human)
    #   tMalaTextNotFound <- paste0(length(tMalaTextNotFound)," genes:",paste(tMalaTextNotFound,collapse=","))
    #   barDfMala <- data.frame(kindof=c("MalaCard","MalaCard"), found=c("found","not_found"),perc=c(tMalapercFound,tMalapercNotFound),txt=c("",tMalaTextNotFound))
    #   barDf <- rbind(barDf,barDfMala)
    # }
    
  output$barPlot <- renderPlotly({
    plot_ly(barDf, x = ~kindof, y = ~perc*100, type = 'bar', color=~found,colors=c("green","grey"), text=~txt,hoverinfo = 'text') %>% 
      layout(title='orthologs found', yaxis = list(title = 'amount of genes [%]'),barmode = "stack", xaxis=list(title="")) 
  } )
  
  
  # pie chart plot, get all intersections and plot
  pieDf <- data.frame(kindof=character(),number=double(),txt=character())
  genesAll <- as.character(tAll[tAll$geneID_mouse!="",]$geneID_mouse)
  interAll <- genesAll
  
  if(!is.null(tOwn)){
    genesOwn <- as.character(tOwn[tOwn$geneID_mouse!="",]$geneID_mouse)
    interOwn <- intersect(genesAll,genesOwn)
    interAll <- interAll[! interAll %in% genesOwn]
    txtInterOwn <- paste0(length(interOwn)," genes:",paste(interOwn,collapse=","))
    b <- data.frame(kindof="own", number= length(interOwn), txt=txtInterOwn)
    pieDf <- rbind(pieDf,b)
  }
  if(!is.null(tMala)){
    genesMala <- as.character(tMala[tMala$geneID_mouse!="",]$geneID_mouse)
    interMala <- intersect(genesAll,genesMala)
    interAll <- interAll[! interAll %in% genesMala]
    txtInterMala <- paste0(length(interMala)," genes:",paste(interMala,collapse=","))
    b <- data.frame(kindof="MalaCard", number= length(interMala), txt=txtInterMala)
    pieDf <- rbind(pieDf,b)
  }
  if(!is.null(tMala) & !is.null(tOwn)){
    pieDf <- pieDf[pieDf$kindof!="own"]
    pieDf <- pieDf[pieDf$kindof!="MalaCard"]
    
    interOwnMala <- intersect(genesMala,genesOwn)
    txtInterOwnMala <- paste0(length(interOwnMala)," genes:",paste(interOwnMala,collapse=","))
    b <- data.frame(kindof="own + MalaCard", number= length(interOwnMala), txt=txtInterOwnMala)
    pieDf <- rbind(pieDf,b)
    
    interOnlyOwn <- interOwn[! interOwn %in% genesMala]
    txtInterOnlyOwn <- paste0(length(interOnlyOwn)," genes:",paste(interOnlyOwn,collapse=","))
    b <- data.frame(kindof="only own", number= length(interOnlyOwn), txt=txtInterOnlyOwn)
    pieDf <- rbind(pieDf,b)
    
    interOnlyMala <- interMala[! interMala %in% genesOwn]
    txtInterOnlyMala <- paste0(length(interOnlyMala)," genes:",paste(interOnlyMala,collapse=","))
    b <- data.frame(kindof="only MalaCard", number= length(interOnlyMala), txt=txtInterOnlyMala)
    pieDf <- rbind(pieDf,b)
  }
  
  if(length(interAll)<100) {
    txtInterAll <- paste0(length(interAll)," genes:",paste(interAll,collapse=","))
  } else {
    txtInterAll <- paste0(length(interAll)," genes:",paste(interAll[1:100],collapse=",")," ...")
  }
  
  b <- data.frame(kindof="other ortholog genes", number= length(interAll), txt=txtInterAll)
  pieDf <- rbind(pieDf,b)
  pieDf$perc <- (pieDf$number/sum(pieDf$number))*100
  
  output$pieChart <- renderPlotly({
    plot_ly(pieDf, labels=~kindof,textinfo='label+percent', values=~perc, type='pie', textposition='inside', text=~txt, hoverinfo="text") %>% 
            layout(title='ortholog genes against atherosclerosic related ones',
                   xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
  } )
  
  }
  
  
  
  outputOptions(output,"geneListLoaded", suspendWhenHidden = FALSE)
  
  
}




