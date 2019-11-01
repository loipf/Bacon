library(shiny)
library(DT)
library(plotly)
library(crosstalk)


scriptPath <- ("/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/predictionTool/")
scriptPath <- ("/Users/stefanloipfinger/Documents/GitHub/Bacon/tools/predictionTool/")
source(paste0(scriptPath,"pred_plsrCADi/plsrCADi.R"))

### TO DO ###
# specify path to gene lists, does not work
# check if input is right for pred
# check if binary or CADi prediction
# add statistic


# execute before running app
# setwd("/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/predictionTool")
# shiny::runApp('/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/predictionTool', port = 7818) # dont know



# global variables
predFile <- NULL
predMatrix <- NULL # without classes
knownPredLabel <- NULL



server <- function(input, output) {
  options(shiny.maxRequestSize=1000*1024^2)  # upload up to 1GB of data

  ##### load all prediction methods #####
  loadAllModels <- function() {
    plsr_initialize()
  }
  loadAllModels() # load everything needed
  
  
  ##### handle pred input #####
  observeEvent(input$fileLoadPred, {
    tryCatch(
      { 
        # reset
        predFile <<- NULL
        knownPredLabel <<- NULL
        
        predFile <<- read.table(input$fileLoadPred$datapath, header = TRUE)
        output$textLoadSuccessful <- renderUI(return("... sucessfully loaded file"))
        cols = colnames(predFile)[2:length(predFile)]
        output$buttonIfLabel <- renderUI({ selectInput("buttonIfLabel","compare with known class?",choices = c("-" ="-",cols))  })
        output$buttonMakePred <- renderUI({ actionButton("buttonMakePred", "start predictions") })
      },
      error = function(e) { 
        output$textLoadSuccessful <- renderUI(return("error loading file"))
    })
  } )
  
  
  
  
  
  ##### show example table #####
  observeEvent(input$buttonShowExampleInput, {
    output$tableExampleInput <- DT::renderDataTable(DT::datatable(
      read.table(paste0(scriptPath,"show_example_pred.txt"), sep="\t", header=TRUE ),
      options = list(searching = FALSE,pageLength = 10,lengthChange = FALSE )  # hide columns: in list columnDefs = list(list(visible=FALSE, targets=2:145))
    ))
  })
  
  observeEvent(input$buttonLoadExampleInput, {
    predFile <<- read.table(paste0(scriptPath,"example_cadi_pred.txt"), sep="\t", header=TRUE )
    output$textLoadSuccessful <- renderUI(return("example loaded, please start prediction"))
    cols = colnames(predFile)[2:length(predFile)]
    output$buttonIfLabel <- renderUI({ selectInput("buttonIfLabel","compare with known class?",choices = c("-",cols))  })
    output$buttonMakePred <- renderUI({ actionButton("buttonMakePred", "start predictions") })
    
  })
  
  
  
  
  ############################
  ##### make predictions #####
  
  ### load test
   predFile <<- read.table(paste0(scriptPath,"example_cadi_pred.txt"), header = TRUE)
  cols = colnames(predFile)[2:length(predFile)]
  output$buttonIfLabel <- renderUI({ selectInput("buttonIfLabel","compare with known class?",choices = c("-",cols))  })
  output$buttonMakePred <- renderUI({ actionButton("buttonMakePred", "start predictions") })
  #
  ###
  
  # main prediction method
  output$predPossible <- eventReactive(input$buttonMakePred, {
    # handle labels if given
    predMatrix <<- predFile
    if(input$buttonIfLabel != "-") {
      knownPredLabel <<- predFile[,input$buttonIfLabel]
      predMatrix <<- predFile[!colnames(predFile) %in% input$buttonIfLabel]
    }
    predGeneList <- colnames(predMatrix)    # modify data and make first row rownames
    rownames(predMatrix) <<- predMatrix[,1]
    predMatrix <<- predMatrix[,-1]
    
    #print(predMatrix)
    # start pred

    # CADi plsr panel
    cadIntroText <- plsr_compareGeneLists(predGeneList)
    output$text_cadi_plsr_overview <- renderText({ cadIntroText })
    #output$plot_cadi_own <- renderPlotly({ plsr_plot_own() })
    createPLSRcadiPanel()
    
    
    
    
    
    
    return("y") # to silly for boolean
  })
  
  outputOptions(output,"predPossible", suspendWhenHidden = FALSE) # show panel
  
  # create prediction output PLSR and plots
    createPLSRcadiPanel <- function(){
      if (cadi_ownPossible>0) {  # own pred possible
        cadi_predTable_own <- plsr_predTable_own(cadi_currOwnModel,predMatrix)
        cadi_predTable_own <- as.data.frame(cadi_predTable_own)
          colnames(cadi_predTable_own) <- c("predicted")
        if(!is.null(knownPredLabel)) {
          cadi_predTable_own <- cbind(cadi_predTable_own,knownPredLabel)
          colnames(cadi_predTable_own) <- c("predicted","measured")
        }
        
        cadi_predTable_own <- as.data.frame(cadi_predTable_own)
        cadi_predTable_own$predicted <- round(cadi_predTable_own$predicted)
        cadi_predTable_own$sampleId <- rownames(predMatrix)
        cadi_predTable_own <- cadi_predTable_own[,c(ncol(cadi_predTable_own),1:(ncol(cadi_predTable_own)-1))]
        output$table_cadi_pred_own <- DT::renderDataTable(DT::datatable(cadi_predTable_own,rownames= FALSE, extensions='Buttons', filter = 'bottom',options = list(lengthChange = FALSE,searching = FALSE,pageLength = 10, dom = "Blfrtip", buttons = list("copy", list( extend = "collection" , buttons = c("csv", "excel", "pdf"), text = "Download" ) ))))
        
        plotDfOwn <- data.frame(sampleId=rownames(cadi_currOwnModel$model), meas=cadi_currOwnModel$model$class, pred=cadi_currOwnModel$fitted.values[,,2], group="training")
        if(!"measured" %in% colnames(cadi_predTable_own))
        { meas_cadi_predTable_own <- data.frame(sampleId=cadi_predTable_own$sampleId,meas=0, pred=cadi_predTable_own$predicted, group="predicted", row.names=rownames(cadi_predTable_own)) }
        else {
          meas_cadi_predTable_own <-  data.frame(sampleId=cadi_predTable_own$sampleId,meas=cadi_predTable_own$measured, pred=cadi_predTable_own$predicted, group="predicted")
        }
        plotDfOwn <- rbind(plotDfOwn, meas_cadi_predTable_own)
        plotDfOwn$txt <- paste0(plotDfOwn$sampleId,"<br>pred: ",round(plotDfOwn$pred),", meas: ",plotDfOwn$meas)
        output$plot_cadi_perf_own <- renderPlotly({
            selectedP <- input$table_cadi_pred_own_rows_selected
            if(!length(selectedP)) {
            p <- plotDfOwn %>% plot_ly() %>% add_trace(x=~meas, y=~pred, color = ~factor(group),colors=c("black","grey","red"), text=plotDfOwn$txt, hoverinfo='text' ) %>% 
                layout(title='prediction of CADi', yaxis = list(title = 'predicted'), xaxis=list(title="measured")) %>%
                add_trace(x=c(0,100), y=c(0,100), mode="lines",line = list(color = 'rgba(26, 26, 26, .4)'), hoverinfo="none", showlegend = FALSE)
                
            #%>% highlight("plotly_selected", color = I('blue'), selected = attrs_selected(name = 'selected'))
            } else if (length(selectedP) ){
              pp <- plotDfOwn %>% plot_ly() %>% add_trace(x=~meas, y=~pred, color = ~factor(group),colors=c("black","red"), text=plotDfOwn$txt, hoverinfo='text' ) %>% 
                layout(title='prediction of CADi', yaxis = list(title = 'predicted'), xaxis=list(title="measured"))
              pp<- add_trace(pp, data = meas_cadi_predTable_own[selectedP, , drop = F], x=~meas, y=~pred,marker = list(size = 15,color = 'rgba(255, 255, 0, .5)',line = list(color = 'rgba(255,0,0, 0.8)', width = 2)),colors = I('blue'), name = 'selected') %>%
                add_trace(x=c(0,100), y=c(0,100), mode="lines",line = list(color = 'rgba(26, 26, 26, .4)'), hoverinfo="none", showlegend = FALSE)
          }
        })
      }
      
      if (cadi_paperPossible>0) {  # paper pred possible
        cadi_predTable_paper <- plsr_predTable_paper(cadi_currPaperModel,predMatrix)
        cadi_predTable_paper <- as.data.frame(cadi_predTable_paper)
        colnames(cadi_predTable_paper) <- c("predicted")
        if(!is.null(knownPredLabel)) {
          cadi_predTable_paper <- cbind(cadi_predTable_paper,knownPredLabel)
          colnames(cadi_predTable_paper) <- c("predicted","measured")
        }
        
        cadi_predTable_paper <- as.data.frame(cadi_predTable_paper)
        cadi_predTable_paper$predicted <- round(cadi_predTable_paper$predicted)
        cadi_predTable_paper$sampleId <- rownames(predMatrix)
        cadi_predTable_paper <- cadi_predTable_paper[,c(ncol(cadi_predTable_paper),1:(ncol(cadi_predTable_paper)-1))]
        output$table_cadi_pred_paper <- DT::renderDataTable(DT::datatable(cadi_predTable_paper,rownames= FALSE, extensions='Buttons', filter = 'bottom',options = list(lengthChange = FALSE,searching = FALSE,pageLength = 10, dom = "Blfrtip", buttons = list("copy", list( extend = "collection" , buttons = c("csv", "excel", "pdf"), text = "Download" ) ))))
        
        plotDfPaper <- data.frame(sampleId=rownames(cadi_currPaperModel$model), meas=cadi_currPaperModel$model$class, pred=cadi_currPaperModel$fitted.values[,,2], group="training")
        if(!"measured" %in% colnames(cadi_predTable_paper))
        { meas_cadi_predTable_paper <- data.frame(sampleId=cadi_predTable_paper$sampleId,meas=0, pred=cadi_predTable_paper$predicted, group="predicted", row.names=rownames(cadi_predTable_paper)) }
        else {
          meas_cadi_predTable_paper <-  data.frame(sampleId=cadi_predTable_paper$sampleId,meas=cadi_predTable_paper$measured, pred=cadi_predTable_paper$predicted, group="predicted")
        }
        plotDfPaper <- rbind(plotDfPaper, meas_cadi_predTable_paper)
        plotDfPaper$txt <- paste0(plotDfPaper$sampleId,"<br>pred: ",round(plotDfPaper$pred),", meas: ",plotDfPaper$meas)
        output$plot_cadi_perf_paper <- renderPlotly({
            selectedP <- input$table_cadi_pred_paper_rows_selected
            if(!length(selectedP)) {
              p <- plotDfPaper %>% plot_ly() %>% add_trace(x=~meas, y=~pred, color = ~factor(group),colors=c("black","grey","red"), text=plotDfPaper$txt, hoverinfo='text' ) %>% 
                layout(title='prediction of CADi', yaxis = list(title = 'predicted'), xaxis=list(title="measured")) %>%
                add_trace(x=c(0,100), y=c(0,100), mode="lines",line = list(color = 'rgba(26, 26, 26, .4)'), hoverinfo="none", showlegend = FALSE)
              
              #%>% highlight("plotly_selected", color = I('blue'), selected = attrs_selected(name = 'selected'))
            } else if (length(selectedP) ){
              pp <- plotDfPaper %>% plot_ly() %>% add_trace(x=~meas, y=~pred, color = ~factor(group),colors=c("black","red"), text=plotDfPaper$txt, hoverinfo='text' ) %>% 
                layout(title='prediction of CADi', yaxis = list(title = 'predicted'), xaxis=list(title="measured"))
              pp<- add_trace(pp, data = meas_cadi_predTable_paper[selectedP, , drop = F], x=~meas, y=~pred,marker = list(size = 15,color = 'rgba(255, 255, 0, .5)',line = list(color = 'rgba(255,0,0, 0.8)', width = 2)),colors = I('blue'), name = 'selected') %>%
                add_trace(x=c(0,100), y=c(0,100), mode="lines",line = list(color = 'rgba(26, 26, 26, .4)'), hoverinfo="none", showlegend = FALSE)
            }
          })
      }
  } # createPLSRcadiPanel end

}


