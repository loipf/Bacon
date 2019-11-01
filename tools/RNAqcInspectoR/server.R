library(shiny)
library(DT)
library(ggplot2)
library(plotly)
library(LSD)
library(pheatmap)
library(limma)
#library(edgeR)
library(DESeq2)
library(hexbin)

# execute before running app
# setwd('/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/RNAqcInspectoR/')
# shiny::runApp('/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/RNAqcInspectoR/', port = 7817) # dont know

path = "/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/RNAqcInspectoR/"
#path = "/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/RNAqcInspectoR/"


# global variables
fileLoaded <- F
dat <- NULL # DGE data object
method <- "rna"
organism <- "human"

server <- function(input,output,session) {
  #options(shiny.maxRequestSize=1000*1024^2)  # upload up to 1GB of data
  
  # load data of an RNAseq experiment and store it to an internal expression data object
  # 
  # dge_file: path to the DGE Matrix File ("DGE_Matrix.txt" by default)
  # targets_file: path to targets file containing information about the experiment (e.g. group descriptors; "Targets.txt" by default)
  # rowMin: minimal number of counts within one row (gene); the default value is 2
  # 
  # value: a list with three elements: "Targets" contains experimental information, rowMin the minimal number of counts within one row,
  # "Matrix" the raw count matrix
  loadData <- function(dge_file="DGE_Matrix.txt",targets_file="Targets.txt",rowMin=2){
    dge = read.csv(dge_file,header=T,sep='\t')
    # replace any '+' in the targets file by '.' ('+' in a name cannot be handled in R)
    targets = as.data.frame(apply(read.csv(targets_file,header=T,sep='\t'),2,function(x){x=gsub("\\+","\\.",x)}),stringsAsFactors=F)
    # if there are any fields with label 'Control', change these to 'control' (the default value for reference samples)
    for(i in 1:nrow(targets)){
      for(j in 1:ncol(targets)){
        if(targets[i,j] == "Control") targets[i,j] = "control"
      }
    }
    colnames(targets)[1] = "Sample"
    rownames(targets) = targets$Sample
    
    colnames(dge)[1] = "GENE"
    rownames(dge) = dge$GENE
    dge = dge[,colnames(dge)!='GENE']
    dge = dge[,targets$Sample]
    dge = dge[apply(dge,1,sum) >= rowMin,colnames(dge) %in% targets$Sample]
    dge = round(dge)
    
    return(list(Targets=targets,rowMin=rowMin,Matrix=dge))
  }
  
  
  # load DGE files
  loadFiles <- eventReactive(input$upload, {
    #reset everything
    fileLoaded <<- F
    dat <<- NULL
    method <<- input$inputType
    organism <<- input$inputOrganism
    
    # try to load files
    tryCatch({
      withProgress(message="loading files...", value = 0, {
        incProgress(1/10) # loading bar
        if(!is.null(input$matrixFile) & !is.null(input$targetsFile)){
          dat <<- loadData(input$matrixFile$datapath,input$targetsFile$datapath)
          fileLoaded <<- T
        }
        
        incProgress(5/10) # loading bar
        
        incProgress(10/10) # loading bar
      })
      
      if(is.null(dat)) {
        fileLoaded <<- F
        return(F)
      } else{
        return(T)
      }
    },
    error = function(e){
      fileLoaded <<- F 
      return(F)}
    )
  })
  
  
  # update upload state
  output$uploadState <- renderUI({
    if(loadFiles()) {
      return("... sucessfully loaded")
    } 
    else{
      return("error loading files")
    }
  })
  
  
  # if files are loaded succesfully, show targets table
  output$filesLoaded <- eventReactive(input$upload, {
    if(fileLoaded){
      output$targetsTable <- DT::renderDataTable(
        DT::datatable(dat$Targets,filter='bottom',options=list(searching=T,pageLength=10,dom="Blfrtip")
      ))
      updateSelectInput(session,"design",choices=colnames(dat$Targets))
      updateSelectInput(session,"corrGroup",choices=colnames(dat$Targets))
      
      return("y") # to show panel
    }
  })
  
  
  output$qcWindow <- eventReactive(input$qcPlots, {
    updateCheckboxGroupInput(session,"densityGroups",choices=unique(dat$Targets[[input$design]]))
    updateCheckboxGroupInput(session,"cumsumGroups",choices=unique(dat$Targets[[input$design]]))
    return("y") # to show panel
  })
  
  
  output$densityPlot <- renderPlotly({
    if(!is.null(input$densityGroups)){
      targets = dat$Targets[dat$Targets[[input$design]]%in%input$densityGroups,]
      
      if(input$densityAverage){
        mat = list()
        for(i in input$densityGroups){
          tmp = targets[targets[[input$design]]==i,]
          mat[[i]] = apply(log2(dat$Matrix[,colnames(dat$Matrix)%in%rownames(tmp)]+1),1,mean)
        }
        mat = as.data.frame(mat)
      } else{
        mat = log2(dat$Matrix[,colnames(dat$Matrix)%in%rownames(targets)]+1)
      }
      mat = stack(mat)
      
      p <- ggplot(mat,aes(x=values)) +
        stat_density(aes(group=ind,color=ind),position="identity",geom="line")
      p <- ggplotly(p)
    } else plotly_empty()
  })
  
  
  output$cumsumPlot <- renderPlotly({
    if(!is.null(input$cumsumGroups)){
      targets = dat$Targets[dat$Targets[[input$design]]%in%input$cumsumGroups,]
      p <- plot_ly()
      
      if(input$cumsumAverage){
        mat = list()
        for(i in input$cumsumGroups){
          tmp = targets[targets[[input$design]]==i,]
          mat[[i]] = apply(log2(dat$Matrix[,colnames(dat$Matrix)%in%rownames(tmp)]+1),1,mean)
        }
        mat = as.data.frame(mat)
      } else{
        for(i in rownames(targets)){
          xvalue = log2(dat$Matrix[,colnames(dat$Matrix)==i]+1)
          xvalue = xvalue[order(xvalue)]
          cumul <- cumsum(xvalue)/sum(xvalue)
          p <- p %>% add_trace(x=xvalue,y=cumul,type="scatter",mode="lines",name=i)
        }
      }
      p %>% layout(xaxis=list(title="log2(readcont)"),yaxis=list(title="cumulative density"))
    } else plotly_empty()
  })
  
  
  output$corrPlot <- renderPlot({
    if(fileLoaded){
      df=data.frame(dat$Targets[,input$corrGroup])
      rownames(df) = rownames(dat$Targets); colnames(df) = input$corrGroup
      pheatmap(cor(dat$Matrix,method=input$corrMethod),clustering_distance_rows='correlation',clustering_distance_cols='correlation',clustering_method='average',annotation_col=df,color=colorRampPalette(rev(c('red','yellow')))(100))
    }
  })
  
  
  output$pcaPlot <- renderPlotly({
    if(fileLoaded){
      mode = input$pcaMode
      col.code = dat$Targets[[input$design]]
      color = rainbow(max(as.numeric(as.factor(col.code))))
      
      pca = prcomp(dat$Matrix,center=T)
      out = data.frame(pca$rotation,col.code,txt=rownames(dat$Targets))
      percentage = signif(pca$sdev^2/sum(pca$sdev^2)*100,2)
      
      if(mode=="2D"){
        plot_ly(out,x=~PC1,y=~PC2,color=~col.code,colors=color,text=~txt,hoverinfo='text') %>%
          layout(title="PCA (2D)",xaxis=list(title=paste0("PC1 (",percentage[1]," % variance)")),yaxis=list(title=paste0("PC2 (",percentage[2]," % variance)")))
      } else{
        plot_ly(out,x=~PC1,y=~PC2,z=~PC3,color=~col.code,colors=color,text=~txt,hoverinfo='text') %>%
          layout(title="PCA (3D)",xaxis=list(title=paste0("PC1 (",percentage[1]," % variance)")),yaxis=list(title=paste0("PC2 (",percentage[2]," % variance)"),zaxis=list(title=paste0("PC3 (",percentage[3]," % variance)"))))
      }
    }
  })
  
  
  outputOptions(output,"filesLoaded",suspendWhenHidden=F)
  outputOptions(output,"qcWindow",suspendWhenHidden=F)
}