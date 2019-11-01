library(shiny)
library(DT)
library(plotly)

ui <- fluidPage(

    h4("load gene list (human or mouse)"),
    wellPanel(
      fluidRow(
        column(4,
          textAreaInput("geneText", label = "paste gene list [comma-separated]", value = "", resize="vertical", height='90px')
        ),
        column(3,
          fileInput("geneFile", HTML("file with gene list<br/>[comma-separeted]"),
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv"))
        ),
        column(3,
               textInput('strainInput', HTML('specify mouse strain<br/>[MGI, mutation or genotype]')),
               checkboxInput("exactMouseName", "exact name matching", FALSE)

          ),
        column(1,
          actionButton("buttonGetTables", "find genes"),
          uiOutput("textFindGenes")
        )
      )
      
    ),
    

    
    # Main panel for displaying outputs
    conditionalPanel(condition="output.geneListLoaded == 'y'",
      wellPanel(
        fluidRow(
          column(5, verbatimTextOutput("textOutputStatistic")),
          column(2, br()),
          column(2, actionButton("buttonShowPlots", "show summary plots")
          ) ),
        tags$hr(),
        
        # plot display
        conditionalPanel(condition="input.buttonShowPlots%2 == 1",
             h5("summary of search:"),
             fluidRow(
               column(6, plotlyOutput("barPlot")),
               column(6, plotlyOutput("pieChart"))
             ),
             tags$hr()
                  
        ),
        
        # table panels
        tabsetPanel(
          id="outputTables",
          tabPanel("all genes", 
                   tags$hr(),
                   DT::dataTableOutput("tableAll")),
          tabPanel("only ath specifc genes [own]",
                   tags$hr(),
                   DT::dataTableOutput("tableOwnAth")),
          tabPanel("only ath specifc genes [MalaCard]", 
                   tags$hr(),
                   DT::dataTableOutput("tableMalaCardAth"))
        )
      )
    )

  
  
  
)


