library(shiny)
library(DT)
library(plotly)

ui <- fluidPage(
    h4("Upload Expression Data (human or mouse)"),
    wellPanel(
      fluidRow(
        column(6,
          fileInput("matrixFile", HTML("Expression Matrix [tab-separeted]"),
            multiple = FALSE,
            accept = c("text/tab-separated-values,text/plain",".tsv",".txt")),
          fileInput("targetsFile", HTML("Targets File [tab-separeted]"),
            multiple = FALSE,
            accept = c("text/tab-separated-values,text/plain",".tsv",".txt"))
        ),
        column(3,
          br()
        ),
        column(3,
          #actionButton("loadEx", "load example files"),
          br(),br(),
          actionButton("upload", "upload files"),
          br(),br(),
          uiOutput("uploadState")
        )
      )
    ),
    
    
  # Main panel for displaying outputs
  conditionalPanel(condition="output.filesLoaded == 'y'",
    wellPanel(
      fluidRow(
        column(9, 
          DT::dataTableOutput("targetsTable")
        ),
        column(3,
          selectInput("design", "Select the column containing the group information:",""),
          br(),br(),br(),
          radioButtons("inputType", "Input data type:",
            c("Microarray" = "micrArr","RNAseq" = "rna")),
          br(),
          radioButtons("inputOrganism", "Organism:",
            c("human","mouse")),
          br(),br(),br(),
          actionButton("qcPlots","QC"),br(),
          actionButton("diffExp","Predict differential expression")
        )
      )
    ),
    tags$hr()
  ),
  conditionalPanel(condition="output.qcWindow == 'y'",
    # plot panels
    tabsetPanel(id="outputPlots",
      tabPanel("Data Dispersion",
        tags$hr(),
        fluidRow(
          column(9,plotlyOutput("dispersionPlot")),
          column(3,
            radioButtons("dispersionTransform", "Transformation:",
              c("rawcounts" = "Matrix","log-transform" = "log","rlog" = "rlog","variant-stabilizing transformation" = "vst","counts per million" = "cpm"))
      ))),
      tabPanel("MA Plot",
        tags$hr(),
        plotlyOutput("MAPlot")),
      tabPanel("QQ Plot", 
        tags$hr(),
        plotlyOutput("QQPlot")),
      tabPanel("Boxplot", 
        tags$hr(),
        fluidRow(
          column(9,plotlyOutput("boxPlot")),
          column(3,
            selectInput("boxplotGroup", "Group:",
              "")
      ))),
      tabPanel("Count Density", 
        tags$hr(),
        fluidRow(
          column(9,plotlyOutput("densityPlot")),
          column(3,
            checkboxGroupInput("densityGroups", "Groups to show:",""),
            br(),
            checkboxInput("densityAverage","group-wise averages")
      ))),
      tabPanel("Cumulative Count Distribution", 
        tags$hr(),
        fluidRow(
          column(9,plotlyOutput("cumsumPlot")),
          column(3,
            checkboxGroupInput("cumsumGroups", "Groups to show:",""),
            br(),
            checkboxInput("cumsumAverage","group-wise averages")
      ))),
      tabPanel("Rankscatter Plot", 
        tags$hr(),
        plotlyOutput("rankPlot")),
      tabPanel("Sample Correlation", 
        tags$hr(),
        fluidRow(
          column(9,plotOutput("corrPlot")),
          column(3,
            radioButtons("corrMethod", "Correlation Method:",
              c("Pearson" = "pearson","Spearman" = "spearman")),
            br(),
            selectInput("corrGroup","Group by:","")
      ))),
      tabPanel("PCA", 
        tags$hr(),
        fluidRow(
          column(9,plotlyOutput("pcaPlot")),
          column(3,
            radioButtons("pcaMode", "PCA Mode:",
              c("2D","3D"))
      )))
    )
  )
)
