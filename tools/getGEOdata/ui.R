library(shiny)
library(DT)

ui <- fluidPage(
  
  
  tabsetPanel(id="importMethods",
  
  tabPanel("load from GEO",
              
  #fluidRow = grid columns add up to 12
  ###### load file ######
  wellPanel(
  fluidRow(
    column(3,
      textInput("GEOacc", "GEO accession:", "GSE18275"), #GSE41571
      actionButton("buttonGetGEO", "load"),
      verbatimTextOutput("textGEOacc") 
    ),
    column(3,
       uiOutput("GEOplatform"),
       uiOutput("buttonGetPlatform")
    ),
    column(6,
      uiOutput("GEOtitle"),
      verbatimTextOutput(""),
      uiOutput("GEOfileLoaded")
    )
    )
  ),
  
  ###### choose appropriate columns ######
  conditionalPanel(condition="output.platformLoaded == 'y'",
                   
                   wellPanel(
                     fluidRow(
                       column(12,
                              verbatimTextOutput("GEOdataInfo")
                       )
                     ),
                     fluidRow(
                       column(2,
                              uiOutput("chooseSampleId")
                       ),
                       column(2,
                              uiOutput("chooseGroup")
                       ),
                       column(2,
                              uiOutput("chooseOtherColumns")
                       ),
                       column(2,
                              uiOutput("chooseGeneSymbol")
                       ),
                       column(4,
                              actionButton("buttonPlatformFinished", "import dataset"),
                              uiOutput("platformFinsishedStatus")
                       )
                     )
                   ),
                   
                   wellPanel(
                     tags$p(tags$b("shortend data preview")),
                     tabsetPanel(
                       id = 'showTables',
                       tabPanel("phenotype table", DT::dataTableOutput('showPhenoTable')),
                       tabPanel("gene feature table", DT::dataTableOutput('showFeatureTable') )
                     )
                   )
                   
                   
  )
  
              ),  #tab end
  
  
  ####### to difficult and takes to long for each format to work ..
  #
  # tabPanel("load whole dataset - impossible !!", # str+shift+c
  #          wellPanel(
  #            fluidRow(
  #              column(3,
  #                     fileInput("datasetFileUp", "choose data in .soft(.gz) format",
  #                               multiple = FALSE,
  #                               accept = c(".soft.gz", ".soft"))
  #              ),
  #              column(3,
  #                     uiOutput("datasetFileLoadedStatus")
  #              )
  #            )
  #          )
  # ),
  
  
  tabPanel("load expression matrix",
           wellPanel(
             fluidRow(
               column(4,
                      fileInput("matrixFileUp", "choose expression matrix file",
                     multiple = FALSE,
                     accept = c(".txt",".csv")),
                     tags$small("[ max 1GB ]")
               ),
               column(4,
                      fileInput("targetsFileUp", "choose targets file",
                                multiple = FALSE,
                                accept = c(".csv",".txt"))
               ),
               column(2,  
                      radioButtons("sepFileUp", "separator",
                                           choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                                           selected = ",")
               ),
               column(2,
                      actionButton("buttonImportFileUp", "import files"),
                      uiOutput("expressionLoadedStatus")
                      )

        )),
        wellPanel(
          h4("example data format"),
          tabsetPanel(
            tabPanel("expression matrix", DT::dataTableOutput("exampleTable_exprs")),
            tabPanel("target file", DT::dataTableOutput("exampleTable_targets") )
          )
          
          
        )
        
        
        ) # tab end
  
  
  
  
  ) # tab end full
  


  

  
  
)


