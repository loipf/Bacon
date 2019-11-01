library(shiny)
library(DT)
library(plotly)

ui <- fluidPage(
  
  #fluidRow = grid columns add up to 12
  ###### load file ######
  wellPanel(
  fluidRow(
    column(4,
           fileInput("fileLoadPred", "choose tab-sep prediction file",
                     multiple = FALSE,
                     accept = c(".csv",".txt",".tsv")),
           uiOutput("textLoadSuccessful")
    ),
    column(3,
           uiOutput("buttonIfLabel"),
           uiOutput("buttonMakePred")
    ),
    column(2,
           br(),
           br()
    ),
    column(3,
      actionButton("buttonShowExampleInput", "show example input"),
      br(),
      br(),
      actionButton("buttonShowGeneLists", "show necessary gene lists",
                   onclick ="/home/proj/biocluster/praktikum/neap_ss18/neapss18_geneexp/finishedScripts/predictionTool/geneLists/', '_blank')")
    )
    
    )
  ),
  
  #### show example input ####
  conditionalPanel(condition="input.buttonShowExampleInput%2 == 1",
       wellPanel(
         fluidRow(
           column(8, h5("example input format (tab-separated, samples in rows, gene symbols in columns): ")),
           column(2, actionButton("buttonLoadExampleInput", "load (extended) example input"))
         ),
         DT::dataTableOutput("tableExampleInput")
     )
  ),
  
  
  #### show prediction output ####
  conditionalPanel(condition="output.predPossible == 'y'", 
       
    # CADi PLSR        
   wellPanel(
     h4(strong("[blood cells] partial least squares regression model of CADi with dataset GSE12288")),  #multivariate
     p("microarray data: GPL96:	[HG-U133A] Affymetrix Human Genome U133A Array"),
     p("overall moderate prediction model, needs expressionrate ~6, not adaptable to other mircoarrays ", a("[paper]", href="http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0007037",target="_blank")),

     hr(),
     verbatimTextOutput("text_cadi_plsr_overview"),
     hr(),
     
     tabsetPanel(
       tabPanel(h5(strong("own analysis results")),
                fluidRow(
                column(4, DT::dataTableOutput("table_cadi_pred_own")    ),
                column(1, br()),
                column(7, plotlyOutput("plot_cadi_perf_own"))
                       )
                
       ),
       tabPanel(h5(strong("paper analysis results")),
                fluidRow(
                  column(5, DT::dataTableOutput("table_cadi_pred_paper")    ),
                  column(7, plotlyOutput("plot_cadi_perf_paper"))
                )
                
       )
       
     )
     
  
     
     # fluidRow(
     #   column(6,align="center",
     #      h5(strong("own analysis results")),
     #      tabsetPanel(
     #        tabPanel("download data", 
     #                 DT::dataTableOutput("table_cadi_pred_own")),
     #        tabPanel("performance plot",
     #          plotOutput("plot_cadi_perf_own")),  # evtl plotly
     #        tabPanel("RMSEP plot",
     #          plotOutput("plot_cadi_rmsep_own"))  # evtl plotly
     #      )
     #   ),
     #   
     #   column(6, align="center",
     #      h5(strong("paper analysis results")),
     #      tabsetPanel(
     #        tabPanel("download data", 
     #                 DT::dataTableOutput("table_cadi_pred_paper")),
     #        tabPanel("performance plot",
     #                 plotOutput("plot_cadi_perf_paper")),  # evtl plotly
     #        tabPanel("RMSEP plot",
     #                 plotOutput("plot_cadi_rmsep_paper"))  # evtl plotly
     #      )
     #   )
     # 
     # )
    
      ), # well panel end
   
   
   # other prediction method
   wellPanel(
     h4("[tissue] prediction model 2")
   )
   
   
   
   
  ) # cond end
  


  

  
  
)


