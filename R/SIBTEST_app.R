#' SIBTEST_app
#' 
#' This function provides a userinterface fo using SIBTEST, Crossing-SIBTEST, and POLYSIBTEST.
#' @aliases SIBTEST_app
#' @export SIBTEST_app
#' @import shiny
#' @references Chalmers, R. P. (2018). Improving the Crossing-SIBTEST statistic for detecting non-uniform DIF. Psychometrika, 83, 2, 376-386.
#' @references Chang, H. H., Mazzeo, J. & Roussos, L. (1996). DIF for Polytomously Scored Items: An Adaptation of the SIBTEST Procedure. Journal of Educational Measurement, 33, 333-353.
#' @references DIF-Pack (2021) Measured Progress. https://psychometrics.onlinehelp.measuredprogress.org/tools/dif/
#' @references Jiang, H., & Stout, W. (1998). Improved Type I Error Control and Reduced Estimation Bias for DIF Detection Using SIBTEST. Journal of Educational and Behavioral Statistics, 23(4), 291–322. https://doi.org/10.3102/10769986023004291
#' @references Li, H.-H. & Stout, W. (1996). A new procedure for detection of crossing DIF. Psychometrika, 61, 647-677.
#' @references Shealy, R. & Stout, W. (1993). A model-based standardization approach that separates true bias/DIF from group ability differences and detect test bias/DTF as well as item bias/DIF. Psychometrika, 58, 159-194.
#' @examples
#' 
#' \dontrun{
#'
#' # Input files should be csv files
#' # When selecting the suspect items a single item can be typed as 21 or 23, but a bundle of items should be typed and separated by commas starting
#' # starting with the lowest item number to highest item number.
#' # Similarly the matching set of items should be typed and separated by commas starting with the lowest item number to highest item number.
#' 
#' SIBTEST_app()
#'
#' }
SIBTEST_app <- function() {
  shinyApp(
    ui = fluidPage(
      titlePanel("SIBTEST"),
      sidebarLayout(
        sidebarPanel(
          fileInput(inputId = "file1",
                    label = "Choose Reference Group File. FILE NEEDS TO Comma Delimited with No headers (CSV FILE)",
                    accept = c(".csv")),
          tags$hr(),
          fileInput(inputId = "file2",
                    label ="Choose Focal Group File. FILE NEEDS TO Comma Delimited with No headers (CSV FILE)",
                    accept = c(".csv")),
          tags$hr(),
          radioButtons("listwise", label = "Delete Missing values (Default is False, this assumes you have complete data)",
                       choices = list("FALSE" = 1,
                                      "TRUE" = 2),selected = 1),
          tags$hr(),
          radioButtons("version", label = "What version do you want to run SIBTEST, Crossing SIBTEST or POLYSIBTEST (SIBTEST is defaulted)",
                       choices = list("SIBTEST" = 1,
                                      "Crossing SIBTEST" = 2,
                                      "POLYSIBTEST" = 3),selected = 1),
          tags$hr(),
          numericInput(inputId = "minc",
                       label = "What is the minimum total test score cell count (2 is smallest recommended by Shealy and Stout (1993)):",
                       value = 2, min = 2, max = 1000, step = 1),
          tags$hr(),
          numericInput(inputId = "cusr",
                       label = "What is average guessing parameter of the items on the test:",
                       value = 0.2, min = 0, max = .5, step = 0.01),
          tags$hr(),
          textInput(inputId = "suspect",
                    label = "What is the suspect item: (e.g. type in 21 for item 21 or type in 21,22,23 for a bundle of items)",
          ),
          tags$hr(),
          textInput(inputId = "matching",
                    label = "What are the matching subtest items: (e.g. type in 1,2,3,4,5,...,20 for items 1 through 20",
          ) ,
          tags$hr(),
          textInput(inputId = "ncat",
                    label = "What are number of categoreis fore each item: (e.g. type in 5,5,5,5,5,...,5 for items 1 through 20",
          ) ,
          tags$head(tags$script(src = "message-handler.js")),
          actionButton("do", "Run")
        ),
        mainPanel(
          tableOutput("contents"),
          tags$hr(),
          textOutput("text"),
          textOutput("text01"),
          textOutput("text02"),
          textOutput("text03"),
          textOutput("text04"),
          textOutput("text05"),
          textOutput("text06"),

          tags$hr(),
          textOutput("text1"),
          textOutput("ref_n"),
          tableOutput("reference"),

          tags$hr(),
          textOutput("text2"),
          textOutput("foc_n"),
          tableOutput("focal")
        )
      )
    ),
    server = function(input,output,session) {
      data_ref <- reactive({
        inFile <- input$file1
        if(is.null(inFile))
          return(NULL)
        df <- read.csv(inFile$datapath, header =FALSE, sep = ",")
        return(df)
      })
      data_foc <- reactive({
        inFile <- input$file2
        if(is.null(inFile))
          return(NULL)
        df <- read.csv(inFile$datapath, header =FALSE, sep = ",")
        return(df)
      })
      listwise <- reactive({as.numeric(input$listwise)})
      matching <- reactive({as.numeric(unlist(strsplit(input$matching,",")))})
      ncat <- reactive({as.numeric(unlist(strsplit(input$ncat,",")))})
      suspect <- reactive({as.numeric(unlist(strsplit(input$suspect,",")))})
      cusr <- reactive({as.numeric(input$cusr)})
      minc <- reactive({as.numeric(input$minc)})
      version <- reactive({as.numeric(input$version)})
      observeEvent(input$do, { if(version() == 1){
        output$contents <- renderTable({DIFSIB::sib(data_ref(),data_foc(),suspect_items = suspect(),matching_items = matching(), cusr = cusr(), minc = minc(),listwise = listwise())
        },digits = 3)}
        if(version() ==2){
          output$contents <- renderTable({DIFSIB::csib(data_ref(),data_foc(),suspect_items = suspect(),matching_items = matching(), cusr = cusr(), minc = minc(),listwise = listwise())
          },digits = 3)}
        if(version() == 3){
          output$contents <- renderTable({DIFSIB::psib(data_ref = data_ref(),data_foc = data_foc(),suspect_items = suspect(),matching_items = matching(),nch=ncat)
          },digits = 3)
        }
      })
      output$text1 <- renderText({
        "Reference Group Data Snapshot if data did not load correct format please check file format need comma delimited csv"})
      output$reference <- renderTable({head(data_ref())  })
      output$ref_n <- renderText({paste0("Number of Participants in Reference Group: ", as.character(nrow(data_ref())))})
      output$text2 <- renderText({
        "Focal Group Data Snapshot if data did not load correct format please check file format need comma delimited csv"})
      output$focal <- renderTable({head(data_foc())  })
      output$foc_n <- renderText({paste0("Number of Participants in Focal Group: ", as.character(nrow(data_foc())))})
      session$onSessionEnded(function() {
        stopApp()
      })
    }
  )
}
