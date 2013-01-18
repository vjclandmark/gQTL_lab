library(shiny)

source("reactomeList.R")
pclasses = reactomeList()
nl = names(pclasses)
pclasses = lapply(pclasses, function(x) c(" ", x))
names(pclasses) = nl

shinyUI( pageWithSidebar(
  headerPanel("gQTL laboratory"),
  sidebarPanel(
    submitButton(),
    radioButtons("ready", "Configuration complete, start analysis..", c("yes", "no"),
       "no"),
    numericInput("plotSelector", "Index for plotting", value=1, min=1, max=1000),
    selectInput("population", "Choose a population:",
     choices = c("CEU", "YRI", "CEU+YRI")),
    numericInput("numPC", "Expression PCs to remove:", 0,
       min=0, max=35),
    numericInput("lbMAF", "lower bound on MAF:", 0.00,
       min=0.00, max=.49, step=.005),
    selectInput("featType", "Feature to be modeled:",
     choices = c("gene", "region", "pathway", "gene set")),
    textInput("geneName", "gene symbol", "CPNE1"),
    conditionalPanel(
       condition = "input.featType == 'pathway'", 
    h4("Note that you must manually clear unused list(s) when re-executing"),
    selectInput("pway2useA", "select the pathway (0..A-G):", 
            choices=pclasses[[1]], "none"),
    selectInput("pway2useH", "select the pathway (H-Q):", 
            choices=pclasses[[2]], "none"),
    selectInput("pway2useR", "select the pathway (R-Z):", 
            choices=pclasses[[3]], "none")
    ),
    sliderInput("radius", "cis radius in bases", value=1000,
      min=1, max=1e6)
    ),
  mainPanel( 
    tableOutput("parms"), br(),
    tableOutput("egids"),
#    uiOutput("activeControls"),
#    textOutput("pickedGene1"),
    plotOutput("pickedGene2")
  )
))

