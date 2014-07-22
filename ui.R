#library(shinyIncubator)
library(shiny)

## ui.R
shinyUI(fluidPage(
    tagList(
        tags$head(
            tags$link(rel="stylesheet", type="text/css",href="style.css"),
            tags$script(type="text/javascript", src = "busy.js")
            )
        ),
    titlePanel("TSS Visualizer"),
    sidebarLayout(
        sidebarPanel(
            helpText(h5("Display coverages around the TSS")),
            uiOutput("selectors")
            ),
        mainPanel(
            textOutput("debug"),
            div(class="test",
                id="dataLoader",
                h5("Loading the data, this may take a while"),
                p("Be patient..."),
                img(src="ajax-loader.gif")
                ),
            div(class = "busy",  
                p("Calculation in progress.."), 
                img(src="ajax-loader.gif")
                ),
            tableOutput("coordinfo"),
            tabsetPanel(type = "tabs",
                        tabPanel("Plot",uiOutput('plotRegion')),
                        tabPanel("Gene Table",uiOutput('selectedGenes'))
                        )
            )
        )
    ))


