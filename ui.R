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
            div(class="test",
                id="dataLoader",
                h5("Loading the data..."),
                p("This may take a few moments"),
                img(src="ajax-loader.gif")
                ),
            div(class = "busy",  
                p("Computation in progress..."), 
                img(src="ajax-loader.gif")
                ),
            tableOutput("coordinfo"),
            tabsetPanel(type = "tabs",
                        tabPanel("Plot",
                                 textOutput('text2'),
                                 textOutput('text'),
                                 plotOutput('plot',clickId="coords",height='800px')),
                        tabPanel("Gene Table",
                                 uiOutput('selectedGenes')
                                 )
                        )
            )
        )
    ))


