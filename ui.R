library(shinyIncubator)
library(shiny)

## ui.R
shinyUI(fluidPage(
    progressInit(),
    titlePanel("TSS Visualizer"),
    
    sidebarLayout(
        sidebarPanel(
            helpText(h5("Display coverages around the TSS")),
            uiOutput("typeSelector"),
            uiOutput("markerSelector"),
            uiOutput("rangeSelector"),
            hr(),
            uiOutput("sampleSelector"),
            uiOutput("plotButton")
            ),
        mainPanel(
            tableOutput("coordinfo"),
            tabsetPanel(type = "tabs",
                        tabPanel("Plot",
                                 textOutput('text'),
                                 plotOutput('plot',clickId="coords",height='800px')),
                        tabPanel("Gene Table",
                                 plotOutput('metaplot'),
                                 dataTableOutput('table')
                                 )
                        )
            )
        )
    ))

