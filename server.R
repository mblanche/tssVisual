library(shinyIncubator)

loading.time <- system.time({
    source("helpers.R")
    source("dataLoader.R")
})

shinyServer(function(input, output, session) {

    ## Render a widget for selecting the sample to display
    observe ({
        if (exists('Relative')){
            output$selectors <- renderUI({
                list(
                    hr(),
                    ## Render the abolute vs relative radio
                    radioButtons('analysisType','Select Coverage Type',c('Relative','Absolute')),
                    ## Render a button to add/remove TSS marker
                    checkboxInput('addCenterMarker','Display Center Marker',TRUE),
                    hr(),
                    ## Add an uiOutput do dynamicaly display ROI files if more than one
                    uiOutput('ROIselector'),
                    ## Render the sample selector
                    selectInput("samples",
                            "Choose sample(s) that will be added to the plot.\n
                              The first one in the list will be used to order the plot.",
                                views,
                                multiple=TRUE),
                    ## The computation is quite extensive, plot on click only...
                    actionButton("goButton","Plot coverages", icon = icon("bar-chart-o")),
                    p(),
                    actionButton('zoom','Zoom',icon=icon('search-plus')),
                    actionButton('resetZoom','Reset Zoom',icon=icon('refresh')),
                    actionButton('resetMarker','Reset Marker',icon=icon('refresh')),
                    uiOutput("refreshPlot"),
                uiOutput("flag")
                    )
            })
        }
    })
    
    ## Register the yclicks for region selection
    y.clicks <- reactiveValues(y1=NULL,y2=NULL)
    observe({
        yvalue <- input$coords$y
        if (!is.null(yvalue)){
            isolate({
                if (!is.null(y.clicks$y2) | is.null(y.clicks$y1)){
                    y.clicks$y1 <- yvalue
                    y.clicks$y2 <- NULL
                } else {
                    y.clicks$y2 <- yvalue
                }
            })
        }
    })
    
    ## Figure out the rows from the original data that makes the selection
    selected <- reactive({
        ## make sure we react only when both clicks are registerd
        if(!is.null(y.clicks$y2)){
            ## Grab the location of the region of interest
            y.vals <- sort(c(y.clicks$y1,y.clicks$y2))
            ## Making sure that the clicks stay within the boundaries of the plot
            if (y.vals[1] < 0) y.vals[1] <- 1e-6
            if (y.vals[2] > 1) y.vals[2] <- 1
            ## Grab either the original slice or the previous slice to keep zooming
            if (is.null(ranks$slice)){
                r <- ranks$r
            } else {
                r <- ranks$slice
            }
            ## Compute the region that need to be resize (applying a range conversion)
            range <- sort(ceiling(length(r) * (1-y.vals)))
            ## Sub-sesting the ranks
            return(r[seq(range[1],range[2])])
        } else {
            return(NULL)
        }
    })

    ## Which data are we ploting
    ## Triggers only upon click to the plotting button
    data <- reactive({
        ## Isolate on to plot click and changes in Rel/Abs
        input$goButton
        input$analysisType
        isolate ({
            ## Make sure the data are loaded before trying to compute
            if(
                !(
                    length(input$samples) == 0 |
                    is.null(input$analysisType) |
                    is.null(Relative) |
                    is.null(Absolute)
                    )
                ){
                samples <- names(views)[match(input$samples,views)]
                d <- switch(input$analysisType,
                            Relative = Relative,
                            Absolute = Absolute)
                return(d[samples])
            } else {
                return(NULL)
            }
        })
    })

    ## Create a reactive that will keep tracks of ranks
    ranks <- reactiveValues()
    ## Set the $r slot to the initial ranks
    observe({
        if (!is.null(data()) ){
            ranks$r <- orderRank(data()[[1]])
        } else {
            ranks$r <- NULL
        }
    })
        
    ## Returning which slices is active
    slice <- reactive({
        if(is.null(ranks$r)){
            return(NULL)
        } else if(is.null(ranks$slice)){
            return(ranks$r)
        } else {
            return(ranks$slice)
        }
    })
    
    ## Acting on clicks to the zoom button
    ## The slot $slice is the current subset
    observe({
        input$zoom
        isolate({
            ranks$slice <- selected()
            ## Removing the blue box
            y.clicks$y1 <- y.clicks$y2 <- NULL
        })
    })
    
    ## Acting on clicks to the reset button
    observe({
        input$resetZoom
        isolate({
            ranks$slice <- NULL
            y.clicks$y1 <- y.clicks$y2 <- NULL
        })
    })
    
    ## Acting on clicks to the reset marker
    observe({
        input$resetMarker
        isolate({
            y.clicks$y1 <- y.clicks$y2 <- NULL
        })
    })
    
    ## Preparing the data needing to be ploted
    ## Need to play isolation... as it keeps getting recomputed every time a sample changes
    toPlot <- reactive({
        s <- slice()
        d <- data()
        if(!is.null(s) & !is.null(d)){
            withProgress(session, {
                setProgress(message = 'Computing the data', detail = "Will take a few seconds")
                setProgress(value = 0.1 )
                ord.data <- lapply(d,function(d) d[slice(),,drop=FALSE])
                setProgress(value = 0.5 )
                data2plot <- lapply(ord.data,downSample)
                setProgress(value = 1 )
                return(data2plot)
            })
        } else {
            return(NULL)
        }
    })

    ## Render the TSS plots
    observe ({
        ## Make sure we have something to plot first
        if ( !is.null( toPlot() ) ){
            ## Render the plot
            output$plot <- renderPlot({
                plotCovs(toPlot(),
                         input$addCenterMarker,
                         yvals=c(y.clicks$y1,y.clicks$y2))
            })
        } 
    })

    ## Create a reactive content to display
    ## Ben's messages on what to do with the clicks act on
    observe({
        if (!is.null(toPlot())){
            if (is.null(y.clicks$y1)){
                output$text <- renderText('Click the plot to select your initial boundary region')
            }
            else if (is.null(y.clicks$y2)){
                output$text <- renderText('Click to finish highlighting your region of interest')
            }
            else {
                output$text <- renderText('The genes from region you have selected are now available!')
            }
        } else {
            output$text <- renderText(NULL)
            output$plot <- renderPlot(NULL)
        }
    })
    
})    
