library(shiny)

source("helpers.R")
source("dataLoader.R")

shinyServer(function(input, output, session) {

    ROI <- function() ROIs$TSS
    
    ## ## Create a ROI reactive that will return the correct ROI
    ## ROI <- reactive({
    ##     if(length(ROIs) > 1 & !is.null(input$ROI)){
    ##         ## If we switch ROI, reset some of the values
    ##         ranks$slice <- NULL
    ##         ranks$r <- NULL
    ##         ##y.clicks <- NULL
    ##         return(ROIs[[input$ROI]])
    ##     } else {
    ##         return(ROIs[[1]])
    ##     }
    ## })
    
    computeView <- reactive({
        if (!is.null( ROI() )){
            mclapply(cov.data, function(cov){
                Views(cov,as(ROI(),'RangesList')[names(cov)])
            },mc.preschedule=FALSE,mc.cores=25)
        } else {
            return(NULL)
        }
    })
    
    Absolute <- reactive({
        if ( !is.null(computeView()) ){
            return(cov2matrix(computeView(),ROI() ))
        } else {
            return(NULL)
        }
    })
    
    Relative <- reactive({
        if(!is.null(Absolute())){
            doRelative(Absolute() )
        } else {
            return(NULL)
        }
    })
    
    ## Render a widget for selecting the sample to display
    observe ({
        if (!is.null(Relative())) {
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
                                covs,
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
            ## ## Display the different ROIs (if more than one)
            ## if(length(ROIs) > 1){
            ##     output$ROIselector <- renderUI({
            ##         list(
            ##             radioButtons('ROI','Select you Regions of Interest',names(ROIs),names(ROIs)[1]),
            ##             hr()
            ##             )
            ##     })
            ## }
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


    ## Create a reactive that will keep tracks of ranks
    ranks <- reactiveValues()


    ## Create a reative values to reset the ranks when the ROI changes
    observe ({
        input$ROI
        isolate({
            ranks <- NULL
        })
    })
    
    ## Set the $r slot to the initial ranks
    observe({
        if (!is.null(input$analysisType) & length(input$samples) != 0){
            d <- switch(input$analysisType,
                        Relative = Relative(),
                        Absolute = Absolute())
            sample2sort <- names(covs)[match(input$samples,covs)]
            ranks$r <- orderRank(d[[sample2sort]])
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
    toPlot <- reactive({
        if(!is.null(slice())){
                samples <- names(covs)[match(input$samples,covs)]
                d <- switch(input$analysisType,
                            Relative = Relative(),
                            Absolute = Absolute())
                d <- lapply(d[samples],function(d) d[slice(),,drop=FALSE])
                d <- lapply(d,downSample)
        } else {
            return(NULL)
        }
    })

    ## Render the TSS plots
    observe ({
        ## Recompute on clicks or on changes of ROIs
        input$goButton
        ## Render a plotUI with variable width, then plot
        isolate({
            if (!is.null( toPlot() )){
                ## Got data to plot, render the plotting UI
                output$plotRegion <- renderUI(list(
                textOutput('text2'),
                    h5(textOutput('text')),
                    plotOutput('plot',clickId="coords",height='800px')
                    ))
                ## Render the plot
                output$plot <- renderPlot({
                    plotCovs(toPlot(),
                             input$addCenterMarker,
                             yval=c(y.clicks$y1,y.clicks$y2)
                             )
                    
                })
            }
        })
    })

    ## Create a reactive content to display
    ## Ben's messages on what to do with the clicks act on
    observe({
        if (is.null(y.clicks$y1)){
            output$text <- renderText('Click the plot to select your initial boundary region')
        }
        else if (is.null(y.clicks$y2)){
            output$text <- renderText('Click to finish highlighting your region of interest')
        }
        else {
            output$text <- renderText('The genes from region you have selected are now available!')
        }
    })
    
    ## Create a reactive context to whipe out the ploting region when there is noghing to plot
    observe({
        if(is.null(toPlot()))
            output$plotRegion <- renderUI(return(NULL))
    })


    ## Rendering the content of the Gene Table panel
    ## I think this block as to be looked over to make sure we don't process it for no good reason...
    observe({
        #& !is.null(selected())
        if(!is.null(y.clicks$y2)){
            ## Ok, I have a selection, I will render a UI to create a plot and a dataTable
            output$selectedGenes <- renderUI({
                list(plotOutput('metaplot'),
                     dataTableOutput('table'))
            })
            ## Then, I will populate with some data
            isolate({
                y.vals <- sort(c(y.clicks$y1,y.clicks$y2))
                ## Making sure that the clicks stay within the boundaries of the plot
                if (y.vals[1] < 0) y.vals[1] <- 1e-6
                if (y.vals[2] > 1) y.vals[2] <- 1

                ## Compute which row to keep
                range <- ceiling(nrow(toPlot()[[1]]) * y.vals)

                ## Massage the data for ggplot
                d.f <- data.frame(x=unlist(lapply(toPlot(),function(x) seq(ncol(x)))),
                                  y=unlist(lapply(toPlot(),function(d) colMeans(d[range[1]:range[2],,drop=FALSE]))),
                                  Exp=rep(names(toPlot()),each=sapply(toPlot(),ncol)))

                ## Render the plot
                output$metaplot <- renderPlot({
                    p <- ggplot(d.f,aes(x,y,colour=Exp)) 
                    p <- p+geom_line()+labs(x="Position",y="Coverage")
                    print(p)
                })
                
                initial.df <- names(ROI())[selected()]
                fb.order <- match(initial.df,ids$ensembl_transcript_id)
                                
                output$table <- renderDataTable({
                    data.frame(Name=initial.df,
                               'Gene ID'=ids$ensembl_gene_id[fb.order],
                               'Gene Symbol'=ids$flybasename_gene[fb.order])
                }, options = list(iDisplayLength = 10))
                
            })
        } else {
            output$selectedGenes <- renderUI(NULL)
        }
    })


})    
    
