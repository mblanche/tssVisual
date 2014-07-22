library(shiny)

source("helpers.R")
source("dataLoader.R")

shinyServer(function(input, output, session) {
    ## Create a reactive that will keep tracks of ranks
    ranks <- reactiveValues()
    
    ## Render a widget for selecting the sample to display
    observe ({
        output$selectors <- renderUI({
            list(
                hr(),
                ## Render the abolute vs relative radio
                radioButtons('analysisType','Select Coverage Type',c('Relative','Absolute')),
                ## Render a button to add/remove TSS marker
                checkboxInput('addTSSmarker','Display the TSS marker',TRUE),
                hr(),
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
                actionButton('resetZoom','Reset Zoom',icon=icon('power-off')),
                actionButton('resetMarker','Reset Marker',icon=icon('refresh'))
                )
        })
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

    ## Returning which slices is active
    slice <- reactive({
        if(is.null(ranks$slice)){
            r <- ranks$r
        } else {
            r <- ranks$slice
        }
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
                
                initial.df <- names(data$ROI)[selected()]
                fb.order <- match(initial.df,ids$ensembl_transcript_id)
                IGV.links <- paste0("<a href=\"",pull_coords(initial.df,data$ROI),"\">Link</a>")
                output$table <- renderDataTable({
                    data.frame(Name=initial.df,
                               'Gene ID'=sapply(ids$ensembl_gene_id[fb.order],function(x) paste0("<a href=\"",linkOut,x,".html\">",x,"</a>")),
                               'Gene Symbol'=ids$flybasename_gene[fb.order],
                               'IGV'=IGV.links)
                }, options = list(iDisplayLength = 10))
                
            })
        } else {
            output$selectedGenes <- renderUI(NULL)
        }
    })

    ## Acting on clicks to the zoom button
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
        samples <- names(covs)[match(input$samples,covs)]
        if (length(samples) != 0){
            d <- lapply(data[[input$analysisType]][samples],function(d) d[slice(),,drop=FALSE])
            d <- lapply(d,downSample)
        } else {
            return(NULL)
        }
    })
    
    observe ({
        ## Do not react on changes to input$samples
        samples <- names(covs)[match(input$samples,covs)]
        if (length(samples) != 0){
            ranks$r <- orderRank((data[[input$analysisType]][samples[1]]))
        }
    })

    ## Render the TSS plots
    observe ({
        ## Do not react on changes of prep data
        isolate({
        })
        ## Recompute on clicks
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
                         input$addTSSmarker,
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

})    
    
