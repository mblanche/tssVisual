library(shinyIncubator)

source("helpers.R")
source("dataLoader.R")

shinyServer(function(input, output, session) {
    ## Create a reactive that will keep tracks of ranks
    ranks <- reactiveValues()
    
    ## Render a widget for selecting the sample to display
    observe ({
        ## Render the abolute vs relative radio
        output$typeSelector <- renderUI({
            radioButtons('analysisType','Select Coverage Type',c('Relative','Absolute'))
        })
        output$markerSelector <- renderUI({
            checkboxInput('addTSSmarker','Display the TSS marker',TRUE)
        })
        ## Render the sample selector
        output$sampleSelector <- renderUI({
                                        #samples <- covs
            samples <- covs
            selectInput("samples",
                        "Choose sample(s) that will be added to the plot.\n
                              The first one in the list will be used to order the plot.",
                        samples,
                         multiple=TRUE)
        })
        ## The computation is quite extensive, plot on click only...
        output$plotButton <- renderUI({
            list(
                actionButton("goButton","Plot coverages", icon = icon("bar-chart-o")),
                p(),
                actionButton('zoom','Zoom',icon=icon('search-plus')),
                actionButton('resetZoom','Reset Zoom',icon=icon('refresh')),
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
    
    selected <- reactive({
        if(!is.null(y.clicks$y2)){
            ## Grab the location of the region of interest
            y.vals <- sort(c(y.clicks$y1,y.clicks$y2))
            ## Making sure that the clicks stay within the boundaries of the plot
            if (y.vals[1] < 0) y.vals[1] <- 0
            if (y.vals[2] > 1) y.vals[2] <- 1
            ## Grab either the original slice or the previous slice to keep zooming
            if (is.null(ranks$slice)){
                r <- ranks$r
            } else {
                r <- ranks$slice
            }
            ## Compute the region that need to be resize (applying a range conversion)
            range <- sort(ceiling(length(r) * (1-y.vals)))
            range[range==0] <- 1 ## if y.vals[1] == 0, set the range to 1
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
    
    observe({
        if(!is.null(selected())){
            ## Ok, I have a selection, I will render a UI to create a plot and a dataTable
            output$selectedGenes <- renderUI({
                list(plotOutput('metaplot'),
                     dataTableOutput('table'))
            })
            ## Then, I will populate with some data
            isolate({
                y.vals <- sort(c(y.clicks$y1,y.clicks$y2))
                ## Making sure that the clicks stay within the boundaries of the plot
                if (y.vals[1] < 0) y.vals[1] <- 0
                if (y.vals[2] > 1) y.vals[2] <- 1
                ## Compute which row to keep
                l <- nrow(toPlot()[[1]])
                range <- sort(ceiling(l * (y.vals)))
                range[range == 0] <- 1 ## if y.vals is 0, set the row to 1

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
                
                ## initial.df <- names(data$ROI)[selected()]
                ## fb.order <- match(initial.df,ids$ensembl_transcript_id)
                ## samples <- names(covs)[match(input$samples,covs)]
                ## data.subset <- toPlot()
                ##                         #if(length(nrow(data.subset[[1]])) > 1){
                ##                         #y.floor <- y.vals[1]*nrow(data.subset[[1]]))
                ## y.floor <-range[1]
                ## y.ceiling <-range[2]

                ##output$text <- renderText({ paste(nrow(data.subset[[1]]),"IM HREE","floor",y.floor,"ceil",y.ceiling) })
                
                initial.df <- names(data$ROI)[selected()]
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

    ## Acting on clicks to the zoom button
    observe({
        input$zoom
        isolate({
            ranks$slice <- selected()
            ## Removing the blue box
            y.clicks$y1 <- y.clicks$y2 <- NULL
            #return( slice() )
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
        #samples <- names(covs)[match(input$samples,covs)]
        samples <- names(covs)[match(input$samples,covs)]
        if (length(samples) != 0){
            d <- lapply(data[[input$analysisType]][samples],function(d) d[slice(),,drop=FALSE])
            d <- lapply(d,downSample)
        }
    })

    observe ({
        ## Do not react on changes to input$samples
        samples <- names(covs)[match(input$samples,covs)]
        if (length(samples) != 0){
            ranks$r <- orderRank((data[[input$analysisType]][samples[1]]))
        }
    })

    counter <- reactiveValues(i=0)
    ## Render the TSS plots
    observe ({
        ## Do not reaction on changes of prep data
        isolate({
            counter$i <- counter$i + 1
            data.plot <- toPlot()
        })
        ## Recompute on clicks
        input$goButton
        ## Render a plotUI with variable width, then plot
        if (!is.null(data.plot)){
            if (is.null(y.clicks$y1)){
                #output$text <- renderText({ paste0('Click the plot to select your initial boundary region')})
            }
            else if (is.null(y.clicks$y2)){
                #output$text <- renderText({ paste0('Click to finish highlighting your region of interest')})
            }
            else {
                #output$text <- renderText({ paste0('The genes from region you have selected are now available!')})
            }
            output$plot <- renderPlot({
                plotCovs(toPlot(),
                         input$addTSSmarker,
                         yval=c(y.clicks$y1,y.clicks$y2)
                         )                                     
            })
        }
    })

})    
    
