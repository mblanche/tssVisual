library(shinyIncubator)

source("helpers.R")

shinyServer(function(input, output, session) {

    ## Assign a reactive in function of the user email address
    ## Then have a selction of data dir from the db of analysis
    dir <- 'data/cov'

    ## Get the names of the data files
    covs <- reactive({
        files <- list.files(dir,"_cov\\.rds$",full=TRUE)
        structure(files,names=sub("_cov.+","",basename(files)))
    })
    
    ## Maybe it makes sense to preload all the data at once...
    data <- reactiveValues()
    ranks <- reactiveValues()

    ## Keep an eye on the coverages
    observe ({
        if(!is.null(covs())){
            if(length(covs()) != 0){
                withProgress(session, {
                    setProgress(message = "Loading coverage data, please be patient",
                                detail = "This may take a few moments...")
                    data$covs <- mclapply(covs(),readRDS,mc.cores=25,mc.preschedule=FALSE)
                    names(data$covs) <- sub("_cov\\.rds","",basename(covs()))
                    
                    data$ROI <- readRDS("data/ROI.rds")

                    data$views <- mclapply(data$covs,function(cov){
                        Views(cov,as(data$ROI,'RangesList')[names(cov)])
                    },mc.preschedule=FALSE,mc.cores=25)
                    
                    data$Absolute <- cov2matrix(data$views,data$ROI)
                    data$Relative <- doRelative(data$Absolute)
                })
            }
        }
    })
    
    ## Render a widget for selecting the sample to display
     observe ({
         if (!is.null(covs())){
             ## Render the abolute vs relative radio
             output$typeSelector <- renderUI({
                 radioButtons('analysisType','Select Coverage Type',c('Relative','Absolute'))
             })
             output$markerSelector <- renderUI({
                 checkboxInput('addTSSmarker','Display the TSS marker',TRUE)
             })
             
             ## output$rangeSelector <- renderUI({
             ##     sliderInput("range", "Region around features to display coverages:",
             ##                 min = -500, max = 500, value = c(-500,500),
             ##                 step=50,ticks=TRUE)
             ## })
             
             ## Render the sample selector
             output$sampleSelector <- renderUI({
                 samples <- covs()
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
                     actionButton('reset','Reset',icon=icon('refresh')))
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

    ## read the ids from GRanges file
    ids <- readRDS("data/martIDs.rds")

    
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
            range <- sort(round(length(r) * (1-y.vals)))
            
            #range[which(range < 1)] <- 1
            #range[which(range > length(r))] <- length(r)
            ## Sub-sesting the ranks
            return(r[seq(range[1],range[2])])
    
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
        yvalue2presence <- input$coords$y
        if(!is.null(y.clicks$y2) & !is.null(selected())){
            isolate({
                if (!is.null(y.clicks$y2)){
                    y.vals <- sort(c(y.clicks$y1,y.clicks$y2))
                    ## Making sure that the clicks stay within the boundaries of the plot
                    if (y.vals[1] < 0) y.vals[1] <- 0
                    if (y.vals[2] > 1) y.vals[2] <- 1
                    ##
                    l <- nrow(toPlot()[[1]])
                    range <- sort(round(l * (y.vals)))
                    range[range == 0] <- 1
                    
                    initial.df <- names(data$ROI)[selected()]
                    fb.order <- match(initial.df,ids$ensembl_transcript_id)
                    samples <- names(covs())[match(input$samples,covs())]
                    data.subset <- toPlot()
                    #if(length(nrow(data.subset[[1]])) > 1){
                    #y.floor <- y.vals[1]*nrow(data.subset[[1]]))
                    y.floor <-range[1]
                    y.ceiling <-range[2]
                    
                    output$text <- renderText({ paste(nrow(data.subset[[1]]),"IM HREE","floor",y.floor,"ceil",y.ceiling) })
                    output$metaplot <- renderPlot({
                        data.t <- metaPrepData(data.subset,(y.floor:y.ceiling))
                        ggplot(data=data.t,aes(x=x,y=y,group=Exp,colour=Exp)) +
                            geom_line() +
                                xlab("Position") +
                                    ylab("Coverage")
                    })
                    
                    output$table <- renderDataTable({
                        data.frame(Name=initial.df,'Gene ID'=ids$ensembl_gene_id[fb.order],'Gene Symbol'=ids$flybasename_gene[fb.order])
                        
                    }, options = list(iDisplayLength = 10))
                }
                else {
                    output$table <- renderDataTable({ NULL })
                }
            })
        }
        else {
            output$table <- renderDataTable({ NULL })
            output$metaplot <- renderPlot({ NULL })
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
        input$reset
        isolate({
            ranks$slice <- NULL
            y.clicks$y1 <- y.clicks$y2 <- NULL
        })
    })
    
    ## Preparing the data needing to be ploted
    toPlot <- reactive({
        samples <- names(covs())[match(input$samples,covs())]
        if (length(samples) != 0){
            d <- lapply(data[[input$analysisType]][samples],function(d) d[slice(),])
            d <- lapply(d,downSample)
        }
    })

    observe ({
        ## Do not react on changes to input$samples
        samples <- names(covs())[match(input$samples,covs())]
        if (length(samples) != 0){
            ranks$r <- orderRank((data[[input$analysisType]][samples[1]]))
        }
    })

    ## Render the TSS plots
    observe ({
        ## Do not reaction on changes of prep data
        isolate({
            data.plot <- toPlot()
        })
        ## Recompute on clicks
        input$goButton
        ## Render a plotUI with variable width, then plot
        if (!is.null(data.plot)){
            #isolate({
            if (is.null(y.clicks$y1)){
                #output$text <- renderText({ paste0('Click the plot to select your initial boundary region')})
            }
            else if (is.null(y.clicks$y2)){
                #output$text <- renderText({ paste0('Click to finish highlighting your region of interest')})
            }
            else {
                #output$text <- renderText({ paste0('The genes from region you have selected are now available!')})
            }

            samples <- names(covs())[match(input$samples,covs())]
      
            output$plot <- renderPlot({
                plotCovs(toPlot(),
                         input$addTSSmarker,
                         yval=c(y.clicks$y1,y.clicks$y2)
                         )                                     
            })
        }
    })
})
