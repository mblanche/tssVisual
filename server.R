library(shinyIncubator)
#library(shiny)

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

    ## Keep an eye on the coverages
    observe ({
        if(!is.null(covs())){
            if(length(covs()) != 0){
                withProgress(session, {
                    setProgress(message = "Loading coverage data, please be patient",
                                detail = "This may take a few moments...")
                    #readCovs(covFeats.name())
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
                 actionButton("goButton","Plot coverages")
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
    
    observe({
        yvalue2presence <- input$coords$y
        if (!is.null(yvalue2presence)){
            isolate({
                if (!is.null(y.clicks$y2)){
                    
                    highlight.y <- c(y.clicks$y1,y.clicks$y2)
                    highlight.y[which(highlight.y > 1)] <- 1
                    highlight.y[which(highlight.y < 0)] <- 0
                    highlight.max <- abs(1-min(highlight.y))
                    highlight.min <- abs(1-max(highlight.y))

                    #orderData <- names(readRDS("data/ROI.rds"))[orderRank((data[[input$analysisType]][samples[1]]))]
                    
                    num.genes <- as.numeric(length(names(data$ROI)[orderData()]))

                    initial.df <- names(data$ROI)[orderData()][floor(highlight.min*num.genes):ceiling(highlight.max*num.genes)]
                    fb.order <- match(initial.df,ids$ensembl_transcript_id)
                    
                    output$metaplot <- renderPlot({
                        all.data.t <- preppedData()
                        all.data.t <- metaPrepData(all.data.t,floor(min(highlight.y)*nrow(all.data.t[[1]])):ceiling(max(highlight.y)*nrow(all.data.t[[1]])))
                        ggplot(data=all.data.t,aes(x=x,y=y,group=Exp,colour=Exp)) +
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
    })

    ## order the data based on coverage at TSS
    orderData <- reactive({
        ## Do not react on changes to input$samples
        samples <- names(covs())[match(input$samples,covs())]
        if (length(samples) != 0){
            orderRank((data[[input$analysisType]][samples[1]]))
        }
    })

    preppedData <- reactive({
        ## Do not react on changes to input$samples
        samples <- names(covs())[match(input$samples,covs())]
        if (length(samples) != 0){
            prepData(data[[input$analysisType]][samples],orderData())
        }
    })
    
    ## Render the TSS plots
    observe ({
        ## Do not reaction on changes of prep data
        isolate({ data <- preppedData() })
        ## Recompute on clicks
        input$goButton
        ## Render a plotUI with variable width, then plot
        if (!is.null(data)){
            #isolate({
            if (is.null(y.clicks$y1)){
                output$text <- renderText({ paste0('Click the plot to select your initial boundary region')})
            }
            else if (is.null(y.clicks$y2)){
                output$text <- renderText({ paste0('Click to finish highlighting your region of interest')})
            }
            else {
                output$text <- renderText({ paste0('The genes from region you have selected are now available!')})
            }
            #})
            
            output$plot <- renderPlot({
                plotCovs(data,
                         input$addTSSmarker,
                         yval=c(y.clicks$y1,y.clicks$y2))
                                      
            })
        }
    })
})
