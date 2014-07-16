library(parallel)
library(GenomicRanges)
#library(biomaRt)
library(ggplot2)

##################################################
## load the data
##################################################

#mart <- useMart("ensembl", dataset="dmelanogaster_gene_ensembl")
#roi <- data$ROI
#ids <-getBM(attributes= c("ensembl_transcript_id","ensembl_gene_id","flybasename_gene"),filters= "ensembl_transcript_id", values=unique(names(readRDS("data/ROI.rds"))),mart)
#saveRDS(ids,file='./data/martIDs.rds')

readCovs <- function(covFeats.name){
    newCovs <- covFeats.name[!sub("_cov.+","",basename(covFeats.name)) %in% ls()]
    if (length(newCovs) == 0) return()
    for (file in newCovs){
        var.name <- sub("_cov.+","",basename(file))
        assign(var.name,readRDS(file),envir=globalenv())
    }
}

getData <- function(exps){
    d <- mclapply(exps,function(cov){
        d <- t(do.call(cbind,as.list(viewApply(get(cov),as.vector))))
    },mc.cores=25,mc.preschedule=FALSE)
    names(d) <- exps
    return(d)
}

cov2matrix <- function(views,ROI){
    d <- mclapply(views,function(cov){
        d <- t(do.call(cbind,as.list(viewApply(cov,as.vector))))
        ## Flip the strand if ROI is on the negative strand
        is.neg <- as.vector(strand(ROI[unlist(lapply(cov,names))]) == '-')
        d[is.neg,] <- d[is.neg,ncol(d):1]
        d
    },mc.cores=25,mc.preschedule=FALSE)
    names(d) <- names(views)
    return(d)
}


orderRank <- function(data){
    center <- ceiling(ncol(data[[1]])/2)
    range <- ceiling(ncol(data[[1]])*0.05)
    rank <- order(rowMeans(data[[1]][,(center-range):(center+range)]),decreasing=TRUE)
}

doRelative <- function(data) {
    lapply(data,function(d) d/rowSums(d))
}

downSample <- function(d) {
    ## Downsampling the rows
    breaks <- cut(1:nrow(d),400,labels=FALSE)
    d <- t(sapply(split(d,breaks),function(t){
        colMeans(matrix(t,ncol=ncol(d)))
    }))
    ## Downsampling the columns
    breaks <- cut(1:ncol(d),200,labels=FALSE)
    d <- sapply(split(t(d),breaks),function(t){
        rowMeans(matrix(t,nrow=nrow(d),byrow=TRUE))
    })
        ## the 0,0 is bottom left, I want to print from top left, reverse the rows
    d <- d[nrow(d):1,]
}


metaPrepData <- function(d, d.t.sub) {
    d <- data.frame(
        x= as.vector(sapply(d,function(x) 1:ncol(x))),
        y= as.vector(sapply(d, function(x) colMeans(x[d.t.sub,]))),
        Exp=rep(names(d), sapply(d,function(x) ncol(x[d.t.sub,])))
        )
}

prepData <- function(d, rank){
    #d <- orderData(d)
    d <- lapply(d,function(d.t) d.t[rank,])
    d <- lapply(d,downSample)
}

plotCovs <- function(d,withTSSmarker=TRUE,yval=NULL) {
    cols <- colorRampPalette(c('black','yellow'))(256)

    if(length(d) < 5){
        ## makes all the printing region 1 in wide
        layout(matrix(c(seq(d),rep(0,5-length(d))),ncol=5))
    } else {
        ## Unless, there is more than five
        layout(matrix(seq(d),ncol=length(d)))
    }
    
    mar <- par('mar')
    mar[3] <- mar[3]+2
    par(mar=mar)
    
    for (exp in names(d)) {
        image(t(d[[exp]]),
              col=cols,
              axes=FALSE,
              main=exp)
        axis(3,c(0,0.5,1),c('-500','TSS','500'))
        ## add a horizontal colored line at y-coordinates defined by click
        if(!is.null(yval)){
            max.y <- max(yval)
            min.y <- min(yval)
            if(length(yval) == 1){
                rect(-500,as.numeric(min.y-0.003),500,(as.numeric(max.y+0.003)),col='#FF000060')
            }
            if(length(yval) == 2){
                rect(-500,as.numeric(min.y),500,(as.numeric(max.y)),col='#0000FF60')
            }
        }
        if (withTSSmarker)
            abline(v=0.5,col='red',lwd=1.5,lty=2)
    }
}


imageTSS <- function(d,withTSSmarker=TRUE,yval) {
  
    cols <- colorRampPalette(c('black','yellow'))(256)
    
    if(length(d) < 5){
        ## makes all the printing region 1 in wide
        layout(matrix(c(seq(d),rep(0,5-length(d))),ncol=5))
    } else {
        ## Unless, there is more than five
        layout(matrix(seq(d),ncol=length(d)))
    }
    
    mar <- par('mar')
    mar[3] <- mar[3]+2
    par(mar=mar)
    
    for (exps in names(d)) {
        image(t(d[[exps]]),
                  col=cols,
              axes=FALSE,
              main=exps)
        axis(3,c(0,0.5,1),c('-500','TSS','500'))
        ## add a horizontal colored rectangle at y-coordinates defined by click
        max.y <- max(yval)
        min.y <- min(yval)

        if(!is.null(yval)){
            if(length(yval) == 1){
                rect(-500,as.numeric(min.y-0.003),500,(as.numeric(max.y+0.003)),col='#FF000060')
            }
            if(length(yval) == 2){
                rect(-500,as.numeric(min.y),500,(as.numeric(max.y)),col='#0000FF60')
            }
        }
        if (withTSSmarker)
            abline(v=0.5,col='red',lwd=1.5,lty=2)
    }
}
