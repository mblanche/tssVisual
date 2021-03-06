library(parallel)
library(GenomicRanges)
library(ggplot2)

##################################################
## Heler functions
##################################################

doAbsolute <- function(view,ROI) {
    cov.list <- as.list(viewApply(view,as.vector))
    d <- t(do.call(cbind,cov.list[sapply(cov.list,length) > 1]))
    ## Flip the strand if ROI is on the negative strand
    is.neg <- as.vector(strand(ROI[unlist(lapply(cov,names))]) == '-')
    d[is.neg,] <- d[is.neg,ncol(d):1]
    return(d)
}


doRelative <- function(data) {
   data/rowSums(data)
}


orderRank <- function(data){
    center <- ceiling(ncol(data)/2)
    range <- ceiling(ncol(data)*0.05)
    rank <- order(rowMeans(data[,(center-range):(center+range)]),decreasing=TRUE)
}

downSample <- function(d) {
    ## Remove any rows with Na
    d <- d[!rowSums(is.na(d)) > 0,]
    
    ## If there is more than 400 row, downsample
    if (nrow(d) > 400){
        breaks <- cut(1:nrow(d),400,labels=FALSE)
        d <- t(sapply(split(d,breaks),function(t){
            colMeans(matrix(t,ncol=ncol(d)))
        }))
    }
    ## Downsampling the columns, there should always be more than 200???
    if (ncol(d) > 200){
        breaks <- cut(1:ncol(d),200,labels=FALSE)
        d <- sapply(split(t(d),breaks),function(t){
            ## Make sure we are returning a matrix (I case we only have one gene...
            rowMeans(matrix(t,nrow=nrow(d),byrow=TRUE))
        })
    }
    ## We have the special case of having only one gene selected, which is now
    ## A vector from the previous step, return a matrix
    if(class(d) != 'matrix') d <- matrix(d,ncol=length(d))
        
    ## the 0,0 is bottom left, I want to print from top left, reverse the rows
    d <- d[nrow(d):1,,drop=FALSE]
}

plotCovs <- function(d,withTSSmarker=TRUE,yvals=NULL) {
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
        if(!is.null(yvals)){
            max.y <- max(yvals)
            min.y <- min(yvals)
            if(length(yvals) == 1){
                rect(-500,as.numeric(min.y-0.003),500,(as.numeric(max.y+0.003)),col='#FF000060')
            }
            if(length(yvals) == 2){
                rect(-500,as.numeric(min.y),500,(as.numeric(max.y)),col='#0000FF60')
            }
        }
        if (withTSSmarker)
            abline(v=0.5,col='red',lwd=1.5,lty=2)
    }
}

metaPlot <- function(data.list){
        ## Downsample the number of column, applying average smoothing
        d <- lapply(data.list,function(d){
            breaks <- split(seq(ncol(d)),cut(seq(ncol(d)),breaks=100,labels=FALSE))
            sapply(breaks,function(b) rowMeans(d[,b]))
        })
        
        ## Massage the data for ggplot
        d.f <- data.frame(x=unlist(lapply(d,function(x) seq(ncol(x)))),
                          y=unlist(lapply(d,colMeans)),
                          exp=rep(names(d),sapply(d,ncol))
                          )
        
        p <- ggplot(d.f,aes(x,y,colour=exp)) 
        p+geom_line()+labs(x="Position",y="Coverage")
    }


addIGVLink <- function(gnModel){
    if(class(gnModel) == 'GRangesList'){
        gene.ids <- names(gnModel)
        gnModel <- unlist(range(gnModel))
        names(gnModel) <- gene.ids
    }
    IGVlinks <- sprintf('http://localhost:60151/goto?locus=%s:%s-%s',
                        seqnames(gnModel),
                        start(gnModel),
                        end(gnModel)
                        )
    hwrite("IGV",link = IGVlinks,target="toIGV", table = FALSE)
}


addlinkOut <- function(){

}
