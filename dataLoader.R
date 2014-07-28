##################################################
## load the data
##################################################
source("helpers.R")
library(BSgenome.Dmelanogaster.BDGP5.73)

##################################################
## Read the different ROI files (GRanges objects)
##################################################
ROIfiles <- list.files('data','_ROI.rds$',full=TRUE)
ROI <- readRDS('data/ROI.rds')


##################################################
## Read in the coverage files (RleLists of coverages)
##################################################
## dir <- 'data/cov'
## files <- list.files(dir,"_cov\\.rds$",full=TRUE)
## covs <- structure(files,names=sub("_cov.+","",basename(files)))

## cov.data <- mclapply(covs,readRDS,mc.cores=25,mc.preschedule=FALSE)
## names(cov.data) <- sub("_cov\\.rds","",basename(covs))

dir <- 'data/views'
files <- list.files(dir,"_view\\.rds$",full=TRUE)
views <- structure(files,names=sub("_view.rds$","",basename(files)))
Absolute <- mclapply(views,function(file,ROI){
        view <- readRDS(file)
        doAbsolute(view,ROI)
    },ROI=ROI,mc.cores=25,mc.preschedule=FALSE)
Relative <- mclapply(Absolute,doRelative,mc.cores=25,mc.preschedule=FALSE)
names(Absolute) <- names(Relative) <- names(views)


##################################################
## read gene information pulled from biomart (data.frame)
##################################################
ids <- readRDS("data/martIDs.rds")
linkout <- readLines("data/linkout.txt")


##################################################
## Computing the sequences under the ROI
##################################################
DNA <- getSeq(BSgenome.Dmelanogaster.BDGP5.73,ROI)
names(DNA) <- names(ROI)
