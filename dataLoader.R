##################################################
## load the data
##################################################
source("helpers.R")

##################################################
## Read the different ROI files (GRanges objects)
##################################################
ROIfiles <- list.files('data','_ROI.rds$',full=TRUE)
ROIs <- lapply(ROIfiles,readRDS)
names(ROIs) <- sub('_ROI.rds$',"",ROIfiles)

##################################################
## Read in the coverage files (RleLists of coverages)
##################################################
dir <- 'data/cov'
files <- list.files(dir,"_cov\\.rds$",full=TRUE)
covs <- structure(files,names=sub("_cov.+","",basename(files)))

data <- list()
data$covs <- mclapply(covs,readRDS,mc.cores=25,mc.preschedule=FALSE)
names(data$covs) <- sub("_cov\\.rds","",basename(covs))


##################################################
## read gene information pulled from biomart (data.frame)
##################################################
ids <- readRDS("data/martIDs.rds")
