##################################################
## load the data
##################################################
source("helpers.R")

dir <- 'data/cov'
files <- list.files(dir,"_cov\\.rds$",full=TRUE)
covs <- structure(files,names=sub("_cov.+","",basename(files)))

data <- list()
data$covs <- mclapply(covs,readRDS,mc.cores=25,mc.preschedule=FALSE)

names(data$covs) <- sub("_cov\\.rds","",basename(covs))

data$ROI <- readRDS("data/ROI.rds")

data$views <- mclapply(data$covs,function(cov){
    Views(cov,as(data$ROI,'RangesList')[names(cov)])
},mc.preschedule=FALSE,mc.cores=25)

data$Absolute <- cov2matrix(data$views,data$ROI)
data$Relative <- doRelative(data$Absolute)


## read the ids from GRanges file
ids <- readRDS("data/martIDs.rds")
