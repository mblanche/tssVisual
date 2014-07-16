library(GenomicFeatures)
library(Rsamtools)
library(simrR)
library(edgeR)
library(rtracklayer)

##################################################
## Loading the location of all transcripts
##################################################
txdb<-loadDb(file.path("/home/mab/genomic_analysis/Conaway/Theo/Dmel_txdb_ensembl_v73.txdb"))
tx <- exonsBy(txdb,'tx',use.names=TRUE)

##################################################
## Loading expression data from some wild type samples
##################################################
BFL.exp <- BamFileList(list.files("/home/mab/genomic_analysis/Conaway/Theo/ssrp_reAnalyse_Dec2013/bam",
                                  patt="Ns.+\\.bam$",full=TRUE))
counts <- BamsGeneCount(BFL.exp,tx,
                        lib.strand='anti',
                        nCores=25)

dge <- DGEList(assay(counts))

##################################################
### Define the subset of tx that are expressed
##################################################
expressed.tx.id <- rownames(dge$counts)[rowMeans(cpm(dge)) > 1]

##################################################
## find the promoter for all tx
##################################################
ROI <- unlist(range(tx[expressed.tx.id]))

##################################################
## Removing tx that are too close to the chromosome edge
##################################################
is.plus.within <- start(ROI[strand(ROI) == "+"]) - 500 > 0
is.minus.within <- end(ROI[strand(ROI)=='-'])+500 < seqlengths(seqinfo(ROI))[as.vector(seqnames(ROI[strand(ROI)=='-']))]
ROI <- c(ROI[strand(ROI) == "+"][is.plus.within],ROI[strand(ROI) == "-"][is.minus.within])

##################################################
## Center our Region Of Interset around the tss, 500 nt on either side
##################################################
ROI[strand(ROI) == '+'] <- resize(shift(ROI[strand(ROI) == '+'],-500),width=2*500,fix='start')
ROI[strand(ROI) == '-'] <- resize(shift(ROI[strand(ROI) == '-'],500),width=2*500,fix='end')

##################################################
## Save ROI as RDS for used in ShinyApps
##################################################
saveRDS(ROI,"data/ROI.rds")

##################################################
## load the coverage for the bam files of the chip data
##################################################
BFL.chip <- BamFileList(list.files("/home/mab/genomic_analysis/Conaway/Theo/ChipSeq_polII_SSRP/bam",
                                   patt="\\.bam$",full=TRUE))

covs <- bams2Covs(BFL.chip,lib.strand='none',nCores=25)

##################################################
### save coverages as bigwigs
##################################################
dir.create("./bigwig",FALSE,TRUE)
mclapply(names(covs),function(exp){
    export(as(covs[[exp]],"GRanges"),file.path("bigwig",paste0(exp,".bw")))
},mc.cores=25,mc.preschedule=TRUE)

##################################################
### save coverages as RDS (need to find the most
### Effective way to process data
##################################################
dir.create("./data/cov",FALSE,TRUE)
mclapply(names(covs),function(exp){
    saveRDS(covs[[exp]],file.path("data/cov",paste0(exp,"_cov.rds")))
},mc.cores=25,mc.preschedule=TRUE)


##################################################
### How fast do we read these
##################################################
files <- list.files("data/cov","_cov\\.rds$",full=TRUE)
system.time({
    covs <- mclapply(files,readRDS,mc.cores=25,mc.preschedule=FALSE)
})
names(covs) <- sub("_cov\\.rds","",basename(files))

##################################################
## Recover an RLEViewList of the coverage
## around the ROI for every bam files
##################################################
ROI.cov <- mclapply(covs,function(cov) Views(cov,as(ROI,'RangesList')[names(cov)]),mc.preschedule=FALSE,mc.cores=25)

##################################################
## Save each individual coverages in its own RDS
## Will see if this make sense in the later run...
##################################################
dir.create('data',FALSE,TRUE)
mclapply(names(ROI.cov),function(n) saveRDS(ROI.cov[[n]],file.path('data',paste0(n,"_covFeats.rds"))),
         mc.preschedule=FALSE,
         mc.cores=25)





names(covs)[[1]]

gr <- ROI[1]

t <- Views(covs[[1]],as(gr,'RangesList')[names(covs[[1]])])
t[elementLengths(t) != 0][[1]]


bw.f <- list.files("bigwig",names(covs)[[1]],full=TRUE)

cov <- import(bw.f)


              , which=gr, as=“RleList”)
