
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(readxl)
library(dplyr)

outPath="/Volumes/meister.data/publicData/genomes/liftover"
setwd(outPath)

ce10toCe11url<-"http://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz"
ce10toCe11<-"ce10Toce11.over.chain"
download.file(ce10toCe11url,paste0(outPath,"/",ce10toCe11,".gz"))
system(paste0("gunzip ",outPath,"/",ce10toCe11,".gz"))
file.remove(paste0(outPath,"/",ce10toCe11,".gz"))
