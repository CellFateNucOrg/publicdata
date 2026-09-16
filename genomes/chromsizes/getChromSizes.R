library(rtracklayer)
library(GenomeInfoDb)

serverPath="/Volumes/meister.data"
workDir=paste0(serverPath,"/publicData/genomes/chromsizes")
setwd(workDir)

# chormsizes
chormsizes_url<-"https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.chrom.sizes"
download.file(chormsizes_url,destfile=basename(chormsizes_url))

system(paste0("gunzip ",basename(chormsizes_url)))

system("sed 's/^chr//' ce11.chrom.sizes | sed 's/^M/MtDNA/' > WBcel235.chrom.sizes")


