library(rtracklayer)
library(GenomeInfoDb)

serverPath="/Volumes/meister.data"
workDir=paste0(serverPath,"/publicData/genomes/blacklist")
setwd(workDir)

# blacklist
blacklist_url<-"https://github.com/Boyle-Lab/Blacklist/raw/refs/heads/master/lists/ce11-blacklist.v2.bed.gz"
download.file(blacklist_url,destfile=basename(blacklist_url))

system(paste0("gunzip ",basename(blacklist_url)))

blackl<-import("ce11-blacklist.v2.bed")

seqlevelsStyle(blackl)<-"Ensembl"

seqlevels(blackl)<-c("I","II","III","IV","V","X")
seqlevels(blackl)<-c("I","II","III","IV","V","X","MtDNA")

export(sort(blackl),"WB235-blacklist.v2.bed")


