library(BSgenome.Celegans.UCSC.ce11)
library(rtracklayer)

#read in the public files' table, always set header so that the header is defined
publicfiles<-read.table("./publicfiles.txt", header=TRUE, stringsAsFactors = FALSE)
publicfiles
for (i in 1:dim(publicfiles)[1]) {
download.file(publicfiles$Sample[i],destfile=paste0(publicfiles$FileName[i],"_ce10.bw"))
}
#import the bigwig files
pub<-import("atac_wild_l3_ce10.bw")
pub<-import(paste0(publicfiles$FileName[i],"_ce10.bw"))
pub

#for the loops
for (i in 1:dim(publicfiles)[1]) {
pub<-import(paste0(publicfiles$FileName[i],"_ce10.bw"))


#import the downloaded chain file
chainFile<-import.chain("./ce10ToCe11.over.chain")

#lift the files from ce10 to ce11
pub_ce11<-unlist(liftOver(pub, chain = chainFile))
#to make sure all chromosomes are present in the seqlevels and seqlengths (get from BSgenome)
seqlevels(pub_ce11)<-seqlevels(Celegans)
seqlengths(pub_ce11)<-seqlengths(Celegans)


#to export the filess
export(pub_ce11, paste0(publicfiles$FileName[i],"_ce11.gtf"),  "gtf")
export(pub_ce11, paste0(publicfiles$FileName[i], "_ce11.bw"), "bw")
}
#to download from modENCODE
myftp=read.table("./modencodefiles.txt", header = TRUE, stringsAsFactors = FALSE)
#pathForFiles="/Users/imac/Desktop/Bolaji/PhD/PhD/Project_experiments/General/analysis/Bioinformatics_scripts/Bigwiggle"
#download.file(myftp,paste0(pathForFiles,basename(myftp)))
#download.file(myftp,paste0(basename(myftp)))
for (i in 1:dim(myftp)[1]) {
  download.file(myftp$Sample[i],destfile=paste0(myftp$FileName[i],".wig.gz"))
}
#to uncompress the files
system("gunzip ./*gz")

a<-wigToBigWig("./AMA-1.wig", seqinfo(Celegans), clip = TRUE)