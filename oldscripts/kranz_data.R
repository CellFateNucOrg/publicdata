library(BSgenome.Celegans.UCSC.ce11)
library(rtracklayer)

if (! file.exists("ce10ToCe11.over.chain")) {
  # download the liftover chain file and unzip it
  download.file("http://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz",
		destfile="ce10ToCe11.over.chain.gz")
  system("gunzip ce10ToCe11.over.chain.gz")
}

#read in the public files' table, always set header so that the header is defined
publicfiles<-read.table("./kranz_files.txt", header=TRUE, stringsAsFactors = FALSE)

for (i in 1:dim(publicfiles)[1]) {
  download.file(publicfiles$Sample[i],destfile=paste0(publicfiles$FileName[i],"_ws220.wig.gz"))
  #unzip the file
  gunzipCMD=paste0("gunzip ",publicfiles$FileName[i],"_ws220.wig.gz")
  system(gunzipCMD)

  #import the wig file
  pub<-import.wig(paste0(publicfiles$FileName[i],"_ws220.wig"))
  # convert the ws220 names to ce10 names
  seqlevels(pub)<-gsub("CHROMOSOME_","chr",seqlevels(pub))
  seqlevels(pub)<-gsub("MtDNA","M",seqlevels(pub))
  #seqlevels(pub)<-seqlevels(Celegans)

  #import the downloaded chain file
  chainFile<-import.chain("ce10ToCe11.over.chain")
  #lift the files from ce10 to ce11
  pub_ce11<-unlist(liftOver(pub, chain = chainFile))
  #to make sure all chromosomes are present in the seqlevels and seqlengths (get from BSgenome)
  seqlevels(pub_ce11)<-seqlevels(Celegans)
  seqlengths(pub_ce11)<-seqlengths(Celegans)

  #to export the filess
  export(pub_ce11, paste0(publicfiles$FileName[i], "_ce11.bw"), "bw")
  file.remove(paste0(publicfiles$FileName[i],"_ws220.wig"))
}


avrSignal<-data.frame(sample=publicfiles$FileName, autosomeAvr=NA,XchrAvr=NA)
#bwList<-list.files(pattern=".*bw")
autosomes<-paste0("chr",c("I","II","III","IV","V"))

for (i in 1:dim(avrSignal)[1]) {
  bw<-import(paste0(publicfiles$FileName[i], "_ce11.bw"),format="bigwig")

  avrSignal$autosomeAvr[i]<-mean(bw$score[as.vector(seqnames(bw) %in% autosomes)])
  avrSignal$XchrAvr[i]<-mean(bw$score[as.vector(seqnames(bw) == "chrX")])
}


avrSignal$difference<-avrSignal$XchrAvr-avrSignal$autosomeAvr
avrSignal$exp2<-round(2^avrSignal$difference,3)

write.csv(avrSignal,"XvsAutosomeAvrSignal.csv",row.names=F)
