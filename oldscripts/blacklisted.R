library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(rtracklayer)

#' Dowload data from the web
#'
#' Uses a url to download a data file from the web and if necessary unzip it
#' then import its it as a data.frame
#' @param dataURL String with url for data file
#' @param destFile String with name of destination file on local directory. If left as
#' default (NULL), the file will be named based on the name in the url
#' @param header Parameter for reading in the data as a data.frame. (default=FALSE)
#' @return Data frame with the downloaded data
#' @export
downloadWebData<-function(dataUrl, destFile=NULL, header=FALSE){
  if (is.null(destFile)){
    print("no destFile")
    destFile<-basename(dataUrl)
    print(destFile)
  }
  try(download.file(dataUrl,destfile=paste0("./",destFile)),silent=T)
  if (grepl(".gz$",destFile)) {
    system(paste0("gunzip ./",destFile))
  }
  if (grepl(".zip$",destFile)) {
    unzip(destFile,junkpaths=T)
  }
  df<-read.delim(gsub(".gz$","",destFile),header=header,stringsAsFactors=F)
  print("downloaded data:")
  print(head(df))
  return(df)
}

#' Convert data frame to genomic ranges object
#'
#' Converts data from a series of vectors to a sorted GRanges object
#' @param seqnames Vector with seqnames
#' @param start Vector with start positions
#' @param end Vector with end positions
#' @param strand Vector with strand specification
#' @param sortGR Sort genomic ranges (defualt=TRUE)
#' @return GenomicRanges object
#' @export
df2gr<-function(seqnames, start, end, strand,sortGR=TRUE) {
  gr<-GenomicRanges::GRanges(seqnames=seqnames,
          IRanges::IRanges(start=start,end=end),strand=strand)
  if (sortGR==TRUE) {
    gr<-sort(GenomeInfoDb::sortSeqlevels(gr))
  }
  print(gr)
  return(gr)
}

#modEncode blacklist https://sites.google.com/site/anshulkundaje/projects/blacklists
blackListed<-"https://github.com/Boyle-Lab/Blacklist/raw/master/lists/ce11-blacklist.v2.bed.gz"
black<-downloadWebData(blackListed)
blackGR<-sort(sortSeqlevels(import(gsub(".gz$","",basename(blackListed)))))
#blackGR<-with(blackdf, df2gr(V1,V2,V3,"*"))
hist(width(blackGR),breaks=100,main="Width of blacklisted regions")
abline(v=median(width(blackGR)),col="red")
file.remove(basename(blackListed))
export(blackGR,gsub(".gz$","",basename(blackListed)),"bed")


# repeats from UCSC
rmsk<-"http://hgdownload.soe.ucsc.edu/goldenPath/ce11/database/rmsk.txt.gz"
rmskdf<-downloadWebData(rmsk)
rmskGR<-with(rmskdf, df2gr(V6,V7+1,V8+1,V10))
rmskdf$V13<-ifelse(rmskdf$V12==rmskdf$V13,"",rmskdf$V13)
rmskGR$name<-with(rmskdf,paste(V11,V12,V13,sep="__"))
w<-width(rmskGR)
w[w>500]<-500
hist(w,breaks=100,main="Width of rmsk regions")
abline(v=median(width(rmskGR)),col="red")
export(rmskGR,"./ce11_repeatMasker.bed","bed")
file.remove((gsub(".gz$","",basename(rmsk))))


simpleRepeats<-"http://hgdownload.soe.ucsc.edu/goldenPath/ce11/database/simpleRepeat.txt.gz"
simple<-downloadWebData(simpleRepeats)
simpleGR<-with(simple, df2gr(V2,V3+1,V4+1,"*"))
hist(width(simpleGR),breaks=100,main="Width of simpleRepeat regions")
w<-width(simpleGR)
w[w>500]<-500
hist(w,breaks=100,main="Width of simpleRepeat regions")
abline(v=median(width(simpleGR)),col="red")
export(simpleGR,"./ce11_simpleRepeats.bed","bed")
file.remove((gsub(".gz$","",basename(simpleRepeats))))
#file.remove(basename(simpleRepeats))

nestedRepeats<-"http://hgdownload.soe.ucsc.edu/goldenPath/ce11/database/nestedRepeats.txt.gz"
nested<-downloadWebData(nestedRepeats)
nestedGR<-with(nested, df2gr(V2,V3+1,V4+1,"*"))
nested$V16<-ifelse(nested$V16==nested$V17,"",nested$V16)
nestedGR$name<-with(nested,paste(V5,V16,V17,sep="__"))
hist(width(nestedGR),breaks=100,main="Width of nestedRepeat regions")
w<-width(nestedGR)
w[w>5000]<-5000
hist(w,breaks=100,main="Width of nestedRepeat regions")
abline(v=median(width(nestedGR)),col="red")
export(nestedGR,"./ce11_nestedRepeats.bed","bed")
file.remove((gsub(".gz$","",basename(nestedRepeats))))


tableColNames<-"http://hgdownload.soe.ucsc.edu/goldenPath/ce11/database/tableList.txt.gz"