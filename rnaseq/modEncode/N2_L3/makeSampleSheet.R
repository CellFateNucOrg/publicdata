library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(stringr)
library(dplyr)

workDir<-"/Volumes/meister.data/publicData/rnaseq/modEncode/N2_L3"
setwd(workDir)
dir.create(workDir,showWarnings=F, recursive=T)

# for fetch_ngs
ss<-read.delim(paste0(workDir,"/sampleinfo.txt"))
write.table(ss$Run,paste0(workDir,"/ids.txt"),col.names=F,row.names=F,quote=F)


# for rnaseq
ss1<-read.csv(paste0(workDir,"/samplesheet/samplesheet.csv"))

df<-data.frame(sample=ss1$library_name,
               fastq_1=ss1$fastq_1,
               fastq_2=ss1$fastq_2,
               strandedness="auto",
               group="N2")

df$fastq_2[is.na(df$fastq_2)]<-""
df$sample<-gsub("L3_ce0120_rw002","N2_L3-0",df$sample)

df
write.csv(df,paste0(workDir,"/samplesheet.csv"),row.names=F,
           quote=F)


