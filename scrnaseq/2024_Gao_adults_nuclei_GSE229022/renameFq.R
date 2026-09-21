library(stringr)
# need to rename fastq files to fit the cellrange convension of:
# prefix_S1_L001_R{1|2}_001.fastq.gz

serverPath="/Volumes/meister.data"
workDir=paste0(serverPath,"/publicData/scrnaseq/2024_Gao_adults_nuclei_GSE229022")
setwd(workDir)

fq<-unlist(list.files(paste0(workDir,"/fastq"), pattern=".fastq.gz$"))

df<-data.frame(oldname=fq,
              sampleid=str_split_i(fq,"_",1),
               runid=str_split_i(fq,"_",2),
               read=gsub("\\.fastq\\.gz","",str_split_i(fq,"_",3)))

df$prefix<-paste(df$sampleid,df$runid,sep="_")
df$newname<-paste0(df$prefix,"_S1_L001_R",df$read,"_001.fastq.gz")

for(i in 1:nrow(df)){
  file.rename(from=paste0(workDir,"/fastq/",df$oldname[i]),
              to=paste0(workDir,"/fastq/",df$newname[i]))
}

libdir=gsub("^\\/Volumes","\\/mnt",paste0(workDir,"/fastq"))

# make one file per sample
for(i in unique(df$sampleid)){
  ss<-df[df$sampleid==i,]
  libs<-data.frame(fastqs=libdir,
                   sample=unique(ss$prefix),
                   library_type="Gene Expression")
  write.table(libs,file=paste0(workDir,"/libraries_",ss$sampleid[1],".csv"),row.names=F,quote=F,sep=",")
}
