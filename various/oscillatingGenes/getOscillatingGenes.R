library(readxl)


serverPath="/Volumes/meister.data"
#serverPath="X:"
workDir=paste0(serverPath,"/publicData/various/oscillatingGenes")

# oscillating genes from Hendriks et al.(2014) -----
## https://doi.org/10.1016/j.molcel.2013.12.013

url<-"https://ars.els-cdn.com/content/image/1-s2.0-S1097276513009039-mmc2.xlsx"
download.file(url,basename(url))

df<-read_excel(basename(url),skip=3)
colnames(df)<-c("wormbaseID","publicID","sequenceID","class","amplitude","phase")
table(df$class)
osc<-df[df$class=="osc",]
write.table(osc,file=paste0(workDir,"/oscillating_Hendriks_2014.tsv"),sep="\t",
            row.names=F,quote=F)
rising<-df[df$class=="rising",]
write.table(rising,file=paste0(workDir,"/rising_Hendriks_2014.tsv"),sep="\t",
            row.names=F,quote=F)
