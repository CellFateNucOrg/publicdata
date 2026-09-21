library(dplyr)

serverPath="/Volumes/meister.data"

workDir=paste0(serverPath,"/publicData/scrnaseq/2024_Gao_adults_nuclei_GSE229022")
setwd(workDir)

df<-read.csv(paste0(workDir,"/SraRunTable.csv"))
head(df)

df<-df[grepl("day1",df$AGE) & grepl("N2",df$strain),]

experiments<-read.csv(paste0(workDir,"/GSE229022_web.csv"),header=F)

experiments<-experiments |> filter(grepl("Day1 ",V2), grepl("OP50",V2), grepl("N2",V2 ))

head(df)
df<-df[df$Library.Name %in% experiments$V1,]

write.table(df$Run[1],"ids.csv",quote=F,row.names=F,col.names=F)
