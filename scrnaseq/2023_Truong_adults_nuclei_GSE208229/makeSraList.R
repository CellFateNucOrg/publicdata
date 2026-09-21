library(dplyr)

serverPath="/Volumes/meister.data"

workDir=paste0(serverPath,"/publicData/scrnaseq/2023_Truong_adults_nuclei_GSE208229")
setwd(workDir)

df<-read.csv(paste0(workDir,"/SraRunTable.csv"))

table(df$experitment_design)

df<-df[df$experitment_design=="normal",]
