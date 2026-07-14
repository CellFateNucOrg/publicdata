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
dim(osc) #2718
write.table(osc,file=paste0(workDir,"/oscillating_Hendriks_2014.tsv"),sep="\t",
            row.names=F,quote=F)
rising<-df[df$class=="rising",]
write.table(rising,file=paste0(workDir,"/rising_Hendriks_2014.tsv"),sep="\t",
            row.names=F,quote=F)


## oscillating genes from Meeuse 2020 -------
## PMID: 32687264 DOI: 10.15252/msb.20209498
url<-"https://static-content.springer.com/esm/art%3A10.15252%2Fmsb.20209498/MediaObjects/44320_2020_BFMSB209498_MOESM3_ESM.xlsx"
download.file(url, basename(url))

df<-read_excel(basename(url))
colnames(df)<-c("wormbaseID","publicID","sequenceID","amplitude","phase","class")
table(df$class)

osc1<-df[df$class=="Osc",]
dim(osc1) #3739
write.table(osc1,file=paste0(workDir,"/oscillating_Meeuse_2020.tsv"),sep="\t",
            row.names=F,quote=F)


## oscillating genes from Latorre 2015 PMID: 25737279
url<-"https://genesdev.cshlp.org/content/suppl/2015/03/03/29.5.495.DC1/Supplemental_TableS7.xlsx"
# need to download manually
#download.file(url,basename(url))
df<-read_excel(basename(url),col_names=F)
colnames(df)<-"sequenceID"
dim(df) #3269
ids<-read.csv(paste0(serverPath,"/publicData/genomes/WS298/c_elegans.PRJNA13758.WS298.geneIDs"),header=F)
colnames(ids)<-c("speciesID","wormbaseID","publicID","sequenceID","status","gene_biotype")
head(ids)
df<-inner_join(df,ids[,c("wormbaseID","publicID","sequenceID")])
dim(df) #3273
write.table(df,file=paste0(workDir,"/oscillating_Latorre_2015.tsv"),sep="\t",
            row.names=F,quote=F)


## plot overlaps
allids<-unique(c(osc$wormbaseID,osc1$wormbaseID,df$wormbaseID))

idmat<-data.frame(Hendriks2014=(allids %in% osc$wormbaseID),
                Meeuse2020=(allids %in% osc1$wormbaseID),
                Latorre2015=(allids %in% df$wormbaseID))
colSums(idmat)
fit<-eulerr::euler(as.matrix(idmat))
p<-plot(fit ,quantities=T)
ggsave("oscillatingGeneSetOverlap.pdf",p)



allids<-unique(c(rising$wormbaseID,osc1$wormbaseID,df$wormbaseID))

idmat<-data.frame(Hendriks2014_rising=(allids %in% rising$wormbaseID),
                  Meeuse2020=(allids %in% osc1$wormbaseID),
                  Latorre2015=(allids %in% df$wormbaseID))
colSums(idmat)
fit<-eulerr::euler(as.matrix(idmat))
p<-plot(fit ,quantities=T)
p
ggsave("oscillatingGeneSetOverlap_rising.pdf",p)
