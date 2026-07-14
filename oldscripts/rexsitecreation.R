library(gdata)
library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(rtracklayer)
library(ggbio)

# BiocManager::install(c("gdata","GenomicRanges","BSgenome.Celegans.UCSC.ce11",
# "rtracklayer","ggbio"))

#########################
# Albritton 2017 rex sites
#########################

# get rex sites from Albritton paper
Albritton_elife2017="https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMjM2NDUvZWxpZmUtMjM2NDUtc3VwcDEtdjEuemlw/elife-23645-supp1-v1.zip?_hash=8w773VoPpu6mx1PUQppgaFNrEniZQDyGrnfBDR6QDt8%3D"
download.file(Albritton_elife2017,destfile="./Albritton_elife2017.zip")
unzip("./Albritton_elife2017.zip",files="Tables/SupplementaryFile1D.xls",junkpaths=T)
rexSites<-read.xls("./SupplementaryFile1D.xls")

# change chr names to UCSC formats
rexSites$Chromosome<-gsub("CHROMOSOME_","chr",rexSites$Chromosome)

# make Granges
rexGR<-GRanges(seqnames=rexSites$Chromosome,
               ranges=IRanges(start=rexSites$Start, end=rexSites$End),
               strand="*", rexSites[,c("Rank","Previous.Name","strength.category")])

# get liftover chain file from UCSC
download.file("http://hgdownload.cse.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz",
              destfile="./ce10ToCe11.over.chain.gz")
system("gunzip ./ce10ToCe11.over.chain.gz")
chainFile<-import.chain("./ce10ToCe11.over.chain")

#lift rex sites from ce10 to ce11
rexGR_ce11<-unlist(liftOver(rexGR,chain=chainFile))
# you need to make sure all chromosomes are present in the seqlevels and seqlengths (get from BSgenome)
seqlevels(rexGR_ce11)<-seqlevels(Celegans)
seqlengths(rexGR_ce11)<-seqlengths(Celegans)

# export as gtf
rexGR_ce11$strength.category<-factor(rexGR_ce11$strength.category,levels=c("weak","intermediate","strong"))
export(rexGR_ce11, "rexsites_Albritton2017_ce11.gtf","gtf")

# export as bed
rexGR_ce11$name<-ifelse(is.na(rexGR_ce11$Previous.Name),
                        paste0("ErcanRex",rexGR_ce11$Rank),
                        paste0("ErcanRex",rexGR_ce11$Rank,"_",rexGR_ce11$Previous.Name))
rexGR_ce11$score<-as.numeric(rexGR_ce11$strength.category)
rexGR_ce11$score<-round(rexGR_ce11$score*1000/max(rexGR_ce11$score),0)
export(rexGR_ce11,"rexsites_Albritton2017_ce11.bed","bed")

# plot rex sites on chromosome layout
autoplot(rexGR_ce11,layout="karyogram",aes(color = strength.category,fill = strength.category)) +
  theme_alignment() +
  ggtitle("Strength of rex sites")

#file.remove("Albritton_elife2017.zip")
file.remove("SupplementaryFile1D.xls")


#########################
# Albritton 2017 rex motifs (genome wide)
#########################

# get rex sites from Albritton paper
Albritton_elife2017="https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMjM2NDUvZWxpZmUtMjM2NDUtc3VwcDEtdjEuemlw/elife-23645-supp1-v1.zip?_hash=8w773VoPpu6mx1PUQppgaFNrEniZQDyGrnfBDR6QDt8%3D"
download.file(Albritton_elife2017,destfile="./Albritton_elife2017.zip")
unzip("./Albritton_elife2017.zip",files="Tables/SupplementaryFile1E.xlsx",junkpaths=T)
rexMotifs<-read.xls("./SupplementaryFile1E.xlsx")

# take only strong motifs (score 7-10 which are enriched on X vs autosomes)
rexMotifs<-rexMotifs[rexMotifs$Score>=7,]
# 218 sites left

# change chr names to UCSC formats
rexMotifs$Chromosome<-gsub("CHROMOSOME_","chr",rexMotifs$Chromosome)
rexMotifs$Chromosome<-gsub("MtDNA","M",rexMotifs$Chromosome)

# make Granges
rexMotifGR<-GRanges(seqnames=rexMotifs$Chromosome,
               ranges=IRanges(start=rexMotifs$Start, end=rexMotifs$End),
               strand=rexMotifs$Strand)

rexMotifGR$name<-paste0(rexMotifs$Chromosome,"_",rexMotifs$Start)
rexMotifGR$score<-rexMotifs$Score
# get liftover chain file from UCSC
#download.file("http://hgdownload.cse.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz", destfile="./ce10ToCe11.over.chain.gz")
#system("gunzip ./ce10ToCe11.over.chain.gz")
chainFile<-import.chain("./ce10ToCe11.over.chain")

#lift rex sites from ce10 to ce11
rexMotifGR_ce11<-unlist(liftOver(rexMotifGR,chain=chainFile))
# you need to make sure all chromosomes are present in the seqlevels and seqlengths (get from BSgenome)
seqlevels(rexMotifGR_ce11)<-seqlevels(Celegans)
seqlengths(rexMotifGR_ce11)<-seqlengths(Celegans)

# export as bed
export.bed(rexMotifGR_ce11,"rexMotifs_Albritton2017_ce11.bed","bed")

# plot rex sites on chromosome layout
# autoplot(rexMotifGR_ce11,layout="karyogram",aes(color = strength.category,fill = strength.category)) +
#   theme_alignment() +
#   ggtitle("Strength of rex sites")

file.remove("Albritton_elife2017.zip")
file.remove("SupplementaryFile1E.xlsx")





#########################
# Crane 2015 rex sites
#########################

# using excell sheet created manually from paper
rex_Crane<-read.xls("rexSites_Crane2015.xlsx")
rexCraneGR<-GRanges(seqnames="chrX",
                    ranges=IRanges(start=rex_Crane$midpoint, end=rex_Crane$midpoint+1),
                    strand="*", rex_Crane[,c("Rank","Name","HOTcold", "boundaryOverlap",                                             "X.Anderson2019")])

rexCraneGR_ce11<-unlist(liftOver(rexCraneGR,chain=chainFile))
rexCraneGR_ce11$name<-rexCraneGR_ce11$Name
rexCraneGR_ce11$score<-ifelse(rexCraneGR_ce11$boundaryOverlap=="notAtBoundary",500,1000)

export(rexCraneGR_ce11,"rexsites_Crane2015_ce11.bed","bed")


#########################
# Anderson 2019 rex sites
#########################

rexAndersonGR_ce11<-rexCraneGR_ce11[rexCraneGR_ce11$X.Anderson2019=="yes"]
rexAndersonGR_ce11$score<-1000

export(rexAndersonGR_ce11,"rexsites_Anderson2019_ce11.bed","bed")

