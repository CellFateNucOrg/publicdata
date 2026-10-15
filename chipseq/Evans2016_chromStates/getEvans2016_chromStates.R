library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)
library(readxl)
library(dplyr)

outPath="/Volumes/meister.data/publicData/chipseq/Evans2016_chromStates"
setwd(outPath)

liftoverPath="/Volumes/meister.data/publicData/genomes/liftover"


#####################-
## Chromatin states Evans et al. (2016) (Ahinger lab)------
#####################-
# https://www.pnas.org/content/pnas/113/45/E7020.full.pdf

# Dataset S1. Coordinates of EE and L3 chromatin states
# File of chromosome, start position, end position, and state number. Coordinates are in WS220 and follow BED conventions (start positions are in zero-based coordinates and end positions in one-based coordinates).


if(!file.exists(paste0(outPath,"/chromDomains_L3_Evans2016_ce11.bed"))){
  chromStatesURL<-"https://www.pnas.org/doi/suppl/10.1073/pnas.1608162113/suppl_file/pnas.1608162113.sd01.xlsx"
  download.file(chromStatesURL,paste0(outPath,"/",basename(chromStatesURL)))


  chrAstates<-readxl::read_excel(paste0(outPath,"/pnas.1608162113.sd01.xlsx"),sheet="L3 autosome states",col_names=c("chr","start","end","state"))
  chrXstates<-readxl::read_excel(paste0(outPath,"/pnas.1608162113.sd01.xlsx"),sheet="L3 chr X states",col_names=c("chr","start","end","state"))


  chrStates<-rbind(chrAstates,chrXstates)
  chrStatesGR<-GRanges(paste0("chr",chrStates$chr,":",
                              (chrStates$start+1),"-",
                              chrStates$end))
  chrStatesGR$score<-c(chrAstates$state,chrXstates$state)
  chainCe10toCe11<-import.chain(paste0(liftoverPath,"/ce10Toce11.over.chain"))
  chrStatesGR_ce11<-unlist(liftOver(chrStatesGR,chain=chainCe10toCe11))

  seqinfo(chrStatesGR_ce11)<-ce11seqinfo
  export(chrStatesGR_ce11,
         con=paste0(outPath,"/chromStates_L3_Evans2016_ce11.bed"),
         format="bed")

  # Dataset S2. Coordinates of EE and L3 domains
  # Excel file of chromosome, start position, end position of EE and L3 active domains, border regions, and regulated domains (each in separate tab). Additionally, border regions have strand information in column six to indicate if active domain is on the left (−) or on the right (+). Coordinates are in WS220 and follow BED conventions.
  chromDomainsURL<-"https://www.pnas.org/highwire/filestream/623778/field_highwire_adjunct_files/1/pnas.1608162113.sd02.xlsx"
  download.file(chromDomainsURL,paste0(outPath,"/",basename(chromDomainsURL)))

  l3active<-readxl::read_excel(paste0(outPath,"/",
                                      basename(chromDomainsURL)),
                               sheet="L3 active domains",
                               col_names=c("chr","start","end"))
  l3active$name<-"active"
  l3active$score<-"."
  l3active$strand<-"*"
  l3regulated<-readxl::read_excel(paste0(outPath,"/",
                                         basename(chromDomainsURL)),
                                  sheet="L3 regulated domains",
                                  col_names=c("chr","start","end"))
  l3regulated$name<-"regulated"
  l3regulated$score<-"."
  l3regulated$strand<-"*"
  l3borders<-readxl::read_excel(paste0(outPath,"/",basename(chromDomainsURL)),sheet="L3 borders",col_names=c("chr","start","end","name","score","strand"))
  l3borders$name<-"border"


  chrDomains<-rbind(l3active,l3regulated,l3borders)
  chrDomainsGR<-GRanges(paste0("chr",chrDomains$chr,":",(chrDomains$start+1),"-",
                               chrDomains$end,":",chrDomains$strand))
  mcols(chrDomainsGR)<-chrDomains[,c("name","score")]
  chrDomainsGR_ce11<-unlist(liftOver(chrDomainsGR,chain=chainCe10toCe11))

  seqinfo(chrDomainsGR_ce11)<-ce11seqinfo
  chrDomainsGR_ce11$score<-NULL
  export(chrDomainsGR_ce11,
         con=paste0(outPath,"/chromDomains_L3_Evans2016_ce11.bed"),
         format="bed")
  file.remove(paste0(outPath,"/",basename(chromDomainsURL)))
}


# ── 1. Define state colors (positions 1–20 match state numbers) ───────────────
stateClrs <- c("#fe0003","#f59745","#008100","#74943a",
               "#c4d69c","#05ff7f","#ceff65","#fd0082",
               "#ff70d0","#ffb6c4","#7f00ff","#1900fe",
               "#528dd4","#814006","#b8cde4","#808080",
               "#dcdac3","#c4bc97","#938b54","#141414")

# ── 2. Helper: hex → "R,G,B" string required by the BED itemRgb field ─────────
hex_to_rgb <- function(hex) {
  m <- col2rgb(hex)
  paste(m[1], m[2], m[3], sep = ",")
}

# ── 3. Read the file ───────────────────────────────────────────────────────────
# The file has a leading row-index column, so we read 7 columns
bed <- read.table(
  "chromStates_L3_Evans2016_ce11.bed",
  header    = FALSE,
  sep       = "\t",
  col.names = c("chrom", "start", "end", "name", "state", "strand")
)

# ── 4. Build BED9 output ───────────────────────────────────────────────────────
# BED9 columns: chrom, chromStart, chromEnd, name, score, strand,
#               thickStart, thickEnd, itemRgb
bed9 <- data.frame(
  chrom      = bed$chrom,
  chromStart = bed$start,
  chromEnd   = bed$end,
  name       = paste0("state", bed$state),   # human-readable label in IGV
  score      = 0,                             # must be 0–1000; not used for color
  strand     = ".",
  thickStart = bed$start,
  thickEnd   = bed$end,
  itemRgb    = sapply(stateClrs[bed$state], hex_to_rgb),
  stringsAsFactors = FALSE
)

# ── 5. Write with required track line ─────────────────────────────────────────
out_file   <- "chromStates_L3_Evans2016_ce11_RGB.bed"
track_line <- 'track name="ChromStates_L3" description="Chromatin States (Evans 2016 ce11)" itemRgb="On" visibility=1'

writeLines(track_line, con = out_file)
write.table(bed9, file = out_file, append = TRUE,
            sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

message("Written to: ", out_file)
