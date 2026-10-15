#!/usr/bin/env bash
# Build one CPM-normalised bigWig per tissue from the GSE157017 (Durham et al. 2021,
# C. elegans L2 sciATAC-seq) pooled Tn5 cut-site files (<tissue>.cuts.bed).
#
# Requires: bedtools, bedGraphToBigWig (UCSC), curl, tar, python3 (chrom sizes only)
# Usage:    [WITH_CONTROL=1] bash build_tissue_bigwigs.sh [workdir] [threads]
#           WITH_CONTROL=1 also fetches the 2.8 GB gDNA negative control (gdna_neg) and builds
#           gdna_neg.blacklist_filtered.cpm.bw. Existing outputs are skipped.
set -euo pipefail

WORKDIR=${1:-tissue_bw}
THREADS=${2:-4}
TISSUES="coelomocyte glia gonad hypodermis intestine muscle neuron pharynx"
[ "${WITH_CONTROL:-0}" = 1 ] && TISSUES="$TISSUES gdna_neg.blacklist_filtered"
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157017/suppl
mkdir -p "$WORKDIR"/{cuts,bigwig,tmp}
cd "$WORKDIR"

# 1. Download missing cut-site BEDs (streams the 2.7 GB tarball once; extracts only the files needed)
MISSING=""
for t in $TISSUES; do [ -s cuts/$t.cuts.bed ] || MISSING="$MISSING $t.cuts.bed"; done
if [ -n "$MISSING" ]; then
  curl -s $GEO/GSE157017_bed_files.tar.gz | tar xz -C cuts $MISSING
fi

# 2. ce11 / WBcel235 chromosome sizes (named chrI..chrX, chrM to match the BEDs / UCSC).
#    UCSC alternative: curl -s https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/ce11.chrom.sizes
if [ ! -s ce11.chrom.sizes ]; then
  curl -s -H "Content-Type: application/json" \
    "https://rest.ensembl.org/info/assembly/caenorhabditis_elegans" | python3 -c '
import sys, json
for r in json.load(sys.stdin)["top_level_region"]:
    if r["coord_system"] == "chromosome":
        n = "chrM" if r["name"] == "MtDNA" else "chr" + r["name"]
        print(n, r["length"], sep="\t")' > ce11.chrom.sizes
fi

# 3. BED (60 bp windows centred on Tn5 cuts) -> coverage bedGraph scaled to CPM -> bigWig
#    CPM scale = 1e6 / (number of cut records in the tissue file)
build_one () {
  t=$1
  [ -s bigwig/$t.cpm.bw ] && { echo "$t: exists, skipping"; return 0; }
  bed=cuts/$t.cuts.bed
  n=$(wc -l < "$bed")
  scale=$(awk -v n="$n" 'BEGIN{printf "%.10g", 1e6/n}')
  # genomecov on integer counts, then scale in awk (bedtools -scale rounds output to ~4 sig. digits)
  LC_ALL=C sort -k1,1 -k2,2n "$bed" | \
    bedtools genomecov -i - -g ce11.chrom.sizes -bg | \
    awk -v s="$scale" 'BEGIN{OFS="\t"} {printf "%s\t%d\t%d\t%.6g\n", $1, $2, $3, $4*s}' > tmp/$t.bedGraph
  LC_ALL=C sort -k1,1 -k2,2n tmp/$t.bedGraph > tmp/$t.sorted.bedGraph
  bedGraphToBigWig tmp/$t.sorted.bedGraph ce11.chrom.sizes bigwig/$t.cpm.bw
  rm tmp/$t.bedGraph tmp/$t.sorted.bedGraph
  echo "$t: $n cut records, scale=$scale"
}
export -f build_one
printf '%s\n' $TISSUES | xargs -P "$THREADS" -I{} bash -c 'build_one {}'
ls -lh bigwig/
