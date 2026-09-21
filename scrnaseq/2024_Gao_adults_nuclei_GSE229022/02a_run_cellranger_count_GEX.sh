#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=2
#SBATCH --job-name=cr_gex

WORK_DIR=/mnt/meister.data/publicData/scrnaseq/2024_Gao_adults_nuclei_GSE229022
OUT_DIR=${WORK_DIR}/results/cellranger_v10.1_gex
CELLRANGER=/mnt/meister.data/sharedSoftware/cellranger-10.1.0/bin/cellranger
mkdir -p $OUT_DIR


ID="day1_N2_OP50_rep2_SRX19885406"
REFERENCE=/mnt/meister.data/publicData/genomes/cellrangerRef/WBcel298/
LIBRARIES=${WORK_DIR}/libraries_SRX19885406.csv

$CELLRANGER count  --id=$ID --transcriptome=$REFERENCE \
  --libraries=$LIBRARIES --create-bam=true --chemistry="ARC-v1" \
  --jobmode=slurm --maxjobs=10 \
  --mempercore=48 --output-dir=${OUT_DIR} 
