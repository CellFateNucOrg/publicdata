#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=2
#SBATCH --job-name=cr_gex

WORK_DIR=/mnt/meister.data/publicData/scrnaseq/2023_Truong_adults_nuclei_GSE208229
ID="P0_normal_SRX16241259"
OUT_DIR=${WORK_DIR}/results/${ID}
CELLRANGER=/mnt/meister.data/sharedSoftware/cellranger-10.1.0/bin/cellranger

REFERENCE=/mnt/meister.data/publicData/genomes/cellrangerRef/WBcel298/
LIBRARIES=${WORK_DIR}/libraries_SRX16241259.csv

$CELLRANGER count  --id=$ID --transcriptome=$REFERENCE \
  --libraries=$LIBRARIES --create-bam=true  \
  --jobmode=slurm --maxjobs=10 \
  --mempercore=48 --output-dir=${OUT_DIR} 
