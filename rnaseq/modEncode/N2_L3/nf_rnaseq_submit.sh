#!/bin/bash
#SBATCH --time=1-17:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

genomeVer=WS298
genomeDir=/mnt/meister.data/publicData/genomes
genomeFile=${genomeDir}/${genomeVer}/c_elegans.PRJNA13758.${genomeVer}.genomic.fa.gz
gtfFile=${genomeDir}/${genomeVer}/c_elegans.PRJNA13758.WS298.canonical_geneset.gtf

WORK_DIR=/mnt/meister.data/publicData/rnaseq/modEncode/N2_L3
CONFIG_FILE=/mnt/meister.data/nf-core/unibe_izb.config


nextflow run nf-core/rnaseq -profile singularity -r 3.26.0 --input ${WORK_DIR}/samplesheet.csv --multiqc_title multiqc_rnaseq --outdir $WORK_DIR -c $CONFIG_FILE --fasta $genomeFile --gtf $gtfFile 

