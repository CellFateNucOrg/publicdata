#!/bin/bash
#SBATCH --time=0-05:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"
export NXF_SYNTAX_PARSER=v1

WORK_DIR=/mnt/meister.data/publicData/scrnaseq/2024_Gao_adults_nuclei_GSE229022
CONFIG_FILE=/mnt/meister.data/nf-core/unibe_izb_apptainer.config


nextflow run nf-core/fetchngs -r 1.12.0  -profile apptainer --input ids.csv --outdir $WORK_DIR -c $CONFIG_FILE --download-method sratools
