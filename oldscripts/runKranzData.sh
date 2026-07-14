#! /bin/bash

## Allocate resources
#SBATCH --time=0-12:00:00
#SBATCH --mem=32G
#SBATCH --partition=all
##SBATCH --tmp=32G

## job name
#SBATCH --job-name="kranzData"

source ${CONDA_ACTIVATE}  MC-HiC-env


Rscript ./kranz_data.R

