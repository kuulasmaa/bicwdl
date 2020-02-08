#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --partition serial
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=END
#SBATCH --mail-user=teemu.kuulasmaa@uef.fi
#SBATCH --cpus-per-task=10
#SBATCH --mem=32000

BIN_DIR="/home/users/kuulasma/bin/"

${BIN_DIR}/plink_1.9b6/plink --bfile NPH_excluded-variants_ref_AF_gnomAD_norm_alt --pca --out NPH_excluded-variants_ref_AF_gnomAD_norm_alt
