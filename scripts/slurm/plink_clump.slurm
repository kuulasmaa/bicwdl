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
PHENO="NPH_code"

${BIN_DIR}/plink_1.9b6/plink --bfile ../imputation/sanger/NPH-2019-10-09.vcfs/all_ifilter --clump $PHENO.assoc --out $PHENO.assoc.clumped
