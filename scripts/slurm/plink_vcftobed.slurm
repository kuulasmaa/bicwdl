#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --partition serial
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=END
#SBATCH --mail-user=teemu.kuulasmaa@uef.fi
#SBATCH --cpus-per-task=10
#SBATCH --mem=64000


DIR_BIN="/home/users/kuulasma/bin/"

${DIR_BIN}/plink_1.9b6/plink --vcf chr${1}_NPH_ifilter_hg19.vcf.gz --make-bed --allow-extra-chr --out chr${1}_NPH_ifilter_hg19

