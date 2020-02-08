#!/bin/bash

# #SBATCH --job-name=plink_array
#SBATCH --ntasks=1
#SBATCH --partition parallel
#SBATCH --time 00:05:00
#SBATCH --no-requeue
# #SBATCH --mail-type=END
# #SBATCH --mail-user=teemu.kuulasmaa@uef.fi
# #SBATCH --mem-per-cpu=2500
#SBATCH --array=1-24

if [ $SLURM_ARRAY_TASK_ID -gt 0 ] & [ $SLURM_ARRAY_TASK_ID -lt 23 ]; then
  ARG="$SLURM_ARRAY_TASK_ID"
elif  [ $SLURM_ARRAY_TASK_ID -eq 23 ]; then
  ARG="X"
elif  [ $SLURM_ARRAY_TASK_ID -eq 24 ]; then
  ARG="Y"
elif  [ $SLURM_ARRAY_TASK_ID -eq 25 ]; then
  ARG="XY"
elif  [ $SLURM_ARRAY_TASK_ID -eq 26 ]; then
  ARG="MT"
fi

#source $HOME/.bashrc
#conda active ngs

INPUT_VCF=$ARG.vcf.gz

$HOME/bin/plink1.9/plink --vcf $INPUT_VCF --freq --out ${INPUT_VCF%%.*}
