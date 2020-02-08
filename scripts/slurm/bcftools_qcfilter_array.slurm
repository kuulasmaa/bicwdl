#!/bin/bash

# #SBATCH --job-name=TEST
#SBATCH --ntasks=1              # -n
#SBATCH --cpus-per-task=10
#SBATCH --partition parallel
#SBATCH --time 02:00:00         # -t
#SBATCH --no-requeue
# #SBATCH --mail-type=END
# #SBATCH --mail-user=teemu.kuulasmaa@uef.fi
# #SBATCH --mem-per-cpu=2500
#SBATCH --array=1-23

module load bcftools

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

#INPUT_VCF=chr${ARG}.dose.vcf.gz
#$HOME/bin/plink1.9/plink --vcf $INPUT_VCF --freq --out ${INPUT_VCF%%.*}

bcftools view -Oz -i 'INFO/TYPED=1 | ( INFO/IMPUTED=1 & INFO/AF>0 & INFO/R2>0.3 )' -o chr${ARG}.dose.ifilter.vcf.gz --threads 10 chr${ARG}.dose.vcf.gz
tabix -p vcf chr${ARG}.dose.ifilter.vcf.gz
