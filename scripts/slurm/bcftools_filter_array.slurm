#!/bin/bash

#SBATCH --ntasks=1              # -n
#SBATCH --partition parallel
#SBATCH --time 15:00:00         # -t
#SBATCH --no-requeue
#SBATCH --mail-type=END
#SBATCH --mail-user=teemu.kuulasmaa@uef.fi
#SBATCH --cpus-per-task=5
#SBATCH --mem=16000
#SBATCH --array=1-24

module load tabix
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

bcftools view -Ou --threads 5 -i 'INFO/TYPED=1 | ( INFO/AC>=1 & INFO/INFO>=0.6 )' ${ARG}.vcf.gz | \
  bcftools annotate -Oz --threads 5 --set-id '%CHROM:%POS:%REF:%ALT' -o ${ARG}_ifilter.vcf.gz;
tabix -p vcf ${ARG}_ifilter.vcf.gz;

#~/bin/plink_1.9b6/plink --vcf ${ARG}_ifilter.vcf.gz --make-bed --out ${ARG}_ifilter;
