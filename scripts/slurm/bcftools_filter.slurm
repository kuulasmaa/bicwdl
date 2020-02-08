#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --partition serial
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=END
#SBATCH --mail-user=teemu.kuulasmaa@uef.fi

module load tabix
module load bcftools

for f in {1..22} X Y ; do
  bcftools view -Ou -i 'INFO/TYPED=1 | ( INFO/AC>=1 & INFO/INFO>=0.6 )' ${f}.vcf.gz | \
    bcftools annotate -Oz --set-id '%CHROM:%POS:%REF:%ALT' -o ${f}_ifilter.vcf.gz ${f}.vcf.gz;
  tabix -p vcf ${f}_ifilter.vcf.gz;
done;
