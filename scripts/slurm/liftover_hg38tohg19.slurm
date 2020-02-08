#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --partition serial
#SBATCH --time 10:00:00
#SBATCH --no-requeue
#SBATCH --mail-type=END
#SBATCH --mail-user=teemu.kuulasmaa@uef.fi
#SBATCH --cpus-per-task=10
#SBATCH --mem=64000


module load htslib
module load openjdk/1.8.0_202-b08


DIR_BIN="/home/groups/bic/bin/"
DIR_REF="/home/groups/bic/ref_data/"

INPUT=$1
OUTPUT=$(basename $INPUT)
OUTPUT=${OUTPUT%%.*}

java -Xms16G  -jar ${DIR_BIN}/picard.jar LiftoverVcf \
  CHAIN=${DIR_REF}/Liftover_Chain_Files/hg38ToHg19.over.chain \
  R=${DIR_REF}/hg19/gatk_bundle/hg19/ucsc.hg19.fasta \
  I=${INPUT} \
  O=${OUTPUT}_hg19.vcf.gz \
  REJECT=${OUTPUT}_hg19-rejected.vcf.gz \
  MAX_RECORDS_IN_RAM=10000

tabix -f -p vcf ${OUTPUT}_hg19.vcf.gz
tabix -f -p vcf ${OUTPUT}_hg19-rejected.vcf.gz
