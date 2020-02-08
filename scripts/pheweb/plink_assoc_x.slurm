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
#PHENOS="NPH_code NPHShunt_code"
PHENOS="NPH_code"

for p in $PHENOS ; do
        #
        # Typed genotypes
        #
        # logistic regression
        ${BIN_DIR}/plink_1.9b6/plink --bfile /home/users/kuulasma/data/neuro/EADB/imputation/NPH_excluded-variants_ref_AF_gnomAD_norm_alt_snpsex --logistic hide-covar --out Typed_${p} \
          --pheno /home/users/kuulasma/data/neuro/EADB/assoc/EADB-NPH_2019-10-13.phe --pheno-name ${p} --remove /home/users/kuulasma/data/neuro/EADB/assoc/NPH_excluded_samples.tsv \
          --covar /home/users/kuulasma/data/neuro/EADB/imputation/NPH_excluded-variants_ref_AF_gnomAD_norm_alt.cov --covar-name PC1,PC2,PC3,PC4,PC5 \
          --missing-phenotype -9 --maf 0.01 --allow-no-sex --threads 10

        tail -n+2 Typed_${p}.assoc.logistic | sed 's/\s\+/\t/g' | sed 's/^\t*//' | sort -k2,2 > Typed_${p}.assoc.logistic.sorted

        echo -e "chrom\tpos\tref\talt\tpval\tor\tmaf" > Typed_${p}.assoc.logistic.pw
        tail -n+2 /home/users/kuulasma/data/neuro/EADB/imputation/NPH_excluded-variants_ref_AF_gnomAD_norm_alt.frq.sorted | join -t $'\t' -1 2 -2 2 Typed_${p}.assoc.logistic.sorted - | \
          awk 'BEGIN {FS="\t"; OFS="\t"} {split($1,a,":"); print $2,a[2],a[3],a[4],$9,$7,$13}' >> Typed_${p}.assoc.logistic.pw

        rm Typed_${p}.assoc.logistic.sorted

        #
        # Imputed genotypes
        #
        ${BIN_DIR}/plink_1.9b6/plink --bfile ../imputation/michigan/NPH-2019-10-09-chr1-22X/chrAll.dose.ifilter.snpsex --logistic hide-covar --out Imputed-michigan_${p} \
          --pheno /home/users/kuulasma/data/neuro/EADB/assoc/EADB-NPH_2019-10-13.phe --pheno-name ${p} --remove /home/users/kuulasma/data/neuro/EADB/assoc/NPH_excluded_samples.tsv \
          --covar /home/users/kuulasma/data/neuro/EADB/imputation/NPH_excluded-variants_ref_AF_gnomAD_norm_alt.cov --covar-name PC1,PC2,PC3,PC4,PC5 \
          --missing-phenotype -9 --maf 0.01 --allow-no-sex --threads 10

        tail -n+2 Imputed-michigan_${p}.assoc.logistic | sed 's/\s\+/\t/g' | sed 's/^\t*//' | sort -k2,2 > Imputed-michigan_${p}.assoc.logistic.sorted

        echo -e "chrom\tpos\tref\talt\tpval\tor\tmaf" > Imputed-michigan_${p}.assoc.logistic.pw
        tail -n+2 ../imputation/michigan/NPH-2019-10-09-chr1-22X/chrAll.dose.ifilter.frq.sorted | join -t $'\t' -1 2 -2 2 Imputed-michigan_${p}.assoc.logistic.sorted - | \
          awk 'BEGIN {FS="\t"; OFS="\t"} {split($1,a,":"); print $2,a[2],a[3],a[4],$9,$7,$13}' >> Imputed-michigan_${p}.assoc.logistic.pw

        rm Imputed-michigan_${p}.assoc.logistic.sorted

        #
        # Typed genotypes
        #
#        ${BIN_DIR}/plink_1.9b6/plink --bfile ../imputation/NPH_excluded-variants_ref_AF_gnomAD_norm_alt --assoc counts --out Typed_${p} \
#          --pheno EADB-NPH_2019-10-13.phe --pheno-name ${p} --remove NPH_excluded_samples.tsv \
#          --covar NPH_excluded-variants_ref_AF_gnomAD_norm_alt.cov --covar-name SNPSEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
#          --missing-phenotype -9 --maf 0.01 --allow-no-sex --threads 10

#        echo -e "chrom\tpos\tref\talt\tpval\tor\tnum_cases\tnum_controls" > Typed_${p}.assoc.pw
#        tail -n+2 Typed_${p}.assoc | sed 's/\s\+/\t/g' | sed 's/^\t*//' | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $4, $7, $9, $10, $5, $6}' - >> Typed_${p}.assoc.pw
#
#        ${BIN_DIR}/plink_1.9b6/plink --bfile ../imputation/NPH_excluded-variants_ref_AF_gnomAD_norm_alt --clump Typed_${p}.assoc --out Typed_${p}.assoc
#
#        tail -n+2 Typed_${p}.assoc | sed 's/\s\+/\t/g' | sed 's/^\t*//' | sort -k2,2 > Typed_${p}.assoc.sorted
#        tail -n+2 Typed_${p}.assoc.clumped | sed 's/\s\+/\t/g' | sed 's/^\t*//' | sort -k3,3 > Typed_${p}.assoc.clumped.sorted
#        join -t $'\t' -1 2 -2 3 -o "1.1,1.3,1.4,1.7,1.9,1.10,1.5,1.6" Typed_${p}.assoc.sorted Typed_${p}.assoc.clumped.sorted > Typed_${p}.assoc.clumped.pw.tmp
#        echo -e "chrom\tpos\tref\talt\tpval\tor\tnum_cases\tnum_controls" > Typed_${p}.assoc.clumped.pw
#        sort -k 1,1n -k 2,2n -k3,3 Typed_${p}.assoc.clumped.pw.tmp >> Typed_${p}.assoc.clumped.pw
#        rm Typed_${p}.assoc.sorted Typed_${p}.assoc.clumped.sorted Typed_${p}.assoc.clumped.pw.tmp

        #
        # Imputed genotypes
        #
#        ${BIN_DIR}/plink_1.9b6/plink --bfile ../imputation/sanger/NPH-2019-10-09.vcfs/all_ifilter --assoc counts --out Imputed_${p} \
#          --pheno EADB-NPH_2019-10-13.phe --pheno-name ${p} --remove NPH_excluded_samples.tsv \
#          --covar NPH_excluded-variants_ref_AF_gnomAD_norm_alt.cov --covar-name SNPSEX,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
#          --missing-phenotype -9 --maf 0.01 --allow-no-sex --threads 10
#
#        echo -e "chrom\tpos\tref\talt\tpval\tor\tnum_cases\tnum_controls" > Imputed_${p}.assoc.pw
#        tail -n+2 Imputed_${p}.assoc | sed 's/\s\+/\t/g' | sed 's/^\t*//' | awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $3, $4, $7, $9, $10, $5, $6}' - >> Imputed_${p}.assoc.pw
#
#        ${BIN_DIR}/plink_1.9b6/plink --bfile ../imputation/sanger/NPH-2019-10-09.vcfs/all_ifilter --clump Imputed_${p}.assoc --out Imputed_${p}.assoc
#
#        tail -n+2 Imputed_${p}.assoc | sed 's/\s\+/\t/g' | sed 's/^\t*//' | sort -k2,2 > Imputed_${p}.assoc.sorted
#        tail -n+2 Imputed_${p}.assoc.clumped | sed 's/\s\+/\t/g' | sed 's/^\t*//' | sort -k3,3 > Imputed_${p}.assoc.clumped.sorted
#        join -t $'\t' -1 2 -2 3 -o "1.1,1.3,1.4,1.7,1.9,1.10,1.5,1.6" Imputed_${p}.assoc.sorted Imputed_${p}.assoc.clumped.sorted > Imputed_${p}.assoc.clumped.pw.tmp
#        echo -e "chrom\tpos\tref\talt\tpval\tor\tnum_cases\tnum_controls" > Imputed_${p}.assoc.clumped.pw
#        sort -k 1,1n -k 2,2n -k3,3 Imputed_${p}.assoc.clumped.pw.tmp >> Imputed_${p}.assoc.clumped.pw
#        rm Imputed_${p}.assoc.sorted Imputed_${p}.assoc.clumped.sorted Imputed_${p}.assoc.clumped.pw.tmp
done
