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
PHENOS="NPH_code NPHShunt_code"
#PHENOS="NPH_code"

for p in $PHENOS ; do
        #
        # Typed genotypes
        #
        # logistic regression
#        ${BIN_DIR}/plink_1.9b6/plink --bfile /home/users/kuulasma/data/neuro/EADB/imputation/NPH_excluded-variants_ref_AF_gnomAD_norm_alt_snpsex --logistic hide-covar --out Typed_${p} \
#          --pheno /home/users/kuulasma/data/neuro/EADB/assoc/EADB-NPH_2019-10-13.phe --pheno-name ${p} --remove /home/users/kuulasma/data/neuro/EADB/assoc/NPH_excluded_samples.tsv \
#          --covar /home/users/kuulasma/data/neuro/EADB/imputation/NPH_excluded-variants_ref_AF_gnomAD_norm_alt.cov --covar-name PC1,PC2,PC3,PC4,PC5 \
#          --missing-phenotype -9 --maf 0.01 --allow-no-sex --threads 10

#        tail -n+2 Typed_${p}.assoc.logistic | sed 's/\s\+/\t/g' | sed 's/^\t*//' | sort -k2,2 > Typed_${p}.assoc.logistic.sorted

#        echo -e "chrom\tpos\tref\talt\tpval\tor\tmaf" > Typed_${p}.assoc.logistic.pw
#        tail -n+2 /home/users/kuulasma/data/neuro/EADB/imputation/NPH_excluded-variants_ref_AF_gnomAD_norm_alt.frq.sorted | join -t $'\t' -1 2 -2 2 Typed_${p}.assoc.logistic.sorted - | \
#          awk 'BEGIN {FS="\t"; OFS="\t"} {split($1,a,":"); print $2,a[2],a[3],a[4],$9,$7,$13}' >> Typed_${p}.assoc.logistic.pw

#        rm Typed_${p}.assoc.logistic.sorted

        #
        # Imputed genotypes
        #
        ${BIN_DIR}/plink_1.9b6/plink --bfile ../imputation/michigan-topmed/chrAll_NPH_ifilter --logistic hide-covar --out Imputed-michigan-topmed_${p} \
          --pheno /home/groups/hiltunen_lab/data/EADB/assoc/EADB-NPH_2019-10-13.phe --pheno-name ${p} --remove /home/groups/hiltunen_lab/data/EADB/assoc/NPH_excluded_samples.tsv \
          --covar /home/groups/hiltunen_lab/data/EADB/imputation/NPH_excluded-variants_ref_AF_gnomAD_norm_alt.cov --covar-name PC1,PC2,PC3,PC4,PC5 \
          --missing-phenotype -9 --maf 0.01 --allow-no-sex --threads 10

        tail -n+2 Imputed-michigan-topmed_${p}.assoc.logistic | sed 's/\s\+/\t/g' | sed 's/^\t*//' | sort -k2,2 > Imputed-michigan-topmed_${p}.assoc.logistic.sorted

        echo -e "chrom\tpos\tref\talt\tpval\tor\tmaf" > Imputed-michigan-topmed_${p}.assoc.logistic.pw
        tail -n+2 ../imputation/michigan-topmed/chrAll_NPH_ifilter.frq.sorted | join -t $'\t' -1 2 -2 2 Imputed-michigan-topmed_${p}.assoc.logistic.sorted - | \
          awk 'BEGIN {FS="\t"; OFS="\t"} {split($1,a,":"); print $2,a[2],a[3],a[4],$9,$7,$13}' >> Imputed-michigan-topmed_${p}.assoc.logistic.pw

        rm Imputed-michigan-topmed_${p}.assoc.logistic.sorted
done
