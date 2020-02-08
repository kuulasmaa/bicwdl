version 1.0

##
## FinnGen demo pipeline
##

workflow test_pipeline {
  File phenoFile = "LIBRARY_RED/finngen-mock-data-library/R1/phenotypes/FINNGEN_R1_phenotypes.phe"
  File Rscript = "SANDBOX_GREEN/extract_phenos.R"

  scatter (chr in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]) {
    call filtter {
      input: chrom=chr
    }
  }

  call mod_pheno {
    input: phenoFile = phenoFile,
    Rscript = Rscript
  }
  
  output {
  }
}

task filtter {
  String chrom
  File chromFile = "LIBRARY_RED/finngen-mock-data-library/R1/FINNGEN_R1_release_CHR_${chrom}.vcf.gz"
  String chromOut = basename(chromFile)

  command {
      #Make a index file for vcf
      bcftools index -t ${chromFile} -o ${chromOut}.tbi

      #Rename SNPs with chr_pos_ref_alt
      bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' ${chromFile} -Oz -o FINNGEN_R1_release_CHR_${chrom}_SNPID.vcf.gz

      #Extract a list of SNPs with INFO higher then 0.9
      bcftools query -f '%ID\t%INFO/INFO\n' FINNGEN_R1_release_CHR_${chrom}_SNPID.vcf.gz | awk '$2 > 0.9 {print $1}' > good_SNP_${chrom}.txt

      #Subset vcf to only include good SNPs
      bcftools view --force-samples -i ID=@good_SNP_${chrom}.txt FINNGEN_R1_release_CHR_${chrom}_SNPID.vcf.gz -Oz -o FINNGEN_R1_release_CHR_${chrom}_SNPID_selected.vcf.gz

      # Index filtered file
      bcftools index -t FINNGEN_R1_release_CHR_${chrom}_SNPID_selected.vcf.gz
      bcftools index -t FINNGEN_R1_release_CHR_${chrom}_SNPID.vcf.gz

      # Count variants stats
      bcftools index -s FINNGEN_R1_release_CHR_${chrom}_SNPID_selected.vcf.gz > FINNGEN_R1_release_CHR_${chrom}_SNPID_selected.stat
      bcftools index -s FINNGEN_R1_release_CHR_${chrom}_SNPID.vcf.gz > FINNGEN_R1_release_CHR_${chrom}_SNPID.stat
  }

  output {
    File out1="${chromOut}.tbi"
    File out2="FINNGEN_R1_release_CHR_${chrom}_SNPID.vcf.gz"
    File out3="FINNGEN_R1_release_CHR_${chrom}_SNPID_selected.vcf.gz"
    File out4="good_SNP_${chrom}.txt"
    File out5="FINNGEN_R1_release_CHR_${chrom}_SNPID_selected.stat"
    File out6="FINNGEN_R1_release_CHR_${chrom}_SNPID.stat"
  }

  runtime {
    docker: "eu.gcr.io/finngen-staging-containers/bioinformatics:1.0.1"
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 50 HDD"
    zones: "europe-west1-b"
  }
}

task mod_pheno {
  File phenoFile = "LIBRARY_RED/finngen-mock-data-library/R1/phenotypes/FINNGEN_R1_phenotypes.phe"
  File Rscript = "SANDBOX_GREEN/extract_phenos.R"

  command {
    Rscript ${Rscript} ${phenoFile} FINNGEN_R1_phenotypes_modified.phe
  }

  output {
    File out1="FINNGEN_R1_phenotypes_modified.phe"
  }

  runtime {
    docker: "eu.gcr.io/finngen-staging-containers/bioinformatics:1.0.1"
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 50 HDD"
    zones: "europe-west1-b"
  }
}

# List the files needed to be backuped
ls /finngen/pipeline/cromwell/workflows/test_pipeline/21ee5b5c-999d-40e4-81da-bd829d57ce5b/call-filtter/shard-*/FINNGEN_R1_release_CHR_*_SNPID_selected.stat

# Make a folder for download request
mkdir -p /finngen/red/tisipila/files_to_download_19-09-19

# Copy the files to be donwloaded into the folder
cp /finngen/pipeline/cromwell/workflows/test_pipeline/21ee5b5c-999d-40e4-81da-bd829d57ce5b/call-filtter/shard-*/FINNGEN_R1_release_CHR_*_SNPID_selected.stat /finngen/red/tisipila/files_to_download_19-09-19

# Pakage the files into single file
cd /finngen/red/tisipila/files_to_download_19-09-19
ls
tar czf tsipila_SNPstat_download_18-09-19.tar.gz FINNGEN_R1_release_CHR_*_SNPID_selected.stat

# Check if successfull
rm FINNGEN_R1_release_CHR_*_SNPID_selected.stat
ls
tar xzf tsipila_SNPstat_download_18-09-19.tar.gz
ls
