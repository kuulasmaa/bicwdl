
###############################################################################
#
# use conda ngs environment
#

conda create --name bic_ngs
conda activate ngs

###############################################################################
#
# parameters
#

DIR_DATA=""
DIR_REF=/research/work/kuulasma/ref_data/rnor60
DIR_WRK=/research/work/kuulasma/data/Isoflurane

###############################################################################
#
# genome indexes
#

bowtie2-build rnor6contam.fa rnor6contam

wget -c ftp://ftp.ensembl.org/pub/release-95/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz
wget -c ftp://ftp.ensembl.org/pub/release-95/fasta/rattus_norvegicus/ncrna/Rattus_norvegicus.Rnor_6.0.ncrna.fa.gz
wget -c ftp://ftp.ensembl.org/pub/release-95/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.95.gtf.gz
wget -c ftp://ftp.ensembl.org/pub/release-95/gff3/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.95.gff3.gz

cat ${DIR_REF}/ensembl/Rattus_norvegicus.Rnor_6.0.cdna.fa.gz ${DIR_REF}/ensembl/Rattus_norvegicus.Rnor_6.0.ncrna.fa.gz > ${DIR_REF}/ensembl/Rattus_norvegicus.Rnor_6.0.cdna.ncrna.fa.gz
kallisto index -i ${DIR_REF}/ensembl/Rattus_norvegicus.Rnor_6.0.cdna.ncrna.fa.gz.idx ${DIR_REF}/ensembl/Rattus_norvegicus.Rnor_6.0.cdna.ncrna.fa.gz

###############################################################################
#
# create folders
#

mkdir fastqc_raw
mkdir fastqc_trim
mkdir fastqc_decon
mkdir raw
mkdir trim
mkdir decon
mkdir align
mkdir R

###############################################################################
#
# raw fastq download
#

wget --user btk --password hE7yvuqP64QJfPs8XaWo --recursive https://bioinfoshare.utu.fi/DataTransfer/190315_J00117_0083_AHYGG5BBXX/180082_Stenroos/
md5sum -c checksums.md5

###############################################################################
#
# raw fastq qc
#

files=$(ls ${DIR_WRK}/raw/*.fastq.gz)
fastqc --extract --threads 10 --nogroup --outdir ${DIR_WRK}/fastqc_raw/ ${files}

multiqc -f ${DIR_WRK}/fastqc_raw/

###############################################################################
#
# trim adapters and low quality bases 
#

DATA_TRIMMOMATIC_ADAPTERS=/research/work/kuulasma/bin/miniconda3/envs/ngs/share/trimmomatic-0.38-1/adapters/TruSeq3-SE.fa
files=$(ls ${DIR_WRK}/raw/*.fastq.gz)
#files=$(ls ${DIR_WRK}/raw/*_L004_*.fastq.gz)
#files=$(ls ${DIR_WRK}/raw/*_L005_*.fastq.gz)
for file in ${files} ; do
	trimmomatic SE ${file} ${DIR_WRK}/trim/$(basename ${file%.fastq.gz})_trim.fastq.gz -threads 10 -phred33 ILLUMINACLIP:$DATA_TRIMMOMATIC_ADAPTERS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35
done

###############################################################################
#
# trimmed fastq qc 
#

files=$(ls ${DIR_WRK}/trim/*.fastq.gz)
fastqc --extract --threads 10 --nogroup --outdir ${DIR_WRK}/fastqc_trim/ ${files}

multiqc -f ${DIR_WRK}/fastqc_trim/

###############################################################################
#
# decontamination: remove ribosomal and mitocondrional reads 
#

files=$(ls ${DIR_WRK}/trim/*.fastq.gz)
#files=$(ls ${DIR_WRK}/trim/*_L004_*.fastq.gz)
#files=$(ls ${DIR_WRK}/trim/*_L005_*.fastq.gz)
for file in ${files} ; do
	bowtie2 -q --end-to-end --sensitive --phred33 --no-unal --threads 10 \
		--un-gz ${DIR_WRK}/decon/$(basename ${file%.fastq.gz})_decon.fastq.gz \
		-x ${DIR_REF}/contam/rnor6contam \
		-U ${file} \
		-S ${DIR_WRK}/decon/$(basename ${file%.fastq.gz})_align.sam

	rm ${DIR_WRK}/decon/$(basename ${file%.fastq.gz})_align.sam
done

###############################################################################
#
# decontamined fastq qc 
#

files=$(ls ${DIR_WRK}/fastq_decon/*_decon.fastq.gz)
fastqc --extract --threads 10 --nogroup --outdir ${DIR_WRK}/fastq_decon/ ${files}

multiqc -f ${DIR_WRK}/fastqc_decon/

###############################################################################
#
# combine fastq files from different lanes 
#

#files=$(ls ${DIR_WRK}/decon/*_L004_R1_*.fastq.gz)
#for file in ${files} ; do
#	cat ${file} ${file/_L004_R1_/_L005_R1_} > ${file/_L004_R1_/_L004-L005_R1_}; 
#done

###############################################################################
#
# kallisto pseudo alignment and counting
#

files=$(ls ${DIR_WRK}/decon/*_L004_R1_*.fastq.gz)
for file in ${files} ; do
	kallisto quant --threads 10 -i ${DIR_REF}/ensembl/Rattus_norvegicus.Rnor_6.0.cdna.ncrna.fa.gz.idx --output ${DIR_WRK}/align/$(basename ${file} | cut -d'_' -f1-4) --single --bootstrap-samples 100 --rf-stranded --fragment-length 350 --sd 30 ${file} ${file/_L004_R1_/_L005_R1_} &> $(echo ${DIR_WRK}/align/$(basename ${file} | cut -d'_' -f1-4).log)
done
