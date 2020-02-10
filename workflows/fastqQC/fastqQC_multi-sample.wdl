version 1.0

##
## 
## 
##
## Main requirements/expectations :
## 
##
## Description of inputs:
##
## ** Runtime **
## 
##
## ** Workflow options **
## 
##
## ** Primary inputs **
## 
##

import "../../tasks/cutadapt.wdl" as cutadapt
import "../../tasks/fastqc.wdl" as fastqc

workflow fastqQC {
    input {
        File sampleInfo
        String outputDirFq = "trim"
        String outputDirQc = "fastqc"
        String? adapterForward = "AGATCGGAAGAG"  # Illumina universal adapter
        String? adapterReverse = "AGATCGGAAGAG"  # Illumina universal adapter
        Array[String]+? contaminations
        Map[String, String] dockerImages = {
          "fastqc": "uefbic/bic_ngs",
          "cutadapt": "uefbic/bic_ngs"
        }

        # 1st column (0) has sample name,
        # 2nd column has path to single-end or paired-end file (.fq.gz)
        # 3nd column has path to the paired file (.fq.gz)
        # If single-end sequencing is used then only the 1st and 2nd columns are filled.
        # No diplicate sample names are allowed!
        Array[Array[String]] samples = read_tsv(sampleInfo)
    }

	# Per sample
	scatter (sample in samples) {

		String sampleName = sample[0]
        File read1
        File read2

        call fastqc.Fastqc as FastqcRead1 {
            input:
                seqFile = read1,
                outdirPath = outputDirQc + "/",
                dockerImage = dockerImages["fastqc"]
        }

        if (defined(read2)) {
            call fastqc.Fastqc as FastqcRead2 {
                input:
                    seqFile = select_first([read2]),
                    outdirPath = outputDirQc + "/",
                    dockerImage = dockerImages["fastqc"]
            }
            String read2outputPath = outputDirFq + "/" + "_R2_trim.fq.gz"
        }

        if (runAdapterClipping) {
            call cutadapt.Cutadapt as Cutadapt {
                input:
                    read1 = read1,
                    read2 = read2,
                    read1output = outputDirFq + "/" + readgroupName + "_R1_trim.fq.gz",
                    read2output = read2outputPath,
                    adapter = select_all([adapterForward]),
                    anywhere = select_first([contaminations, []]),
                    adapterRead2 = adapterReverseDefault,
                    anywhereRead2 = if defined(read2)
                        then select_first([contaminations, []])
                        else [],
                    reportPath = outputDirFq + "/" + readgroupName +  "_cutadapt_report.txt",
                    dockerImage = dockerImages["cutadapt"]
            }

            call fastqc.Fastqc as FastqcRead1After {
                input:
                    seqFile = Cutadapt.cutRead1,
                    outdirPath = outputDirQc + "/",
                    dockerImage = dockerImages["fastqc"]
            }

            if (defined(read2)) {
                call fastqc.Fastqc as FastqcRead2After {
                    input:
                        seqFile = select_first([Cutadapt.cutRead2]),
                        outdirPath = outputDirQc + "/",
                        dockerImage = dockerImages["fastqc"]
                }
            }
        }

    output {
        File qcRead1 = if runAdapterClipping
            then select_first([Cutadapt.cutRead1])
            else read1
        File? qcRead2 = if runAdapterClipping
            then Cutadapt.cutRead2
            else read2
        File read1htmlReport = FastqcRead1.htmlReport
        File read1reportZip = FastqcRead1.reportZip
        File? read2htmlReport = FastqcRead2.htmlReport
        File? read2reportZip = FastqcRead2.reportZip
        File? read1afterHtmlReport = FastqcRead1After.htmlReport
        File? read1afterReportZip = FastqcRead1After.reportZip
        File? read2afterHtmlReport = FastqcRead2After.htmlReport
        File? read2afterReportZip = FastqcRead2After.reportZip
        File? cutadaptReport = Cutadapt.report
        Array[File] reports = select_all([
            read1htmlReport,
            read1reportZip,
            read2htmlReport,
            read2reportZip,
            read1afterHtmlReport,
            read1afterReportZip,
            read2afterHtmlReport,
            read2afterReportZip,
            cutadaptReport
            ])
    }
    
    parameter_meta {
        read1: {description: "The first or single end fastq file to be run through cutadapt.", category: "required"}
        read2: {description: "An optional second end fastq file to be run through cutadapt.", category: "common"}
        outputDir: {description: "The directory to which the outputs will be written.", category: "common"}
        adapterForward: {description: "The adapter to be removed from the reads first or single end reads.", category: "common"}
        adapterReverse: {description: "The adapter to be removed from the reads second end reads.", category: "common"}
        contaminations: {description: "Contaminants/adapters to be removed from the reads.", category: "common"}
        readgroupName: {description: "The name of the readgroup.", category: "common"}
        dockerImages: {description: "The docker images used. Changing this may result in errors which the developers may choose not to address.",
                       category: "advanced"}
        runAdapterClipping: {description: "Whether or not adapters should be removed from the reads.", category: "advanced"}
    }
} # end workflow



