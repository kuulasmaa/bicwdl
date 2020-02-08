version 1.0

task Pseudocounting {
    input {
        Array[File]+ inputFqs
        String index
        String output_dir = "."
        Boolean? bias
        Int? bootstrap_samples
        Int? seed
        Boolean? plaintext
        Boolean? fusion
        Boolean? single
        Boolean? single_overhang
        Boolean? fr_stranded
        Boolean? rf_stranded
        Float? fragment_length
        Float? sd
        String? pseudobam
        String? genomebam
        String? gtf
        String? chromosomes

        Int threads = 4
        String memory = "8G"
        String dockerImage = select_first([dockerImage,"quay.io/biocontainers/star:2.7.3a--0"])
    }

    command {
        set -e
        kallisto quant \
            --index ~{index} \
            --output-dir ~{output_dir} \
            ~{true="--bias" false="" bias} \
            ~{"--bootstrap_samples " + bootstrap_samples} \
            ~{"--seed " + seed} \
            ~{true="--plaintext" false="" plaintext} \
            ~{true="--fusion" false="" fusion} \
            ~{true="--single" false="" single} \
            ~{true="--single-overhang" false="" single_overhang} \
            ~{true="--fr-stranded" false="" fr_stranded} \
            ~{true="--rf-stranded" false="" rf_stranded} \
            ~{"--fragment-length " + fragment_length} \
            ~{"--sd " + sd} \
            ~{"--threads " + threads} \
            ~{"--pseudobam " + pseudobam} \
            ~{"--genomebam " + genomebam} \
            ~{"--gtf " + gtf} \
            ~{"--chromosomes " + chromosomes} \
            ~{sep=" " inputFqs}
    }

    output {
        #File abundanceHDF5 = basename(select_first(inputFqs), ".fq.gz") + "counts"
		File abundanceHDF5 = "counts"
    }

    runtime {
        cpu: threads
        memory: memory
        docker: dockerImage
    }

    parameter_meta {
        inputFqs: {description: "Filename for the single-end or paired-end fastq.gz files.", category: "required"}
        index: {description: "Filename for the kallisto index to be used for quantification.", category: "required"}
        output_dir: {description: "Directory to write output to. Equals to kallisto's --output-dir parameter.", category: "optional"}
        bias: {description: "Perform sequence based bias correction.", category: "optional"}
        bootstrap_samples: {description: "Number of bootstrap samples (default: 0). Equals to kallisto's --bootstrap-samples parameter.", category: "optional"}
        seed: {description: "Seed for the bootstrap sampling (default: 42).", category: "optional"}
        plaintext: {description: "Output plaintext instead of HDF5.", category: "optional"}
        fusion: {description: "Search for fusions for Pizzly.", category: "optional"}
        single: {description: " Quantify single-end reads.", category: "optional"}
        single_overhang: {description: "Include reads where unobserved rest of fragment is predicted to lie outside a transcript. Equals to kallisto's --single-overhang parameter.", category: "optional"}
        fr_stranded: {description: "Strand specific reads, first read forward. Equals to kallisto's --fr-stranded parameter.", category: "optional"}
        rf_stranded: {description: "Strand specific reads, first read reverse. Equals to kallisto's --rf-strandedparameter.", category: "optional"}
        fragment_length: {description: "Estimated average fragment length. Equals to kallisto's --fragment-length parameter.", category: "optional"}
        sd: {description: "Estimated standard deviation of fragment length (default: -l, -s values are estimated from paired end data, but are required when using --single).", category: "optional"}
        threads: {description: "Number of threads to use (default: 1).", category: "optional"}
        pseudobam: {description: "Save pseudoalignments to transcriptome to BAM file.", category: "optional"}
        genomebam: {description: "Project pseudoalignments to genome sorted BAM file.", category: "optional"}
        gtf: {description: "GTF file for transcriptome information (required for --genomebam).", category: "optional"}
        chromosomes: {description: "Tab separated file with chromosome names and lengths (optional for --genomebam, but recommended).", category: "optional"}
    }
}
