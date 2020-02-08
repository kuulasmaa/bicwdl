version 1.0

##
## Generic plink workflow
##

workflow plink {
  input {
    File genoFile
    String params
    File? phenoFile
    String? out
  }

  call run {
    input:
      genoFile = genoFile,
      phenoFile = phenoFile,
      params = params,
      out = out
  }

  output {
    Array[File] outFiles = run.outFiles
  }
}

task run {
  input {
    File genoFile
    String params
    File? phenoFile
    String? out
    String? outDir = select_first([outDir, "output"])
  }

  command {
    set -e
    mkdir -p ~{outDir}
    plink --bfile ~{genoFile} \
        ~{"--out " + outDir + "/" + out} \
        ~{"--pheno " + phenoFile} \
        ~{params}
  }

  output {
    Array[File] outFiles = glob(outDir + "/*")
  }

  runtime {
    docker: "uefbic/bic_gwas"
    cpu: 2
    memory: "8 GB"
  }
}
