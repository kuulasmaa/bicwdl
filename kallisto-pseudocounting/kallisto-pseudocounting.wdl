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

import "tasks/kallisto.wdl" as kallisto

workflow KallistoPseudcounting {
    input {
    }

    call kallisto.Pseudocounting as Pseudocounting {
    }

    output {
        File abundanceHDF5s = Pseudocounting.abundanceHDF5
    }
}
