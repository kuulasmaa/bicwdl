#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --partition serial
#SBATCH --time 00:05:00
#SBATCH --no-requeue
#SBATCH --mail-type=END
#SBATCH --mail-user=teemu.kuulasmaa@uef.fi

echo $0 > output.txt
