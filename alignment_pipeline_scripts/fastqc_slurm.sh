#!/bin/bash
#SBATCH --job-name=FastQC
#SBATCH --ntasks=2
#SBATCH --mem=16000


fastqc $1 -o $2

