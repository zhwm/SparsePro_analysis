#!/bin/bash
#
#SBATCH --ntasks 40            # number of tasks
#SBATCH --mem-per-cpu 10G            # mem-per-cpuory pool per process
#SBATCH -o TRAIT_fastGWA.out    # STDOUT
#SBATCH -t 12:00:00            # time (D-HH:MM)

~/utils/gcta_1.93.2beta/gcta64 --mbfile geno_chrs_bgen_files.txt --fastGWA-lr --pheno TRAIT_INToutcome.txt --threads 40 --out TRAIT_fastgwalr --maf 0.001
