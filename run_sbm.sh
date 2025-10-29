#!/bin/bash
#SBATCH --job-name=sbm-nsbm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=Output/run_sbm/logs/%x-%j.out
#SBATCH --error=Output/run_sbm/logs/%x-%j.err

ml purge
ml fhR/4.4.1-foss-2023b-R-4.4.1

mkdir -p Output/run_sbm
mkdir -p Output/run_sbm/logs
unset DISPLAY

export OUTDIR="$PWD/Output/run_sbm"

Rscript --vanilla Bayes_DPMMSBM.R
Rscript --vanilla -e 'source("Bayes_DPMMSBM_Visualization.R"); SBM_visualization("_rslurm_sbm", "Output/run_sbm")'

rm -rf _rslurm_sbm