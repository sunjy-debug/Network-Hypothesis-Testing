#!/bin/bash
#SBATCH --job-name=mfmsbm-seq-nsbm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=Output/run_mfmsbm_seq/logs/%x-%j.out
#SBATCH --error=Output/run_mfmsbm_seq/logs/%x-%j.err

ml purge
ml fhR/4.4.1-foss-2023b-R-4.4.1

mkdir -p Output/run_mfmsbm_seq
mkdir -p Output/run_mfmsbm_seq/logs
unset DISPLAY

export OUTDIR="$PWD/Output/run_mfmsbm_seq"

Rscript --vanilla MFMSBM_sequential.R
Rscript --vanilla -e 'source("SBM_visualization.R"); SBM_visualization("_rslurm_mfmsbmseq", "Output/run_mfmsbm_seq")'

rm -rf _rslurm_mfmsbmseq