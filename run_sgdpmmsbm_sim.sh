#!/bin/bash
#SBATCH --job-name=sgdpmmsbm-sim-nsbm
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=Output/run_sgdpmmsbm_sim/logs/%x-%j.out
#SBATCH --error=Output/run_sgdpmmsbm_sim/logs/%x-%j.err

ml purge
ml fhR/4.4.1-foss-2023b-R-4.4.1

mkdir -p Output/run_sgdpmmsbm_sim/logs
unset DISPLAY

export OUTDIR="$PWD/Output/run_sgdpmmsbm_sim"

Rscript --vanilla SGDPMMSBM_simultaneous.R