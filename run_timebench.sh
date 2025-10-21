#!/bin/bash
#SBATCH --job-name=timebench
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=Output/timebench/logs/%x-%j.out
#SBATCH --error=Output/timebench/logs/%x-%j.err

ml purge
ml fhR/4.4.1-foss-2023b-R-4.4.1

mkdir -p Output/timebench/logs
unset DISPLAY

export OUTDIR="$PWD/Output/timebench"

Rscript --vanilla time_consumption.R