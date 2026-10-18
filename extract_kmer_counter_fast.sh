#!/bin/bash
#SBATCH --job-name=kmer_fast
#SBATCH --output=kmer_fast_%j.out
#SBATCH --error=kmer_fast_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=standard
set -euo pipefail
module load miniforge3
source activate TargetMosaic-env

python TargetMosaic/extract_protein_kmer_fast.py \
    --input WT.fasta \
    --dir ./faa_files/ \
    --k 5 \
    --thresholds 0.10 0.25 0.50 0.75 0.90 \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --output-dir .
