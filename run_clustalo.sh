#!/bin/bash
#SBATCH --job-name=clustalo_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --output=clustalo_%j.out
#SBATCH --error=clustalo_%j.err

# Activate your environment
module load miniforge3
source activate myenv

# Input & output
INPUT="unique_seq.fasta"
OUTPUT_FASTA="unique_seq_aligned.fasta"
OUTPUT_MATRIX="identity_matrix.txt"

# Run Clustal Omega
clustalo \
  -i $INPUT \
  -o $OUTPUT_FASTA \
  --threads=$SLURM_CPUS_PER_TASK \
  --outfmt=fa \
  --distmat-out=$OUTPUT_MATRIX \
  --percent-id \
  -v
