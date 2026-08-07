#!/bin/bash
#SBATCH --job-name=clustalo_job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --output=clustalo_%j.out
#SBATCH --error=clustalo_%j.err

set -euo pipefail

# Move to the directory from which sbatch was submitted.
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Load and activate the environment.
module load miniforge3

# Conda activation in non-interactive SLURM shells may require conda.sh.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate myenv

# Input and output files.
INPUT="unique_seq.fasta"
OUTPUT_FASTA="unique_seq_aligned.fasta"
OUTPUT_MATRIX="identity_matrix.txt"

# Check dependencies and input before starting.
command -v clustalo >/dev/null 2>&1 || {
    echo "ERROR: clustalo is not available in environment: ${CONDA_DEFAULT_ENV:-unknown}" >&2
    exit 1
}

[[ -s "$INPUT" ]] || {
    echo "ERROR: Input file not found or empty: $(pwd)/$INPUT" >&2
    exit 1
}


clustalo \
    --infile="$INPUT" \
    --outfile="$OUTPUT_FASTA" \
    --threads="${SLURM_CPUS_PER_TASK:-1}" \
    --outfmt=fa \
    --distmat-out="$OUTPUT_MATRIX" \
    --percent-id \
    --full \
    --force \
    --verbose

echo "Alignment written to: $OUTPUT_FASTA"
echo "Identity matrix written to: $OUTPUT_MATRIX"
