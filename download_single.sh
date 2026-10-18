#!/bin/bash
#SBATCH --job-name=faa_download
#SBATCH --output=download_%j.out
#SBATCH --error=download_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=standard
##SBATCH --mail-user=user@example.com
##SBATCH --mail-type=BEGIN,END

module load miniforge3
# module load anaconda3
source activate myenv

# Stop script if a command fails
set -e

# Download directory in current folder
OUTDIR=./faa_files/
mkdir -p "$OUTDIR"

# Download all files using a single BakRep process
echo "Starting download: $(date)"

bakrep download \
    --tsv ../SAMID.tsv \
    --match filetype:faa \
    --directory "$OUTDIR"

echo "All downloads finished: $(date)"
