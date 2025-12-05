#!/bin/bash
#SBATCH --job-name=flatten_copy
#SBATCH --output=flatten_copy_%j.out
#SBATCH --error=flatten_copy_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
##SBATCH --mail-user=user@example.com
##SBATCH --mail-type=BEGIN,END 

DEST="../../faa_files"

mkdir -p "$DEST"

# Flatten everything into one folder, skip if file already exists
find . -type f -exec rsync --ignore-existing {} "$DEST" \;
