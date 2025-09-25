#!/bin/bash
#SBATCH --job-name=faa_download
#SBATCH --output=download_%j.out
#SBATCH --error=download_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16     # adjust for parallelism
#SBATCH --mem=32G
#SBATCH --partition=standard

module load anaconda3
source activate myenv

# Download directory in current folder
OUTDIR=./faa_files/
mkdir -p "$OUTDIR"

# Number of chunks = number of parallel jobs (should not exceed CPUs)
NCHUNKS=16

# Split SAMID.tsv into NCHUNKS roughly equal parts
split -n l/$NCHUNKS SAMID.tsv chunk_

# Loop over each chunk and launch download in background
for CHUNK in chunk_*; do
    echo "Launching download for $CHUNK"
    bakrep download --tsv "$CHUNK" --match filetype:faa --directory "$OUTDIR"/ &
done

# Wait for all background jobs to finish
wait
echo "All downloads finished!"
