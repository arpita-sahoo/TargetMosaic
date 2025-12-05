#!/bin/bash
#SBATCH --job-name=extract_kmer_counter
#SBATCH --output=extract_kmer_counter_%j.out
#SBATCH --error=extract_kmer_counter_%j.err
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --partition=standard
##SBATCH --mail-user=user@example.com
##SBATCH --mail-type=BEGIN,END 

module load miniforge3
# module load anaconda3
source activate myenv
python TargetMosaic/extract_protein_kmer.py --input WT.fasta --dir ./faa_files/ --out_fasta matches.fasta --out_counts counts.tsv --k 5 --threshold 0.75

echo "Extraction and counting finished!"
