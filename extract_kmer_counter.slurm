#!/bin/bash
#SBATCH --job-name=extract_kmer_counter
#SBATCH --output=extract_kmer_counter_%j.out
#SBATCH --error=extract_kmer_counter_%j.err
#SBATCH --partition=standard
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
##SBATCH --mail-user=user@example.com
##SBATCH --mail-type=BEGIN,END 

module load miniforge3
# module load anaconda3
source activate myenv

# thresholds to test
thresholds=(0.10 0.25 0.50 0.75 0.90)

for t in "${thresholds[@]}"
do
    # remove decimal point for cleaner filenames (e.g., 0.75 -> 75)
    t_label=$(echo $t | awk '{printf "%d", $1*100}')

    echo "Running threshold $t"

    python TargetMosaic/extract_protein_kmer.py \
        --input WT.fasta \
        --dir ./faa_files/ \
        --out_fasta matches_${t_label}.fasta \
        --out_counts counts_${t_label}.tsv \
        --k 5 \
        --threshold $t

    echo "Finished threshold $t"
done

echo "All runs completed!"
echo "Extraction and counting finished!"
