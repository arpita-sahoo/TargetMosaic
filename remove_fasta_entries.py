import sys
from Bio import SeqIO

if len(sys.argv) != 4:
    print("Usage: python3 remove_fasta_entries.py input.fasta exclude.tsv aligned.fasta")
    sys.exit(1)

input_fasta = sys.argv[1]
exclude_tsv = sys.argv[2]
output_fasta = sys.argv[3]


# Read IDs from first column of exclude.tsv
exclude_ids = set()

with open(exclude_tsv, "r") as f:
    for line in f:
        line = line.strip()

        if not line:
            continue

        exclude_ids.add(line.split("\t")[0])


# Read FASTA and keep only sequences not in exclude list
kept_records = []
removed = 0

for record in SeqIO.parse(input_fasta, "fasta"):

    if record.id in exclude_ids:
        removed += 1
    else:
        kept_records.append(record)


# Write filtered FASTA
SeqIO.write(kept_records, output_fasta, "fasta")

print(f"Excluded IDs loaded: {len(exclude_ids)}")
print(f"Sequences removed: {removed}")
print(f"Sequences retained: {len(kept_records)}")
print(f"Output written to: {output_fasta}")
