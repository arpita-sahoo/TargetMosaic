import csv
import argparse
from Bio import SeqIO

def main(fasta_file, output_csv_file="unique_seq.csv"):
    # Load sequences to compare
    records = list(SeqIO.parse(fasta_file, "fasta"))

    # Write to seq_dict
    seq_dict = {record.id: str(record.seq) for record in records}

    # Count occurrences of each unique sequence
    count = {}
    for seq in seq_dict.values():
        count[seq] = count.get(seq, 0) + 1

    sorted_count = dict(sorted(count.items(), key=lambda x: x[1], reverse=True))

    # Save counts to CSV
    with open(output_csv_file, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(["Sequence", "Count"])
        w.writerows(sorted_count.items())

    print("Occurrences of each sequence in Descending order")
    print(sorted_count.values())

    total_seq = len(seq_dict)
    print("\nTotal Number of sequences entered:", total_seq)

    unique_seq = len(count)
    print("\nNumber of unique sequences:", unique_seq)


    # Output fasta file
    output_fasta = "unique_seq.fasta"

    # Open output FASTA file
    with open(output_fasta, "w") as fasta:
        var=1
        for sequence, count in sorted_count.items():
            # Write in FASTA format
            fasta.write(f">Var{var}_{count}\n{sequence}\n")
            var+=1
    
    print(f"FASTA file '{output_fasta}' created successfully!")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute unique sequences and output CSV.")
    parser.add_argument("fasta", help="Input FASTA file of sequences to compare")
    parser.add_argument("-o", "--output", default="unique_seq.csv", help="Output CSV file name")

    args = parser.parse_args()
    main(args.fasta, args.output)

"""
python unique_seq.py matches.fasta

"""