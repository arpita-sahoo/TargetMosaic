import gzip
import os
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write


def parse_faa(file_path):
    """Parse a local FAA file, auto-detecting gzip or plain text."""
    with open(file_path, "rb") as test_f:
        start = test_f.read(2)

    if start == b"\x1f\x8b":  # gzip magic number
        with gzip.open(file_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield record
    else:
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield record


def read_fasta(file_path):
    """Read the input protein FASTA (assumes single sequence)."""
    records = list(SeqIO.parse(file_path, "fasta"))
    if len(records) != 1:
        raise ValueError("Input FASTA should contain exactly one sequence")
    return records[0]


def kmer_similarity(seq1, seq2, k=5):
    """Compute Jaccard similarity between two sequences using k-mers."""
    kmers1 = {str(seq1[i:i+k]) for i in range(len(seq1) - k + 1)}
    kmers2 = {str(seq2[i:i+k]) for i in range(len(seq2) - k + 1)}
    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)
    return intersection / union if union > 0 else 0.0


def search_similar_sequences(input_fasta, faa_file, output_fasta="matches.fasta", output_counts="counts.tsv", k=5, threshold=0.9):
    query = read_fasta(input_fasta)
    query_seq = str(query.seq)

    sam_id = os.path.basename(faa_file).split(".")[0]

    print(f"\nSearching in {faa_file} for sequences ≥{threshold*100:.1f}% similar to query {query.id}\n")

    matches = []
    for record in parse_faa(faa_file):
        target_seq = str(record.seq)
        sim = kmer_similarity(query_seq, target_seq, k=k)
        if sim >= threshold:
            new_id = f"{sam_id}_{record.id}|similarity={sim:.2f}"
            # print(f">{new_id}")
            # print(record.seq)
            matches.append(SeqRecord(record.seq, id=new_id, description=""))

    if matches:
        with open(output_fasta, "a") as out_f:
            write(matches, out_f, "fasta")
        with open(output_counts, 'a', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
            writer.writerow([sam_id, len(matches)])
    else:
        with open(output_counts, 'a', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t', lineterminator='\n')
            writer.writerow([sam_id, "0"])

def search_directory(input_fasta, directory, output_fasta, output_counts, k=5, threshold=0.9):
    for file in os.listdir(directory):
        if file.endswith(".faa") or file.endswith(".faa.gz"):
            faa_path = os.path.join(directory, file)
            search_similar_sequences(input_fasta, faa_path, output_fasta=output_fasta, output_counts=output_counts, k=k, threshold=threshold)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="K-mer based protein similarity search")
    parser.add_argument("--input", required=True, help="Input protein FASTA file")
    parser.add_argument("--faa", help="Single FAA file (.faa or .faa.gz)")
    parser.add_argument("--dir", help="Directory containing multiple FAA files")
    parser.add_argument("--out_fasta", required=True, help="Output FASTA file for matches(default=matches.fasta)")
    parser.add_argument("--out_counts", required=True, help="Output Counts file for SAM IDs (default=counts.tsv)")
    parser.add_argument("--k", type=int, default=5, help="K-mer size (default=5)")
    parser.add_argument("--threshold", type=float, default=0.9, help="Similarity threshold (default=0.9)")

    args = parser.parse_args()

    # Clear output file before writing
    open(args.out_fasta, "w").close()
    open(args.out_counts, "w").close()

    if args.faa:
        search_similar_sequences(args.input, args.faa, args.out_fasta, args.out_counts, k=args.k, threshold=args.threshold)
    elif args.dir:
        search_directory(args.input, args.dir, args.out_fasta, args.out_counts, k=args.k, threshold=args.threshold)
    else:
        raise ValueError("You must provide either --faa (file) or --dir (directory)")



"""
## usage
python extract_protein_kmer.py --input 	Tsx_WT.fasta --dir ./faa_test/ --out_fasta matches.fasta --out_counts counts.tsv --k 5 --threshold 0.75

"""