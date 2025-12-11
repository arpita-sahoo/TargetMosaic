import gzip
import os
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO import write
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import tempfile
import shutil


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


def search_similar_sequences(input_fasta, faa_file, k=5, threshold=0.9):
    """Search a single FAA file, return results instead of writing to disk."""
    query = read_fasta(input_fasta)
    query_seq = str(query.seq)

    sam_id = os.path.basename(faa_file).split(".")[0]

    matches = []
    count = 0

    for record in parse_faa(faa_file):
        target_seq = str(record.seq)
        sim = kmer_similarity(query_seq, target_seq, k=k)
        if sim >= threshold:
            new_id = f"{sam_id}_{record.id}|similarity={sim:.2f}"
            matches.append(SeqRecord(record.seq, id=new_id, description=""))
            count += 1

    return sam_id, matches, count


def process_file(args):
    """Wrapper for parallel processing."""
    input_fasta, faa_file, k, threshold = args
    return search_similar_sequences(input_fasta, faa_file, k=k, threshold=threshold)


def search_directory(input_fasta, directory, output_fasta, output_counts, k=5, threshold=0.9):
    """Process all FAA files in parallel and write outputs safely."""
    faa_files = [
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.endswith(".faa") or f.endswith(".faa.gz")
    ]

    if not faa_files:
        print("No FAA files found in directory.")
        return

    tasks = [(input_fasta, faa_path, k, threshold) for faa_path in faa_files]

    num_workers = multiprocessing.cpu_count()
    print(f"Running in parallel with {num_workers} workers...")

    # Use temporary files to avoid race conditions
    temp_dir = tempfile.mkdtemp()

    results = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for sam_id, matches, count in executor.map(process_file, tasks):
            # Save temporary fasta
            if matches:
                temp_fasta = os.path.join(temp_dir, f"{sam_id}_matches.fasta")
                with open(temp_fasta, "w") as out_f:
                    write(matches, out_f, "fasta")
            # Save temporary counts
            temp_count = os.path.join(temp_dir, f"{sam_id}_counts.tsv")
            with open(temp_count, "w", newline="") as tsvfile:
                writer = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")
                writer.writerow([sam_id, count])

    # Merge temporary files into final outputs
    with open(output_fasta, "w") as out_f:
        for f in os.listdir(temp_dir):
            if f.endswith("_matches.fasta"):
                with open(os.path.join(temp_dir, f)) as tmp_f:
                    out_f.write(tmp_f.read())

    with open(output_counts, "w", newline="") as out_counts:
        for f in os.listdir(temp_dir):
            if f.endswith("_counts.tsv"):
                with open(os.path.join(temp_dir, f)) as tmp_f:
                    out_counts.write(tmp_f.read())

    # Clean up
    shutil.rmtree(temp_dir)
    print(f"Processing complete. Results written to {output_fasta} and {output_counts}.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="K-mer based protein similarity search")
    parser.add_argument("--input", required=True, help="Input protein FASTA file")
    parser.add_argument("--faa", help="Single FAA file (.faa or .faa.gz)")
    parser.add_argument("--dir", help="Directory containing multiple FAA files")
    parser.add_argument("--out_fasta", required=True, help="Output FASTA file for matches")
    parser.add_argument("--out_counts", required=True, help="Output Counts file for SAM IDs")
    parser.add_argument("--k", type=int, default=5, help="K-mer size (default=5)")
    parser.add_argument("--threshold", type=float, default=0.9, help="Similarity threshold (default=0.9)")

    args = parser.parse_args()

    if args.faa:
        sam_id, matches, count = search_similar_sequences(args.input, args.faa, k=args.k, threshold=args.threshold)
        # Write results
        if matches:
            with open(args.out_fasta, "w") as out_f:
                write(matches, out_f, "fasta")
        with open(args.out_counts, "w", newline="") as tsvfile:
            writer = csv.writer(tsvfile, delimiter="\t", lineterminator="\n")
            writer.writerow([sam_id, count])
        print(f"Processing complete. Results written to {args.out_fasta} and {args.out_counts}.")

    elif args.dir:
        search_directory(args.input, args.dir, args.out_fasta, args.out_counts, k=args.k, threshold=args.threshold)
    else:
        raise ValueError("You must provide either --faa (file) or --dir (directory)")
