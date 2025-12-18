import argparse
import os
from collections import defaultdict
from Bio import SeqIO
from Bio.Align import PairwiseAligner


def find_domain_positions(aligned_fasta, domains_fasta, reference_id):
    """
    Map domain sequences (from reference) to alignment column positions.
    Returns {domain_name: (alignment_start, alignment_end)}.
    """
    records = {rec.id: rec for rec in SeqIO.parse(aligned_fasta, "fasta")}
    if reference_id not in records:
        raise ValueError(f"Reference ID '{reference_id}' not found in alignment.")
    ref_seq_aligned = str(records[reference_id].seq).upper()

    domain_positions = {}
    ungapped_ref = ref_seq_aligned.replace("-", "")

    # Configure aligner
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5

    for domain in SeqIO.parse(domains_fasta, "fasta"):
        domain_seq = str(domain.seq).upper()

        # Perform local alignment
        alignments = aligner.align(ungapped_ref, domain_seq)
        aln = alignments[0]  # best alignment

        # Extract start and end in ungapped reference
        start_in_ref = aln.aligned[0][0][0]  # first block start
        end_in_ref = aln.aligned[0][-1][1]   # last block end

        # Map ungapped ref positions → alignment indices
        mapping = {}
        ungapped_pos = 0
        for aln_idx, aa in enumerate(ref_seq_aligned, 1):  # 1-based
            if aa != "-":
                ungapped_pos += 1
                mapping[ungapped_pos] = aln_idx

        aln_start = mapping[start_in_ref + 1]  # convert to 1-based
        aln_end = mapping[end_in_ref]

        domain_positions[domain.id] = (aln_start, aln_end)

    return domain_positions


def extract_domains_from_alignment(aligned_fasta, domain_positions, fasta_out, csv_out):
    """
    Extract domains from all sequences in an aligned FASTA,
    group identical domain sequences, and merge headers.
    Sorts variants alphabetically within each header, 
    and final FASTA by total occurrences (descending).
    """
    grouped = defaultdict(list)

    # Parse aligned FASTA
    for record in SeqIO.parse(aligned_fasta, "fasta"):
        header = record.id
        seq = str(record.seq)

        # Extract all domain regions
        extracted = "".join(seq[start - 1:end] for _, (start, end) in domain_positions.items())
        grouped[extracted].append(header)

    results = []
    for seq, headers in grouped.items():
        variants = []
        total_occurrences = 0
        for h in headers:
            parts = h.split("_")
            variant = parts[0]
            count = int(parts[1]) if len(parts) > 1 else 1
            variants.append(variant)
            total_occurrences += count

        variants.sort()
        combined_header = "_".join(variants) + f"_{total_occurrences}"
        results.append((seq, combined_header, total_occurrences))

    results.sort(key=lambda x: x[2], reverse=True)

    with open(fasta_out, "w") as out:
        for seq, header, total_occurrences in results:
            out.write(f">{header}\n{seq}\n")
    with open(csv_out, "w") as out:
        for seq, header, total_occurrences in results:
            out.write(f"\t{total_occurrences}\t{seq}\t{header}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract domain sequences from an aligned FASTA.")
    parser.add_argument("--aligned", required=True, help="Aligned FASTA file (input)")
    parser.add_argument("--domains", required=True, help="FASTA file with reference domain sequences")
    parser.add_argument("--reference", required=True, help="Reference sequence ID in alignment")
    parser.add_argument("--fasta_out", help="Output FASTA file (optional)")
    parser.add_argument("--csv_out", help="Output CSV file (optional)")

    args = parser.parse_args()

    # Automatically derive suffix from domains file name
    domain_basename = os.path.splitext(os.path.basename(args.domains))[0]

    fasta_out = args.fasta_out or f"extracted_domains_{domain_basename}.fasta"
    csv_out = args.csv_out or f"extracted_domains_{domain_basename}.csv"

    # Step 1: find domain positions in alignment
    domain_positions = find_domain_positions(args.aligned, args.domains, args.reference)
    print("Domain positions found:")
    for name, (s, e) in domain_positions.items():
        print(f"  {name}: {s}-{e}")

    # Step 2: extract, merge, and sort domain sequences
    extract_domains_from_alignment(args.aligned, domain_positions, fasta_out, csv_out)
    print(f"\nExtracted domains written to {fasta_out} and {csv_out}")
