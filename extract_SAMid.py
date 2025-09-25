#!/usr/bin/env python3
"""
This script extracts SAM IDs from the Bakta Export file.
Only assemblies with correct Species identification by both bakta and gtdbtk.classification will be extracted.
"""

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Extract SAM IDs for a given genus and species.")
    parser.add_argument("genus", help="Genus name (e.g., Pseudomonas)")
    parser.add_argument("species", help="Species name (e.g., aeruginosa)")
    parser.add_argument(
        "-i", "--input", default="bakrep-export.tsv",
        help="Input Bakta export TSV file (default: bakrep-export.tsv)"
    )
    parser.add_argument(
        "-o", "--output", default="SAMID.tsv",
        help="Output file to write SAM IDs (default: SAMID.tsv)"
    )
    
    args = parser.parse_args()

    # Read file
    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    # Filter rows
    filtered_values = df.loc[
        (df["gtdbtk.classification.genus"] == args.genus)
        & (df["gtdbtk.classification.species"] == f"{args.genus} {args.species}"),
        "#id",
    ]

    # Write results
    with open(args.output, "w") as f:
        for value in filtered_values:
            f.write(f"{value}\n")

if __name__ == "__main__":
    main()
