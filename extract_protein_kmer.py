#!/usr/bin/env python3
"""Parallel, single-pass TargetMosaic k-mer Jaccard scan.

Each FAA is processed once. The parent process writes the same FASTA IDs and
TSV count rows as the original script, for every requested threshold.
"""
import argparse
import gzip
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from Bio import SeqIO

_QUERY_KMERS = None
_K = None
_MIN_THRESHOLD = None


def init_worker(query_kmers, k, min_threshold):
    global _QUERY_KMERS, _K, _MIN_THRESHOLD
    _QUERY_KMERS = query_kmers
    _K = k
    _MIN_THRESHOLD = min_threshold


def parse_faa(path):
    with open(path, 'rb') as f:
        zipped = f.read(2) == b'\x1f\x8b'
    opener = gzip.open if zipped else open
    with opener(path, 'rt') as handle:
        yield from SeqIO.parse(handle, 'fasta')


def scan_one(path):
    sam_id = os.path.basename(path).split('.')[0]
    query = _QUERY_KMERS
    nquery = len(query)
    hits = []
    for record in parse_faa(path):
        seq = str(record.seq)
        kmers = {seq[i:i + _K] for i in range(max(0, len(seq) - _K + 1))}
        intersection = len(query.intersection(kmers))
        union = nquery + len(kmers) - intersection
        score = intersection / union if union else 0.0
        if score >= _MIN_THRESHOLD:
            hits.append((record.id, seq, score))
    return sam_id, hits


def fasta_write(handle, identifier, seq):
    handle.write('>' + identifier + '\n')
    for i in range(0, len(seq), 60):
        handle.write(seq[i:i + 60] + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', required=True)
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument('--dir', help='Directory with .faa and .faa.gz files')
    grp.add_argument('--faa', help='One .faa or .faa.gz file')
    ap.add_argument('--k', type=int, default=5)
    ap.add_argument('--thresholds', type=float, nargs='+', default=[.10, .25, .50, .75, .90])
    ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--output-dir', default='.')
    args = ap.parse_args()
    if args.k < 1 or args.workers < 1 or not args.thresholds or any(not 0 <= t <= 1 for t in args.thresholds):
        ap.error('k and workers must be positive; thresholds must lie between 0 and 1')
    thresholds = sorted(set(args.thresholds))
    labels = [str(round(t * 100)) for t in thresholds]
    if len(set(labels)) != len(labels):
        ap.error('Two thresholds map to the same output filename percentage')
    records = list(SeqIO.parse(args.input, 'fasta'))
    if len(records) != 1:
        ap.error('--input must contain exactly one FASTA record')
    seq = str(records[0].seq)
    query_kmers = {seq[i:i + args.k] for i in range(max(0, len(seq) - args.k + 1))}
    if args.faa:
        files = [args.faa]
    else:
        files = [str(Path(args.dir) / f) for f in os.listdir(args.dir) if f.endswith(('.faa', '.faa.gz'))]
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    from contextlib import ExitStack
    with ExitStack() as stack:
        output = {}
        for t, label in zip(thresholds, labels):
            fa = stack.enter_context(open(Path(args.output_dir) / f'matches_{label}.fasta', 'w'))
            count = stack.enter_context(open(Path(args.output_dir) / f'counts_{label}.tsv', 'w'))
            output[t] = (fa, count)
        with ProcessPoolExecutor(max_workers=args.workers,
                                 initializer=init_worker,
                                 initargs=(query_kmers, args.k, thresholds[0])) as pool:
            for number, (sam_id, hits) in enumerate(pool.map(scan_one, files, chunksize=1), start=1):
                for threshold, (fa, count) in output.items():
                    matches = 0
                    for record_id, protein_seq, similarity in hits:
                        if similarity >= threshold:
                            matches += 1
                            fasta_write(fa, f'{sam_id}_{record_id}|similarity={similarity:.2f}', protein_seq)
                    count.write(f'{sam_id}\t{matches}\n')
                if number % 1000 == 0:
                    print(f'Processed {number}/{len(files)} FAA files', flush=True)
    print(f'Finished {len(files)} FAA files using {args.workers} workers; thresholds: {thresholds}')


if __name__ == '__main__':
    main()
