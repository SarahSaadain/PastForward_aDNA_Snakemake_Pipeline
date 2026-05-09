#!/usr/bin/env python3
"""
Compute per-sequence coverage and SNP statistics from a single .so file.

Output is a TSV with one row per sequence, intended as input for compare_stats.py.

Usage
-----
    python so2stats.py --so sample1.so --sample-id Dmel1933_SL07 --outfile Dmel1933_SL07.stats.tsv

Authors
-------
    (adapted from run_plotable.py structure by Robert Kofler / Sarah Saadain)
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from modules import SeqEntryReader

log = logging.getLogger(__name__)


def compute_stats(filepath: Path, sample_id: str) -> pd.DataFrame:
    """
    Stream through the .so file and compute per-sequence statistics in one pass.
    Never materialises per-position rows — memory is O(number of sequences).

    Coverage metrics
    ----------------
    seq_len     : number of positions in the sequence
    median_cov  : median coverage across all positions; used as copy number proxy
                  when data is SCG-normalised (1 = single copy, 2 = duplicated, etc.)
    mad_cov     : median absolute deviation of coverage; robust spread measure
    cv_cov      : MAD / median; scale-independent coverage variation
    max_cov     : peak coverage (useful for spotting sharp spikes)
    frac_low    : fraction of positions with coverage < 0.1 (proxy for absent/deleted)

    SNP metrics
    -----------
    n_snps      : total number of alt allele observations across all positions
    snp_density : n_snps per 100 bp (length-normalised)
    median_alt  : median alt allele count across SNP positions
    """
    records = []

    for se in SeqEntryReader(str(filepath)):
        vals = np.asarray(se.cov, dtype=np.float32)
        seq_len = len(vals)

        median_cov = float(np.median(vals))
        mad_cov    = float(np.median(np.abs(vals - median_cov)))
        cv_cov     = mad_cov / median_cov if median_cov > 0 else np.nan
        max_cov    = float(np.max(vals))
        frac_low   = float(np.mean(vals < 0.1))

        alt_counts = [
            count
            for snp in se.snplist
            for base, count in (("A", snp.ac), ("T", snp.tc), ("C", snp.cc), ("G", snp.gc))
            if base != snp.refc and count > 0
        ]
        n_snps      = len(alt_counts)
        snp_density = (n_snps / seq_len * 100) if seq_len > 0 else 0.0
        median_alt  = float(np.median(alt_counts)) if alt_counts else 0.0

        records.append({
            "seqid":       se.seqname,
            "sampleid":    sample_id,
            "seq_len":     seq_len,
            "median_cov":  round(median_cov,  3),
            "mad_cov":     round(mad_cov,     3),
            "cv_cov":      round(cv_cov,      3),
            "max_cov":     round(max_cov,     3),
            "frac_low":    round(frac_low,    3),
            "n_snps":      n_snps,
            "snp_density": round(snp_density, 3),
            "median_alt":  round(median_alt,  3),
        })

    return pd.DataFrame(records) if records else pd.DataFrame()


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    parser = argparse.ArgumentParser(
        description="Compute per-sequence stats from a single .so file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--so", required=True, type=Path, metavar="FILE",
        help="Input .so file.",
    )
    parser.add_argument(
        "--sample-id", required=True, dest="sample_id", metavar="ID",
        help="Sample identifier written into the output (e.g. Dmel1933_SL07).",
    )
    parser.add_argument(
        "--outfile", "-o", type=Path, metavar="FILE", default=None,
        help="Output TSV file (default: <sample_id>.stats.tsv).",
    )

    args = parser.parse_args()

    if not args.so.is_file():
        parser.error(f"File not found: {args.so}")

    outfile = args.outfile if args.outfile else Path(f"{args.sample_id}.stats.tsv")

    log.info("Parsing %s (sample_id=%s)", args.so.name, args.sample_id)
    stats = compute_stats(args.so, args.sample_id)

    if stats.empty:
        log.error("No coverage data found in %s. Exiting.", args.so)
        sys.exit(1)

    stats.to_csv(outfile, sep="\t", index=False)
    log.info("Stats written to %s  (%d sequences)", outfile, len(stats))

    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_colwidth", 50)
    print(stats.to_string(index=False))


if __name__ == "__main__":
    main()
