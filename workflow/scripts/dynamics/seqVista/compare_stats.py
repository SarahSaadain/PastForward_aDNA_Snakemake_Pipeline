#!/usr/bin/env python3
"""
Compare per-sequence stats across samples and flag sequences of interest.

Input is two or more .stats.tsv files produced by so2stats.py.
Output has one row per sequence with per-sample columns side by side,
cross-sample metrics, and flags.

Copy number flags (median_cov used as SCG-normalised copy number proxy)
-----------------------------------------------------------------------
CN_FC   log2(max/min) of median_cov across samples exceeds --cn-fc threshold.
        Catches relative shifts at low copy number (e.g. 1 -> 5, log2FC = 2.32).
CN_ABS  max - min of median_cov across samples exceeds --cn-abs threshold.
        Catches large absolute shifts at high copy number (e.g. 50 -> 70).

A sequence is flagged if ANY flag applies. Multiple flags are pipe-separated.

Usage
-----
    python compare_stats.py --stats s1.stats.tsv s2.stats.tsv --cn-fc 1.0 --cn-abs 10

    python compare_stats.py --stats *.stats.tsv --cn-fc 1.0 --cn-abs 10 --outfile results.tsv

    # only show flagged sequences in terminal and output
    python compare_stats.py --stats *.stats.tsv --cn-fc 1.0 --cn-abs 10 --flagged-only

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

log = logging.getLogger(__name__)

# Per-sample columns to carry into the wide output
# These are written as {metric}__{sampleid} in the final table
PER_SAMPLE_COLS = [
    "median_cov", "mad_cov", "cv_cov", "max_cov", "frac_low",
    "n_snps", "snp_density", "median_alt",
]


# ── loading ───────────────────────────────────────────────────────────────────

def load_stats(files: list[Path]) -> pd.DataFrame:
    """
    Load and concatenate stats TSVs produced by so2stats.py.
    Each file must contain a 'sampleid' column; duplicate (seqid, sampleid)
    combinations across files raise an error.
    """
    frames = []
    for fp in files:
        if not fp.is_file():
            log.error("File not found: %s", fp)
            sys.exit(1)
        log.info("Loading %s", fp.name)
        df = pd.read_csv(fp, sep="\t")
        frames.append(df)

    combined = pd.concat(frames, ignore_index=True)

    # catch duplicate (seqid, sampleid) pairs early
    dupes = combined.duplicated(subset=["seqid", "sampleid"])
    if dupes.any():
        log.error(
            "Duplicate (seqid, sampleid) pairs found:\n%s",
            combined[dupes][["seqid", "sampleid"]].to_string(index=False)
        )
        sys.exit(1)

    return combined


# ── wide format ───────────────────────────────────────────────────────────────

def pivot_to_wide(long_df: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot long-format stats (one row per seqid x sampleid) to wide format
    (one row per seqid).

    Column naming: {metric}__{sampleid}

    Cross-sample columns added for copy number (median_cov):
      cn_min      : lowest median_cov across samples
      cn_max      : highest median_cov across samples
      cn_abs      : cn_max - cn_min  (absolute copy number range)
      cn_log2fc   : log2(cn_max / cn_min)  (relative shift; NaN if cn_min == 0)
    """
    if long_df.empty:
        return pd.DataFrame()

    samples = sorted(long_df["sampleid"].unique())

    # pivot all per-sample metrics at once
    wide = long_df.pivot(index="seqid", columns="sampleid", values=PER_SAMPLE_COLS)
    wide.columns = [f"{metric}__{sample}" for metric, sample in wide.columns]
    wide = wide.reset_index()

    # seq_len is identical across samples for the same sequence
    seq_len_df = long_df.groupby("seqid")["seq_len"].max().reset_index()
    wide = wide.merge(seq_len_df, on="seqid", how="left")

    # --- cross-sample copy number metrics ---
    cn_cols = [f"median_cov__{s}" for s in samples if f"median_cov__{s}" in wide.columns]

    wide["cn_min"]    = wide[cn_cols].min(axis=1)
    wide["cn_max"]    = wide[cn_cols].max(axis=1)
    wide["cn_abs"]    = (wide["cn_max"] - wide["cn_min"]).round(3)

    # log2FC is undefined when cn_min is 0; set to NaN in that case
    wide["cn_log2fc"] = np.where(
        wide["cn_min"] > 0,
        np.log2(wide["cn_max"] / wide["cn_min"]),
        np.nan
    ).round(3)

    return wide, samples


# ── flagging ──────────────────────────────────────────────────────────────────

def add_flags(wide_df: pd.DataFrame,
              cn_fc_threshold: float,
              cn_abs_threshold: float) -> pd.DataFrame:
    """
    Add 'flag' (pipe-separated reasons) and 'flagged' (bool) columns.

    CN_FC   cn_log2fc >= cn_fc_threshold
            Catches relative copy number shifts (sensitive at low copy number).
            Example: 1 -> 5  gives log2FC 2.32; threshold of 1.0 means a 2x change.

    CN_ABS  cn_abs >= cn_abs_threshold
            Catches large absolute shifts (sensitive at high copy number).
            Example: 50 -> 70 gives cn_abs 20; set threshold accordingly.
    """
    df = wide_df.copy()

    flags = pd.DataFrame(index=df.index)
    flags["CN_FC"]  = df["cn_log2fc"].notna() & (df["cn_log2fc"] >= cn_fc_threshold)
    flags["CN_ABS"] = df["cn_abs"]            >= cn_abs_threshold

    df["flag"]    = flags.apply(lambda row: "|".join(k for k, v in row.items() if v), axis=1)
    df["flagged"] = df["flag"].str.len() > 0
    return df


# ── column ordering ───────────────────────────────────────────────────────────

def reorder_columns(df: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    """
    Readable column order:
      seqid | seq_len | flag | flagged
      cross-sample CN metrics
      per-metric blocks: all samples for median_cov, then mad_cov, ...
    """
    head    = ["seqid", "seq_len", "flag", "flagged"]
    cn_summary = ["cn_min", "cn_max", "cn_abs", "cn_log2fc"]
    per_sample = [
        f"{metric}__{s}"
        for metric in PER_SAMPLE_COLS
        for s in samples
        if f"{metric}__{s}" in df.columns
    ]
    ordered = head + cn_summary + per_sample
    rest    = [c for c in df.columns if c not in ordered]
    return df[ordered + rest]


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    parser = argparse.ArgumentParser(
        description="Compare per-sequence stats across samples and flag sequences of interest.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--stats", nargs="+", required=True, type=Path, metavar="FILE",
        help="Two or more .stats.tsv files produced by so2stats.py.",
    )
    parser.add_argument(
        "--outfile", "-o", type=Path, metavar="FILE", default=None,
        help="Output TSV file (default: comparison.tsv).",
    )
    parser.add_argument(
        "--flagged-only", action="store_true",
        help="Only output sequences that have at least one flag.",
    )

    thresh = parser.add_argument_group("flagging thresholds")
    thresh.add_argument(
        "--cn-fc", type=float, default=2, metavar="FLOAT",
        help=(
            "log2 fold-change threshold for CN_FC flag. "
            "1.0 = 2x change, 2.32 = 5x change. "
            "Example: --cn-fc 1.0"
        ),
    )
    thresh.add_argument(
        "--cn-abs", type=float, default=10, metavar="FLOAT",
        help=(
            "Absolute copy number difference threshold for CN_ABS flag. "
            "Example: --cn-abs 10  flags sequences where max-min median_cov >= 10."
        ),
    )

    args = parser.parse_args()

    if len(args.stats) < 2:
        parser.error("--stats requires at least two files to compare.")

    outfile = args.outfile if args.outfile else Path("comparison.tsv")

    # --- load and process ---
    long_df = load_stats(args.stats)
    wide, samples = pivot_to_wide(long_df)
    wide = add_flags(wide, args.cn_fc, args.cn_abs)
    wide = reorder_columns(wide, samples)

    # sort: flagged sequences first, then by cn_abs descending within each group
    wide = wide.sort_values(["flagged", "cn_abs"], ascending=[False, False]).reset_index(drop=True)

    if args.flagged_only:
        n_before = len(wide)
        wide     = wide[wide["flagged"]]
        log.info("--flagged-only: %d / %d sequences retained", len(wide), n_before)

    wide.to_csv(outfile, sep="\t", index=False)
    log.info("Comparison written to %s  (%d sequences)", outfile, len(wide))

    # terminal preview: key columns only
    preview_cols = ["seqid", "flag", "cn_min", "cn_max", "cn_abs", "cn_log2fc"]
    # add one median_cov column per sample so copy number is visible per sample
    preview_cols += [f"median_cov__{s}" for s in samples if f"median_cov__{s}" in wide.columns]
    preview_cols  = [c for c in preview_cols if c in wide.columns]

    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_colwidth", 50)
    print(f"\nTop sequences (flagged first, then by cn_abs; full details in {outfile}):")
    print(wide[preview_cols].head(20).to_string(index=False))

    # flag summary
    print(f"\nFlagged: {wide['flagged'].sum()} / {len(wide)} sequences")
    all_flags = wide["flag"].str.split("|").explode()
    all_flags = all_flags[all_flags.str.len() > 0]
    if not all_flags.empty:
        print("Reason counts:")
        print(all_flags.value_counts().to_string())


if __name__ == "__main__":
    main()