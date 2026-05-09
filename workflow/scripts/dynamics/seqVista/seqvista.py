#!/usr/bin/env python3
"""
SeqVista — central command for the SeqVista analysis pipeline

Subcommands
-----------
  bam2so        Convert BAM/SAM to sequence-overview (.so) format
  normalize     Normalize .so coverage to single-copy genes
  estimate      Estimate per-entry coverage statistics from a .so file
  so2plotable   Convert .so to R-plottable format
  plot          Render .plotable files to PNG via R

Pass --help after any subcommand for its full usage.
"""

import argparse
import subprocess
import sys
from pathlib import Path
from version import __version__

_SRC = Path(__file__).parent

_SUBCOMMANDS = {
    "bam2so":      (_SRC / "bam2so.py",      "Convert BAM/SAM → .so (sequence overview)"),
    "normalize":   (_SRC / "normalize-so.py", "Normalize .so coverage to single-copy genes"),
    "estimate":    (_SRC / "estimate-so.py",  "Estimate per-entry coverage statistics"),
    "so2plotable": (_SRC / "so2plotable.py",  "Convert .so → R-plottable format"),
    "plot":        (_SRC / "run_plotable.py", "Render .plotable files to PNG via R"),
    "summarize":   (_SRC / "so2summary.py",     "Summarize .so files into a variation statistics table"),
}


def main():
    parser = argparse.ArgumentParser(
        prog="SeqVista",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument("--version", action="version", version=f"SeqVista {__version__}")

    sub = parser.add_subparsers(dest="subcommand", metavar="<subcommand>")
    for name, (_, help_text) in _SUBCOMMANDS.items():
        sub.add_parser(name, help=help_text, add_help=False)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args, remaining = parser.parse_known_args()

    if args.subcommand is None:
        parser.print_help()
        sys.exit(1)

    script, _ = _SUBCOMMANDS[args.subcommand]
    result = subprocess.run([sys.executable, str(script)] + remaining)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
