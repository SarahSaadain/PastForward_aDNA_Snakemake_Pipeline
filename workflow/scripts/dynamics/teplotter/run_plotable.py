#!/usr/bin/env python3
"""
Run visualize-plotable.R on .plotable files from one or multiple sample folders.

Single folder mode:
    python run_plotable.py --folder Dmel01_plottable
    python run_plotable.py --folder Dmel01_plottable --output results/
    python run_plotable.py --folder Dmel01_plottable --log
    python run_plotable.py --folder Dmel01_plottable --log 1000

Multi-folder (merged/facet) mode:
    python run_plotable.py --folders Dmel01_plottable Dmel02_plottable --output merged_results/
    python run_plotable.py --folders Dmel01_plottable Dmel02_plottable --output merged_results/ --log

Authors
-------
    Robert Kofler
    Sarah Saadain
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path
import tempfile

log = logging.getLogger(__name__)


def find_plotables(folder: Path) -> dict[str, Path]:
    """Return a dict of {filename: path} for all .plotable files in folder."""
    return {f.name: f for f in folder.glob("*.plotable")}


R_SCRIPT = Path(__file__).parent / "visualize-plotable.R"


def run_rscript(input_path: Path, output_path: Path, log_arg: str | None = None):
    """Call the R script with input and output paths, and optional log flag."""
    cmd = ["Rscript", str(R_SCRIPT), str(input_path), str(output_path)]
    if log_arg is not None:
        cmd.append(log_arg)
    log.info("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=False)
    if result.returncode != 0:
        log.warning("Rscript exited with code %d for %s", result.returncode, input_path.name)


def merge_plotables(file_paths: list[Path], dest: Path):
    """Concatenate multiple .plotable files into dest."""
    with open(dest, "wb") as out_fh:
        for i, fp in enumerate(file_paths):
            with open(fp, "rb") as in_fh:
                content = in_fh.read()
                if i > 0 and content and not content.startswith(b"\n"):
                    out_fh.write(b"\n")
                out_fh.write(content)


def single_folder_mode(folder: Path, output: Path, log_arg: str | None = None):
    """Plot each .plotable file in folder independently."""
    output.mkdir(parents=True, exist_ok=True)
    plotables = find_plotables(folder)

    if not plotables:
        log.error("No .plotable files found in %s", folder)
        sys.exit(1)

    log.info("Found %d .plotable file(s) in %s", len(plotables), folder)
    for name, src_path in sorted(plotables.items()):
        out_path = output / (src_path.stem + ".png")
        log.info("[%s]", name)
        run_rscript(src_path, out_path, log_arg)

    log.info("Done.")


def multi_folder_mode(folders: list[Path], output: Path, log_arg: str | None = None):
    """Merge same-named .plotable files across folders and plot each merged file."""
    output.mkdir(parents=True, exist_ok=True)

    all_names: dict[str, list[Path]] = {}
    for folder in folders:
        for name, path in find_plotables(folder).items():
            all_names.setdefault(name, []).append(path)

    if not all_names:
        log.error("No .plotable files found in any of the provided folders.")
        sys.exit(1)

    n_folders = len(folders)
    for name, paths in sorted(all_names.items()):
        if len(paths) < n_folders:
            log.warning("'%s' missing from %d folder(s) — will merge available copies only.",
                        name, n_folders - len(paths))

    log.info("Found %d unique .plotable file name(s) across %d folder(s).", len(all_names), n_folders)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)

        for name, paths in sorted(all_names.items()):
            merged_file = tmp_path / name
            log.info("[%s] — merging %d file(s)", name, len(paths))
            merge_plotables(paths, merged_file)

            out_path = output / (Path(name).stem + ".png")
            run_rscript(merged_file, out_path, log_arg)

    log.info("Done.")


def build_log_arg(log_value: str | None) -> str | None:
    """Convert the --log argparse value to the R script flag string."""
    if log_value is None:
        return None
    if log_value == "":          # plain --log with no value
        return "--log"
    return f"--log={log_value}"  # --log=N threshold form


def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    parser = argparse.ArgumentParser(
        description="Run visualize-plotable.R on .plotable files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--folder",
        type=Path,
        metavar="DIR",
        help="Single sample folder. Each .plotable file is plotted independently.",
    )
    mode.add_argument(
        "--folders",
        nargs="+",
        type=Path,
        metavar="DIR",
        help="Multiple sample folders. Same-named .plotable files are merged before plotting.",
    )

    parser.add_argument(
        "--outdir", "-o",
        type=Path,
        metavar="DIR",
        default=None,
        help=(
            "Output directory for plots. "
            "Required when --folders is used. "
            "Defaults to the source folder when --folder is used."
        ),
    )

    parser.add_argument(
        "--log",
        nargs="?",        # 0 or 1 argument: bare --log or --log N
        const="",         # value when --log is given with no argument
        default=None,     # value when --log is not given at all
        metavar="N",
        help=(
            "Use logarithmic y-axis. "
            "Without a value (--log): always use log scale. "
            "With a value (--log 1000): auto-switch to log if max coverage exceeds N."
        ),
    )

    args = parser.parse_args()

    log_arg = build_log_arg(args.log)

    # Validate
    if args.folders is not None and args.outdir is None:
        parser.error("--outdir is required when using --folders")

    if args.folder is not None:
        if not args.folder.is_dir():
            parser.error(f"Folder not found: {args.folder}")
        output = args.outdir if args.outdir else args.folder
        single_folder_mode(args.folder, output, log_arg)

    else:
        for f in args.folders:
            if not f.is_dir():
                parser.error(f"Folder not found: {f}")
        multi_folder_mode(args.folders, args.outdir, log_arg)


if __name__ == "__main__":
    main()