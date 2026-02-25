#!/usr/bin/env python3
"""
Run visualize-plotable.R on .plotable files from one or multiple sample folders.

Single folder mode:
    python run_plotable.py --folder Dmel01_plottable
    python run_plotable.py --folder Dmel01_plottable --output results/

Multi-folder (merged/facet) mode:
    python run_plotable.py --folders Dmel01_plottable Dmel02_plottable --output merged_results/
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


def run_rscript(input_path: Path, output_path: Path):
    """Call the R script with input and output paths."""
    cmd = ["Rscript", str(R_SCRIPT), str(input_path), str(output_path)]
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
                # Avoid double newlines between files; ensure trailing newline
                if i > 0 and content and not content.startswith(b"\n"):
                    out_fh.write(b"\n")
                out_fh.write(content)


def single_folder_mode(folder: Path, output: Path):
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
        run_rscript(src_path, out_path)

    log.info("Done.")


def multi_folder_mode(folders: list[Path], output: Path):
    """Merge same-named .plotable files across folders and plot each merged file."""
    output.mkdir(parents=True, exist_ok=True)

    # Collect all unique file names across all folders
    all_names: dict[str, list[Path]] = {}
    for folder in folders:
        for name, path in find_plotables(folder).items():
            all_names.setdefault(name, []).append(path)

    if not all_names:
        log.error("No .plotable files found in any of the provided folders.")
        sys.exit(1)

    # Warn about files not present in all folders
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
            run_rscript(merged_file, out_path)

    log.info("Done.")


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

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

    args = parser.parse_args()

    # Validate
    if args.folders is not None and args.outdir is None:
        parser.error("--outdir is required when using --folders")

    if args.folder is not None:
        if not args.folder.is_dir():
            parser.error(f"Folder not found: {args.folder}")
        output = args.outdir if args.outdir else args.folder
        single_folder_mode(args.folder, output)

    else:  # multi-folder mode
        for f in args.folders:
            if not f.is_dir():
                parser.error(f"Folder not found: {f}")
        multi_folder_mode(args.folders, args.outdir)


if __name__ == "__main__":
    main()