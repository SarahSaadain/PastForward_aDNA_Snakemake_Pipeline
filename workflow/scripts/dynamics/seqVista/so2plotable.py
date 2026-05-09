#!/usr/bin/env python
import argparse
import logging
import os
from modules import PlotableFormater, SeqEntryReader, Writer, load_bed
from version import __version__


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

padto=9
def format_col(topr:list):
    if padto==0:
        return "\t".join(topr)
    tf=[]
    for i in range(padto):
        if(i< len(topr)):
            tf.append(str(topr[i]))
        else:
            tf.append("")
    return "\t".join(tf)


parser = argparse.ArgumentParser(description="""
converts sequence overview (so) files to plottable format
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
    Sarah Saadain
""")
parser.add_argument('--so', type=str, default=None, dest="so", required=True, help="A sequence overview (*.so) file")
parser.add_argument("--seq-ids", type=str, required=False, dest="seqids", default="ALL", help="IDs of the entries that should be plotted; separated by comma; can also be 'ALL'")
parser.add_argument("--sample-id", type=str, required=False, dest="sampleid", default="x", help="the ID of current sample")
parser.add_argument("--prefix", type=str, required=False, dest="prefix", default="", help="the prefix for the output file; only valid with --outdir")
parser.add_argument("--outdir", type=str, required=False, dest="outputdir", default=None, help="the output directory; a plotable will be written for each fasta entry")
parser.add_argument("--outfile", type=str, required=False, dest="outfile", default=None, help="output file in plotable format")
parser.add_argument("--mask-bed", type=str, required=False, dest="maskbed", default=None, help="a BED file for masking; regions in the file will be masked (0-based coordinates)")
parser.add_argument("--mask-ymax", type=int, required=False, dest="ymax", default=None, help="mask positions with coverage exceeding this value")
parser.add_argument("--log-level", type=str, required=False, dest="loglevel", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)")
parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

args = parser.parse_args()
logging.getLogger().setLevel(args.loglevel)

if args.outfile is not None and args.outputdir is not None:
    parser.error("invalid parameters; either provide --outfile or --outdir; not both")

if args.prefix != "" and args.outputdir is None:
    parser.error("invalid parameters; --prefix only works if --outdir is provided")

writer = Writer(args.outfile)
tomask = load_bed(args.maskbed)

#if no output file is provided, don't write log to screen, otherwise it will mess up the output
if args.outfile is None:
    logging.getLogger().setLevel("ERROR")

if "," in args.seqids:
    seqset = set(args.seqids.split(","))
else:
    seqset = set([args.seqids])

is_print_all_requested = (args.seqids.upper() == "ALL")

if args.outputdir is not None:
    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

prefix = args.prefix

for se in SeqEntryReader(args.so):
    if is_print_all_requested or se.seqname in seqset:
        tmp = PlotableFormater.prepareForPrint(se, args.sampleid, tomask, args.ymax)
        tr = [format_col(i) for i in tmp]
        tp = "\n".join(tr)

        if args.outputdir is not None:
            filename = se.seqname
            filename = filename.replace("/", "_")
            filename = filename.replace(" ", "_")
            filename = prefix + filename + ".plotable"
            full_path = os.path.join(args.outputdir, filename)
            with open(full_path, "w") as f:
                f.write(tp)
        else:
            writer.write(tp)
