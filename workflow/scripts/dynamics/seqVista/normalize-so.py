#!/usr/bin/env python
import argparse
from modules import SeqEntryReader, NormFactor, Writer
from version import __version__
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = argparse.ArgumentParser(description="""
normalizes the coverage for seqentries
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
    Sarah Saadain
""")
parser.add_argument('--so', type=str, default=None,dest="seqentry", required=True, help="A sequence overview (so) file")
parser.add_argument("--scg-end", type=str, required=False, dest="scgend", default="_scg", help="the ending by which to recognize single copy gens (scg)")
parser.add_argument("--end-distance", type=int, required=False, dest="enddist", default=100, help="distance from ends for normalizing")
parser.add_argument("--exclude-quantile", type=int, required=False, dest="quantile", default=25, help="exclude the most extreme coverage quantiles for normalizing")
parser.add_argument("--outfile", type=str, required=False, dest="outfile", default=None, help="output file in so format; if none is provided output will be screen")
parser.add_argument("--log-level", type=str, required=False, dest="loglevel", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)")
parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

args = parser.parse_args()
logging.getLogger().setLevel(args.loglevel)

#if no output file is provided, don't write log to screen, otherwise it will mess up the output
if args.outfile is None:
    logging.getLogger().setLevel("ERROR")

writer = Writer(args.outfile)
# first get the normalization factor
normfactor = NormFactor.getNormalizationFactor(args.seqentry, args.scgend ,args.enddist, args.quantile)

if normfactor == 0:
    logging.error("Normalization factor is zero. This may indicate that no single copy genes were found or that the coverage is zero or very low. Please check your input data.")
    exit(1)

# than normalize each entry
for se in SeqEntryReader(args.seqentry):
    sen = se.normalize(normfactor)
    writer.write(str(sen))
