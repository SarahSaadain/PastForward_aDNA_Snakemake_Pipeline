#!/usr/bin/env python
import argparse
import logging
from modules import SeqEntryReader, Writer, NormFactor
from version import __version__


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


parser = argparse.ArgumentParser(description="""
estimates the average coverage for each sequence overview entry;
notably, the average coverage corresponds to the copy number if the coverage is normalized to the coverage of single copy genes
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
    Sarah Saadain
""")
parser.add_argument('--so', type=str, default=None,dest="seqentry", required=True, help="A sequence overview (*.so) file")
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

# than normalize each entry
for se in SeqEntryReader(args.seqentry):
    cost=NormFactor.getCovStat(se,args.enddist,args.quantile)
    topr=[se.seqname,str(len(se.cov))]
    form=[]
    for c in cost:
        if c is not None:
            form.append(f"{c:.2f}",)
        else:
            form.append("na")
    topr.extend(form)
    tp="\t".join(topr)
    writer.write(tp)
