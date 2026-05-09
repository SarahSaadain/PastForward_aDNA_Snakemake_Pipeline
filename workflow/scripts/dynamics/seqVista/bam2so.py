#!/usr/bin/env python
import argparse
import logging
import pysam
from modules import SeqBuilder, Writer, load_fasta
from version import __version__


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = argparse.ArgumentParser(description="""
summarize coverage for diverse features
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
    Sarah Saadain
""")
parser.add_argument('--infile', type=str, dest="infile", required=True, help="Input BAM or SAM file path")
parser.add_argument("--fasta", type=str, required=True, dest="fasta", default=None, help="the fasta file to which reads were mapped")
parser.add_argument("--mapqth", type=int, required=False, dest="mapqth", default=5, help="mapping quality threshold; below ambiguous")
parser.add_argument("--mc-snp", type=int, required=False, dest="mcsnp", default=5, help="minimum count of SNPs")
parser.add_argument("--mf-snp", type=float, required=False, dest="mfsnp", default=0.1, help="minimum frequency of SNPs")
parser.add_argument("--mc-indel", type=int, required=False, dest="mcindel", default=3, help="minimum count of indels")
parser.add_argument("--mf-indel", type=float, required=False, dest="mfindel", default=0.01, help="minimum frequency of indels")
parser.add_argument("--outfile", type=str, required=False, dest="outfile", default=None, help="output file in so format; if none is provided output will be screen")
parser.add_argument("--log-level", type=str, required=False, dest="loglevel", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)")
parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

args = parser.parse_args()
logging.getLogger().setLevel(args.loglevel)

#if no output file is provided, don't write log to screen, otherwise it will mess up the output
if args.outfile is None:
    logging.getLogger().setLevel("ERROR")

writer = Writer(args.outfile)

# load fasta from file into dict
reference_dict = load_fasta(args.fasta)

builder=None
seen_sequences = set()

infile_path = args.infile

mode = 'rb' if infile_path.lower().endswith('.bam') else 'r'
logging.info(f"Processing file: {infile_path} with mode: {mode}")

samfile = pysam.AlignmentFile(infile_path, mode)

for read in samfile:

    if read.is_unmapped:
        continue
    if read.is_secondary:
        continue
    if read.is_supplementary:
        continue

    ref_name = read.reference_name
    pos = read.reference_start  # 0-based
    mapq = read.mapping_quality if read.mapping_quality is not None else 0
    cigar = read.cigarstring
    read_sequence = read.query_sequence.upper() if read.query_sequence is not None else ''

    if cigar is None or cigar == '*':
        continue

    if builder is None:
        ref_sequence = reference_dict[ref_name]
        builder = SeqBuilder(ref_sequence, ref_name, args.mapqth)

    if ref_name != builder.seqname:
        seq_entry = builder.toSeqEntry(
            mcsnp=args.mcsnp,
            mfsnp=args.mfsnp,
            mcindel=args.mcindel,
            mfindel=args.mfindel)

        writer.write(str(seq_entry))
        seen_sequences.add(builder.seqname)

        ref_sequence = reference_dict[ref_name]
        builder = SeqBuilder(ref_sequence, ref_name, args.mapqth)

    builder.add_read(pos, cigar, mapq, read_sequence)

samfile.close()

# process the last one as well
seq_entry = None

if builder is not None:
    seq_entry = builder.toSeqEntry(args.mcsnp, args.mfsnp, args.mcindel, args.mfindel)
    writer.write(str(seq_entry))
    seen_sequences.add(seq_entry.seqname)

for ref_name, ref_sequence in reference_dict.items():
    if ref_name not in seen_sequences:
        empty_builder = SeqBuilder(ref_sequence, ref_name, args.mapqth)
        empty_entry = empty_builder.toSeqEntry(args.mcsnp, args.mfsnp, args.mcindel, args.mfindel)
        writer.write(str(empty_entry))
