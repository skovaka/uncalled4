#!/usr/bin/env python

from uncalled4 import PoreModel
import argparse
import sys

parser = argparse.ArgumentParser(description="Convert binary pore model (.npz) to TSV file")
parser.add_argument("model", type=str, help="Model filename or preset")
parser.add_argument("-o", "--outfile", required=False, type=str, default=None, help="")
args = parser.parse_args(sys.argv[1:])

m = PoreModel(args.model)
out = sys.stdout if args.outfile is None else args.outfile

try:
    m.to_tsv(out)
except (BrokenPipeError,KeyboardInterrupt):
    pass
