#!/usr/bin/env python

import sys, os
import numpy as np
import pandas as pd
import argparse
from gtfparse import read_gtf

parser = argparse.ArgumentParser(description="Translates transcriptome to genome coordinates")
parser.add_argument("gtf", type=str, help="GTF annotation file")
parser.add_argument("infile", type=str, help="Tab-delimited file (or set by --delim) containing transcript names and coordinates")
parser.add_argument("-t", "--trans-col", type=str, default="transcript_id", help="Name of column with transcript IDs (matching GTF transcript_id)")
parser.add_argument("-p", "--pos-col", type=str, default="transcript_position", help="Name of column with transcript coordinates")
#parser.add_argument("-s", "--strand-col", type=str, default="strand")
parser.add_argument("-d", "--delim", type=str, default="\t", help="Field seperator")
#parser.add_argument("-g", "--gene-to-chr", action="store_true")
parser.add_argument("-o", "--outfile", required=False, type=str, default=sys.stdout, help="Output file")
args = parser.parse_args(sys.argv[1:])

anno = read_gtf(args.gtf)

if not isinstance(anno, pd.DataFrame):
    anno = anno.to_pandas()

input_iter  = pd.read_csv(
    args.infile, 
    sep=args.delim, 
    index_col=[args.trans_col,args.pos_col],
    chunksize=10000,
)

anno = anno[anno["feature"] == "exon"]
rev = anno.index[anno["strand"] == "-"]
tmp = anno.loc[rev,"start"]
anno.loc[rev, "start"] = -anno.loc[rev, "end"]
anno.loc[rev, "end"] = -tmp
anno.sort_values(["gene_id","transcript_id","start"], inplace=True)#.set_index(["gene_id","transcript_id"])
anno["len"] = anno["end"] - anno["start"] + 1
anno["tloc"] = anno.groupby("transcript_id")["len"].cumsum() - anno["len"]
anno = anno.reset_index().set_index("transcript_id").sort_values("tloc").sort_index(kind="stable")
anno.sort_values(["seqname","start"], inplace=True)

#trans_df[args.trans_col] = trans_df[args.trans_col].str.slice(0,15)
anno_trans = anno.index.unique().sort_values()

header = True
for chunk in input_iter:
    for tid,df in chunk.groupby(level=args.trans_col):
        if not tid in anno_trans: continue
        trans = anno.loc[[tid]]

        tlocs = df.index.get_level_values(args.pos_col).to_numpy()
        exons = trans.iloc[ trans["tloc"].searchsorted(tlocs, side="right") - 1 ].copy()

        offs = (exons["start"] + tlocs - exons["tloc"]).abs() - 1

        coords = pd.DataFrame({
            "gene_id" : exons["gene_id"],
            "chr" : exons["seqname"],
            "chr_loc" : offs,
            args.pos_col : tlocs
        }).set_index(args.pos_col, append=True)

        df_out = pd.concat([df, coords], axis=1)
        df_out.to_csv(args.outfile, sep=args.delim, mode="append", header=header)
        header = False
