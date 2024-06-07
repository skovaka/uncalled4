import sys, os
import numpy as np
import argparse
import re
import time
import types
import pandas as pd
import scipy.stats
from collections import defaultdict

from .. import config
from ..tracks import RefstatsSplit, ALL_REFSTATS, parse_layers
from ..argparse import Opt, comma_split
from ..ref_index import str_to_coord

def refstats_single(tracks):
    pass

def refstats(conf):
    """Calculate per-reference-coordinate statistics"""
    from ..tracks import Tracks

    t0 = time.time()

    conf.tracks.shared_refs_only = True
    conf.shared_refs_only = True
    conf.read_index.load_signal = False

    tracks = Tracks(conf=conf)
    conf = tracks.conf

    if getattr(conf, "outfile", None) is not None:
        outfile = open(conf.outfile,"w")
    else:
        outfile = sys.stdout

    #if conf.tracks.io.processes == 1

    stats = RefstatsSplit(conf.refstats, len(tracks.alns))
    layers = list(parse_layers(conf.tracks.layers))

    #if conf.verbose_refs:
    columns = ["ref_name", "ref", "strand"] + [".".join([track.name, "cov"]) for track in tracks.alns]
    #else:
    #    columns = ["ref"]

    for track in tracks.alns:
        name = track.name
        for group, layer in layers:
            for stat in stats.layer:
                columns.append(".".join([track.name, group, layer, stat]))

    for group,layer in layers:
        for stat in stats.compare:
            columns.append(".".join([stat, group, layer, "stat"]))
            columns.append(".".join([stat, group, layer, "pval"]))

    #columns.append("kmer")

    outfile.write("\t".join(columns)+"\n")

    strand_chrs = np.array(["-","+"])

    for chunk in tracks.iter_refs():
        chunk.prms.refstats = conf.refstats
        chunk.prms.refstats_layers = layers
        stats = chunk.refstats#
        if stats is None: 
            stats = chunk.calc_refstats()

        if len(stats) == 0: continue

        stats = pd.concat({chunk.coords.name : stats}, axis=0) #\
        strands = strand_chrs[stats.index.get_level_values("seq.fwd").astype(int)]
        stats = stats.reset_index(level="seq.fwd",drop=True) \
                     .set_index(strands,append=True)

        #          .reset_index(level="seq.fwd")
        #stats["seq.strand"] = stats["seq.fwd"].replace({True:"+",False:"-"})
        #stats = stats.set_index("seq.strand",append=True)
        #del stats["seq.fwd"]
        outfile.write(stats.to_csv(sep="\t",header=False,na_rep=0))

    outfile.close()
