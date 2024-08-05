"""Compute, compare, and query alignment layers

subcommand options:
compare  Compare different alignment methods on the same set of reads"""

import sys, os
import numpy as np
import argparse
import re
import time
import scipy.stats.mstats as mstats
import types
import pandas as pd
import scipy.stats

from .. import config
from ..tracks import Tracks, parse_layers


def compare(conf):
    """Compute distance between alignments of the same reads"""
    t = time.time()

    group_b = "moves" if conf.moves else "dtw"

    if conf.moves:
        conf.tracks.layers += [("mvcmp", "dist")] 
    else:
        conf.tracks.layers += [("dtwcmp", "dist")] 
    conf.read_index.load_signal = False

    tracks = Tracks(conf=conf)

    t = time.time()
    t_all = time.time()

    t = time.time()

    for read_id,alns in tracks.iter_reads():
        if not isinstance(alns, list) or len(alns) != 2:
            raise RuntimeError("compare must be run on two BAM files")
        a,b = alns
        if a is None or b is None: continue
        a.calc_dtwcmp(b.instance)
        tracks.write_alignment(a)

def dump(conf):
    """Output DTW alignment paths and statistics"""

    tracks = Tracks(conf=conf)
    tracks.set_layers(conf.layers)

    layer_groups = {group for group,_ in parse_layers(conf.layers)}

    need_cmp = "cmp" in layer_groups
    need_mvcmp = "mvcmp" in layer_groups

    header = True


    for read_id,tracks in tracks.iter_reads():
        for track in tracks:
            if header:
                columns = track.layers.index.names + [
                    ".".join([c for c in col if len(c) > 0]) 
                    for col in track.layers.columns]
                header = False
            layers = track.layers.dropna(axis=0, how="any")
            sys.stdout.write(layers.to_csv(sep="\t", header=False))
        t = time.time()
