import sys, os
import numpy as np
import argparse
import re
import time
import scipy.stats.mstats as mstats
import types
import pandas as pd
import scipy.stats
from sklearn.decomposition import PCA

from .. import config
from ..tracks import Tracks


class Readstats:
    STATS = {"abs_diff", "pca",} #"speed", "hierclust", "kmeans"

    def __call__(self, *args, **kwargs):
        conf, prms = config._init_group("readstats", *args, **kwargs)

        if isinstance(prms.layers, list):
            layers = prms.layers
        else:
            layers = [prms.layers]

        for stat in layers:
            conf.tracks.layers.append(stat)
        
        io = Tracks(conf=conf)

        reads = list()

        header = True
        
        for read_id, aln in io.iter_reads():
            df = aln.to_pandas(layers, index=["aln.read_id","aln.id"])
            
            df = df.groupby(level=[0,1]).agg(prms.stats)

            sys.stdout.write(df.to_csv(sep="\t", header=header))
            header = False


def readstats(*args, **argv):
    """Perform per-read analyses of DTW alignments"""
    df = Readstats()(*args, **argv)
