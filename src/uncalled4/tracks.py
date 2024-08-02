from . import config
from .params import BASECALLER_PROFILES
from .pore_model import PoreModel, Sequence
from .ref_index import load_index, RefCoord, str_to_coord
from .read_index import ReadIndex, Fast5Reader, Slow5Reader

import re
from . import config, ref_index 

#from .layers import LAYER_META
#from .tracks import LAYER_META

import sys, os
import numpy as np
import pandas as pd
import collections
from collections import defaultdict, namedtuple
import scipy.stats

from time import time

import _uncalled4

from collections import namedtuple
_Layer = namedtuple("_Layer", ["dtype", "label", "default"], defaults=[None,None,False])

ALN_LAYERS = {
    "start" : _Layer("Int32", "Sample Start", True),
    "length" : _Layer("Int32", "Sample Length", True),
    "current" : _Layer(np.float32, "Current (norm)", True),
    "current_sd" : _Layer(np.float32, "Stdv (norm)", True),
    "events" : _Layer(np.float32, "Event Count", True),
    "start_sec" : _Layer(np.float32, "Time (s)", False),
    "length_sec" : _Layer(np.float32, "Dwell (s)", False),
    "dwell" : _Layer(np.float32, "Dwell (ms)", False),
    "end" : _Layer("Int32", "Sample End", False),
    "middle" : _Layer(np.float32, "Sample Middle", False),
    "middle_sec" : _Layer(np.float32, "Sample Middle", False),
    "model_diff" : _Layer(np.float32, "Model Norm. Diff.", False),
    "abs_diff" : _Layer(np.float32, "Abs. Model Diff.", False),
}

LAYERS = {
    "seq" : {
        "mpos" : _Layer("Int64", "Mirror Ref.", True),
        "pos" : _Layer("Int64", "Reference Coord.", False),
        "pac" : _Layer("Int64", "Packed Ref. Coord.", False),
        "kmer" : _Layer("Int32", "Kmer", True),
        "current" : _Layer(str, "Model Mean (norm)", False),
        "name" : _Layer(str, "Reference Name", False),
        "fwd" : _Layer(bool, "Forward", False),
        "strand" : _Layer(str, "Strand", False),
        "base" : _Layer(str, "Base", False),
    }, "aln" : {
        "track" : _Layer(str, "Track ID"),
        "id" : _Layer("Int64", "Aln. ID"),
        "read_id" : _Layer(str, "Read ID"),
    }, "dtw" : ALN_LAYERS, "moves" : ALN_LAYERS,
    "dtwcmp" : {
        "aln_a" : _Layer("Int32", "Compare alignment A", False),
        "aln_b" : _Layer("Int32", "Compare alignment B", False),
        "group_b" : _Layer(str, "Compare type", False),
        "jaccard" : _Layer(np.float32, "Jaccard Distance", True),
        "dist" : _Layer(np.float32, "Sig-to-Ref Dist", True),
        "current" : _Layer(np.float32, "Current Cmp (norm)", True),
        "current_sd" : _Layer(np.float32, "Stdv Cmp (norm)", True),
    }, "mvcmp" : {
        "jaccard" : _Layer(np.float32, "Jaccard Distance", True),
        "dist" : _Layer(np.float32, "Sig-to-Ref Dist", True),
        "current" : _Layer(np.float32, "Current Cmp (norm)", True),
        "current_sd" : _Layer(np.float32, "Stdv Cmp (norm)", True),
    }
}

LAYER_META = pd.concat([
    pd.concat({
        group : pd.DataFrame(layers, index=_Layer._fields).transpose()
    }, names=("group","layer"))  
    for group, layers in LAYERS.items()
])

DEFAULT_LAYERS = LAYER_META.index[LAYER_META["default"]]

def get_layer_offs(name):
    spl = re.split("[-+]", name)
    if len(spl) > 2:
        raise ValueError(f"Invalid layer name: {name}")
    layer = spl[0]
    if len(spl) == 2:
        offs = int(spl[1])
        if '-' in name:
            offs = -offs
    else:
        offs = 0
    return layer,offs

def parse_layer(layer):
    if isinstance(layer, str):
        spl = layer.split(".")
    elif isinstance(layer, tuple):
        spl = layer
    else:
        raise ValueError("Layer must be string or tuple")

    if len(spl) == 2:
        group,layer = spl
    elif len(spl) == 1:
        if layer not in LAYER_META.index.get_level_values(0):
            raise ValueError(f"Unknown layer or group '{layer}'")
        group = layer
        _,layers = DEFAULT_LAYERS.get_loc_level(layer)
        for layer in layers:
            yield (group, layer)
        return
    else:
        raise ValueError(f"Invalid layer: \"{layer}\"")

    #find upstream or downstream layers
    spl = re.split("[-+]", layer)
    name,_ = get_layer_offs(layer)
    if not (group, name) in LAYER_META.index:
        raise ValueError(f"Invalid layer \"{group}.{layer}\"")

    #must include centered layer to compute up/downstream versions
    if len(spl) == 2:
        yield (group,spl[0])

    yield (group, layer)

def parse_layers(layers):
    if layers is None:
        return pd.Index([])

    if isinstance(layers, str):
        layers = layers.split(",")

    ret = list()

    for layerstr in layers:
        for layer in parse_layer(layerstr):
            ret.append(layer)

    return pd.Index(ret).unique()


PANDAS_DTYPES = {
    "int16" : pd.Int16Dtype(),
    "int32" : pd.Int32Dtype(),
    "int64" : pd.Int64Dtype(),
    "uint16" : pd.UInt16Dtype(),
    "uint32" : pd.UInt32Dtype(),
    "uint64" : pd.UInt64Dtype(),
}

class AlnDF:
    #sample_rate = 4000

    def __init__(self, seq, start=None, length=None, current=None, current_sd=None):
        self.seq = seq
        self._extra = dict()

        if isinstance(start, _uncalled4._AlnDF):
            self.instance = start
        else:
            if current is None:
                current = np.zeros(0, np.float32)
            if current_sd is None:
                current_sd = np.zeros(0, np.float32)

            self.instance = _uncalled4._AlnDF(self.seq.index, start, length, current, current_sd)

        self.instance.mask(self.na_mask)

    def set_layer(self, name, vals):
        if not len(vals) == len(self): 
            raise ValueError(f"'{name}' values not same length as alignment")
        self._extra[name] = np.array(vals)
        
    @property
    def na_mask(self):
        return self.samples.mask.to_numpy()

    @property
    def start(self):
        return (self.samples.starts.to_numpy())

    @property
    def end(self):
        return (self.samples.ends.to_numpy())

    @property
    def length(self):
        return (self.samples.lengths.to_numpy())

    @property
    def start_sec(self):
        return (self.start / self.seq.model.PRMS.sample_rate)

    @property
    def length_sec(self):
        return (self.length / self.seq.model.PRMS.sample_rate)

    @property
    def middle(self):
        return (self.start + self.length / 2).astype(np.float32)

    @property
    def middle_sec(self):
        return self.middle / self.seq.model.PRMS.sample_rate

    @property
    def model_diff(self):
        return np.array(self.current) - self.seq.current

    @property
    def model_diff(self):
        return np.array(self.current) - self.seq.current

    @property
    def dwell(self):
        return self.length #1000 * self.length / self.seq.model.PRMS.sample_rate
    
    def _get_series(self, layer):
        vals = getattr(self, layer, None)
        if vals is None or len(vals) == 0:
            return None
        ret = pd.Series(vals, copy=True)
        dtype = PANDAS_DTYPES.get(ret.dtype.name, None)
        na = np.nan
        if dtype is not None:
            ret = ret.astype(dtype)
            na = pd.NA
        if len(self.na_mask) > 0:
            ret[~self.na_mask] = na
        return ret

    def to_pandas(self, layers=None):
        if layers is None:
            layers = ["start", "length", "current", "current_sd"]

        cols = dict()
        for name in layers:
            vals = self._get_series(name)
            if vals is not None:
                cols[name] = vals
        cols["index"] = self.index.expand()
        df = pd.DataFrame(cols)
        return df.set_index("index")

    def __len__(self):
        return len(self.instance)

    def __getattr__(self, name):
        ret = getattr(self.instance, name, None)

        if ret is None:
            layer,offs = get_layer_offs(name)
            vals = getattr(self, layer, None)
            if vals is None:
                raise AttributeError(f"AlnDF has no attribute '{name}'")
            if offs < 0:
                ret = np.pad(vals,(-offs,0),constant_values=0)[:offs]
            else:
                ret = np.pad(vals,(0,offs),constant_values=0)[offs:]
        return ret

class CmpDF:
    def __init__(self, seq, instance):
        self.seq = seq
        self.instance = instance

    def to_pandas(self, layers=None):
        if layers is None:
            df = pd.DataFrame({
                "dist" : self.instance.dist, 
                "jaccard" : self.instance.jaccard,
                "current" : self.instance.current,
                "current_sd" : self.instance.current_sd,
            })
        else:
            df = pd.DataFrame({
                l : getattr(self, l) for l in layers if len(getattr(self,l,[])) > 0
            })
        
        df["index"] = self.index.expand()

        return df.set_index("index")

    def __len__(self):
        return len(self.instance)

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise AttributeError(f"AlnDF has no attribute '{name}'")
        return self.instance.__getattribute__(name)

class Alignment:
    def __init__(self, aln_id, read, seq, sam, track_name):
        self.seq = seq
        self.sam = sam
        self.track = track_name

        if isinstance(read, str):
            read_id = read
            self.read = None
        else:
            read_id = read.id
            self.read = read

        if isinstance(self.seq.model, _uncalled4.PoreModelU16):
            Super = _uncalled4._AlignmentU16
        elif isinstance(self.seq.model, _uncalled4.PoreModelU32):
            Super = _uncalled4._AlignmentU32
        else:
            raise ValueError(f"Unknown PoreModel type: {model.instance}")

        self.instance = Super(aln_id, read_id, seq.instance)

        self.dtw = AlnDF(seq, self.instance._dtw)
        self.moves = AlnDF(seq, self.instance._moves)
        self.mvcmp = CmpDF(seq, self.instance._mvcmp)
        self.dtwcmp = CmpDF(seq, self.instance._dtwcmp)

        self.norm_scale = self.norm_shift = None

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise AttributeError(f"Alignment has no attribute '{name}'")
        return self.instance.__getattribute__(name)
    
    def set_norm(self, scale, shift):
        self.norm_scale = scale
        self.norm_shift = shift
        
    def set_dtw(self, df):
        if isinstance(df, AlnDF):
            df = df.instance
        self.instance.set_dtw(df)
        
    def set_moves(self, df):
        if isinstance(df, AlnDF):
            df = df.instance
        self.instance.set_moves(df)

    def to_pandas(self, layers=None, index=None, join_index=True, bounds=None):
        vals = dict()

        layers = parse_layers(layers)
        if index is not None:
            index = parse_layers(index)
            #layers = layers.union(index)
            layers = layers.append(index)

        idx = self.seq.mpos#.expand().to_numpy()
        if bounds is not None:
            pos = self.seq.pos
            idx = idx[(pos >= bounds.start) & (pos < bounds.end)]

        for name in layers.unique(0):
            _,group_layers = layers.get_loc_level(name)

            if name == "aln":
                vals[name] = pd.DataFrame({
                    l : getattr(self, l) for l in group_layers
                }, index=idx)
            else:
                group = getattr(self, name, [])
                if len(group) > 0:
                    vals[name] = group.to_pandas(group_layers).reindex(idx)#.reset_index(drop=True)

        df = pd.concat(vals, axis=1, names=["group", "layer"]).reset_index(drop=True)
        
        if index is not None:
            df.set_index(list(index), inplace=True)
            if join_index:
                df.index.names = [".".join(idx) for idx in index]

        return df#.set_index(index)

    Attrs = namedtuple("Attrs", [
        "track", "id", "read_id", "ref_name", "ref_start", "ref_end", 
        "fwd", "sample_start", "sample_end", "coord", "norm_scale", "norm_shift"
    ])

    #@property
    def attrs(self):
        samp_start = 1000000000
        samp_end = 0
        for df in [self.dtw, self.moves]:
            if not df.empty():
                samp_start = min(samp_start, df.samples.start)
                samp_end = max(samp_end, df.samples.end)

        return self.Attrs(
            self.track, self.id, self.read_id, self.seq.coord.name, self.seq.coord.get_start(), self.seq.coord.get_end(),
            self.seq.fwd, samp_start, samp_end, self.seq.coord, self.norm_scale, self.norm_shift
        )

class AlnTrack:
    def __init__(self, *args, **kwargs):
        if isinstance(args[0], AlnTrack):
            self._init_slice(*args, **kwargs)
        else:
            self._init_new(*args, **kwargs)

        self.mat = None

    def _init_new(self, track_id, name, desc, conf, model=None, fast5s=None):
        self.id = track_id
        self.name = name
        self.desc = desc
        self.conf = conf
        self.coords = None

        self.fast5s = fast5s #TODO get rid of this

        if model is not None:
            self.model = model 
        elif len(conf.pore_model.name) > 0:
            try:
                self.model = PoreModel(params=conf.pore_model)
            except FileNotFoundError:
                sys.stderr.write(f"Warning: pore model not found '{conf.pore_model.name}'\n")
                conf.pore_model.name = ""
                self.model = PoreModel(params=conf.pore_model)
        else:
            self.model = None

    def _init_slice(self, p, coords=None, alignments=None, layers=None, order=["fwd", "ref_start"]):
        self._init_new(p.id, p.name, p.desc, p.conf, p.model, p.fast5s)

    def _parse_layers(self, df):
        if df.index.names[0] == "seq.pac":
            df = df.rename(index=self.coords.pac_to_pos, level=0)
            df.index.names = ("seq.pos", "aln.id")

            refs = df.index.get_level_values(0)
            if len(df) > 1 and refs[0] > refs[1]:
                df = df.iloc[::-1]

        if self.conf.tracks.mask_indels is not None and ("moves","indel") in df.columns:
            df = df[df["moves","indel"].abs() < self.conf.tracks.mask_indels]


        if self.conf.tracks.mask_skips is not None and ("dtw","events") in df.columns:
            skips = df["dtw","events"] < 1
            if self.conf.tracks.mask_skips == "all":
                df = df[~skips]

            elif self.conf.tracks.mask_skips == "keep_best":
                sdf = df[skips].set_index(("dtw","start"), append=True)
                sdf["diff"] = (sdf["dtw","current"] - self.model[sdf["dtw","kmer"]]).abs()
                grp = sdf.groupby(level=2, group_keys=False)
                keep_idx = grp["diff"].idxmin()
                sdf = sdf.drop(keep_idx).reset_index(level=2)
                df = df.drop(sdf.index)
            else:
                raise ValueError("Unknown masking mode: " + self.conf.tracks.mask_skips)
            
        return df

    def set_data(self, coords, alignments, layers, order=["fwd", "ref_start"]):

        self.coords = coords
        self.alignments = alignments
        self.layers = self._parse_layers(layers)

        if not self.coords.stranded and (self.all_fwd or self.all_rev):
            self.coords = self.coords.ref_slice(fwd=self.all_fwd)

        isnone = [coords is None, alignments is None, layers is None]
        if np.all(isnone) or len(alignments) == 0:
            return
        elif np.any(isnone):
            raise ValueError("Must specify AlnTrack coords, alignments, and layers")
        self.alignments = self.alignments.sort_values(order)


        self.has_fwd = np.any(self.alignments['fwd'])
        self.has_rev = not np.all(self.alignments['fwd'])

    @property
    def all_fwd(self):
        return np.all(self.alignments["fwd"])

    @property
    def all_rev(self):
        return not np.any(self.alignments["fwd"])

    def slice(self, coords=slice(None), aln_ids=None, reads=None, order=["fwd","ref_start"]):
        if self.empty:
            return AlnTrack(self, coords, self.alignments, self.layers)
        layer_refs = self.layers.index.get_level_values("seq.pos")


        layers = self.layers.loc[(layer_refs >= coords.refs.min()) & (layer_refs <= coords.refs.max())]

        if reads is not None:
            if aln_ids is not None:
                raise ValueError("Only one of 'aln_ids' and 'reads' can be specified")
            aln_ids = self.alignments.index[self.alignments["read_id"].isin(reads)]

        if aln_ids is not None: 
            layer_alns = layers.index.get_level_values("aln.id")
            aln_ids = layer_alns.unique().intersection(aln_ids)
            layers = layers.loc[layer_alns.isin(aln_ids)] 
        else:
            aln_ids = slice(None)

        alignments = self.alignments.loc[aln_ids]

        return AlnTrack(self, coords, alignments, layers, order=order)
        

    @property
    def empty(self):
        return self.coords is None or len(self.alignments) == 0

    def _group_layers(self, group, layers):
        return pd.concat({group : layers}, names=["group", "layer"], axis=1)

    def _aln_id_or_default(self, aln_id):
        if aln_id is None:
            if len(self.alignments) == 1:
                return self.alignments.index[0]
            raise ValueError("Must specify aln_id for Track with more than one alignment loaded")
        return aln_id
        
    def add_layer_group(self, group, layers, aln_id, overwrite):
        aln_id = self._aln_id_or_default(aln_id)

        layers["aln.id"] = aln_id
        df = layers.set_index("aln.id", append=True).reorder_levels(["seq.pos","aln.id"])

        meta = LAYER_META.loc[group]
        df = df[meta.index.intersection(df.columns)]
        df = df.astype(meta.loc[df.columns, "dtype"], copy=False)

        df = self._parse_layers(pd.concat({group : df}, names=["group", "layer"], axis=1))

        if self.layers is None or overwrite:
            self.layers = df
        else:
            self.layers = pd.concat([self.layers, df], axis=1)

        return df

    def aln_ref_coord(self, aln_id):
        return RefCoord(*self.alignments[["ref_name","ref_start","ref_end","fwd"]].loc[aln_id])

    def has_group(self, group):
        return group in self.layers.columns.get_level_values(0)


    def calc_layers(self, layers):
        for group, layer in layers:
            if group in self.layers and (group, layer) not in self.layers.columns:
                meta = LAYER_META.loc[(group,layer)]

                #Make sure layer dependencies exist
                if not self.empty and (meta["deps"] is None or len(self.layers.columns.intersection(meta["deps"])) == len(meta["deps"])):

                    fn = meta["fn"]
                    if fn is None:
                        continue
                    vals = fn(self)
                    self.layers[group,layer] = vals

    def load_mat(self):
        df = self.layers.copy()
        df = df.reset_index()

        self.mat = df.pivot(index="aln.id", columns=["seq.pos"]) \
                     .rename_axis(("group","layer","seq.pos"), axis=1) \
                     .reindex(self.coords.refs, axis=1, level=2) \
                     .sort_index(axis=1)

        self.mat = self.mat.reindex(self.alignments.index, copy=False)

        self.width = len(self.coords.refs)
        self.height = len(self.alignments)

    def sort_coord(self, layer, pos):
        order = (-self.mat[layer,pos].fillna(0)).argsort()
        self.sort(order)

    def sort(self, order):
        self.mat = self.mat.iloc[order]
        self.alignments = self.alignments.iloc[order]

    def cmp(self, other, calc_jaccard, calc_dist):
        groups_b = other.alignments.groupby("read_id")

        df = pd.DataFrame(
            columns=["aln_b", "group_b", "dist", "jaccard"],
            index = self.layers.index
        )

        for id_a, aln_a in self.alignments.iterrows():
            read_id = aln_a["read_id"]
            self.calc_layers([("dtw","end")])
            dtw_a = self.layers.loc[(slice(None),id_a),"dtw"][["start","end"]]

            for id_b, aln_b in groups_b.get_group(read_id).iterrows():
                self._compare_alns(dtw_a, other, id_b, "dtw", df, calc_jaccard, calc_dist)

        df["group_b"] = "dtw"
        return df#.set_index(["aln_b", "group_b"], append=True)

    def mvcmp(self, other, calc_jaccard, calc_dist):
        if other != self:
            groups_b = other.alignments.groupby("read_id")

        df = pd.DataFrame(
            columns=["aln_a", "aln_b", "group_b", "dist", "jaccard"],
            index = self.layer_refs
        )

        for id_a, aln_a in self.alignments.iterrows():
            self.calc_layers([("dtw","end")])
            dtw = self.layers.loc[(slice(None),id_a),"dtw"][["start","end"]]

            if dtw is None:
                continue

            if other == self:
                self._compare_alns(dtw, self, id_a, "moves", df, calc_jaccard, calc_dist)
            else:
                read_id = aln_a["read_id"]
                for id_b, aln_b in groups_b.get_group(read_id).iterrows():
                    self._compare_alns(dtw, other, id_b, "moves", df)

        df["group_b"] = "moves"
        return df#.set_index(["aln_b", "group_b"], append=True)

    def _compare_alns(self, aln_a, other, id_b, group, df, calc_jaccard=True, calc_dist=True):

        other.calc_layers([(group,"end")])

        aln_b = other.layers \
                   .loc[(slice(None),id_b),group][["start","end"]] \
                   .reset_index(level="aln.id") \
                   .rename(columns={"aln.id" : "aln_b"})
                

        alns_a = aln_a.index.unique(1)
        if len(alns_a) != 1:
            raise ValueError("Can only compare two alignments at a time")
        
        has_pac = "seq.pac" in aln_b.index.name
        flip = self.all_fwd == self.conf.is_rna

        def coords(df, track):
            df = df.dropna().sort_index().reset_index()
            if flip:
                end = -df["start"]
                df["start"] = -df["end"]
                df["end"] = end
            return AlnCoords(df)


        coords_a = coords(aln_a, self)
        coords_b = coords(aln_b, other)

        compare = Compare(coords_a, coords_b)

        cmp_df = pd.DataFrame(compare.to_numpy()).dropna(how="all")

        if has_pac:
            cmp_df["seq.pac"] = self.coords.pos_to_pac(pd.Index(cmp_df["seq.pos"]))
            cmp_df = cmp_df.set_index("seq.pac")
        else:
            cmp_df = cmp_df.set_index("seq.pos")
        df["jaccard"] = cmp_df["jaccard"]
        df["dist"] = cmp_df["dist"]
        df["aln_a"] = alns_a[0]
        df["aln_b"] = id_b



    @property
    def read_ids(self):
        return self.alignments["read_id"].unique()

    @property
    def layer_aln_ids(self):
        return self.layers.index.get_level_values("aln.id")

    @property
    def layer_pacs(self):
        return self.coords.pos_to_pac(self.layer_refs)

    @property
    def layer_fwds(self):
        return self.alignments.loc[self.layer_aln_ids, "fwd"]

    @property
    def layer_strands(self):
        return self.layer_fwds.map({True : "+", False : "-"})

    @property
    def layers_pac_index(self):
        index = pd.MultiIndex.from_arrays([self.coords.pos_to_pac(self.layer_refs), self.layers.index.get_level_values("aln.id")], names=["seq.pac","aln.id"])
        return self.layers.set_index(index, drop=True)

    @property
    def layers_desc_index(self):
        index = pd.MultiIndex.from_arrays(
            [self.layer_refs, self.layer_strands, self.layer_aln_ids])
        return pd.concat({self.coords.ref_name : self.layers.set_index(index, drop=True)}, names=[("seq.pos","name"),("seq.pos","coord"),("seq.pos","strand"),"aln.id"])

        stats.index = pd.MultiIndex.from_product([
            [chunk.coords.ref_name], stats.index, ["+" if chunk.coords.fwd else "-"]
        ])

    @property
    def layer_refs(self):
        return self.layers.index.get_level_values("seq.pos")

_REFSTAT_AGGS = {
    "cov" : "count",
    "mean" : "mean", 
    "median" : "median", 
    "q5" : (lambda x: np.quantile(x, 0.05)),
    "q95" : (lambda x: np.quantile(x, 0.95)),
    "q25" : (lambda x: np.quantile(x, 0.25)),
    "q75" : (lambda x: np.quantile(x, 0.75)),
    "stdv" : "std", 
    "var"  : "var",
    "skew" : scipy.stats.skew,
    "kurt" : scipy.stats.kurtosis,
    "min"  : "min", 
    "max"  : "max",
}

_REFSTAT_CMPS = {
    "ks" : lambda df: scipy.stats.stats.ks_2samp(df.loc[df.index.unique(0)],df.loc[df.index.unique(1)],mode="asymp")
}

LAYER_REFSTATS = {"cov", "min", "max", "mean", "median", "stdv", "var", "skew", "kurt", "q25", "q75", "q5", "q95"}
COMPARE_REFSTATS = {"ks"}
ALL_REFSTATS = LAYER_REFSTATS | COMPARE_REFSTATS

REFSTAT_LABELS = {
    "cov" : "Coverage", 
    "min" : "Minimum", 
    "max" : "Maximum", 
    "mean" : "Mean", 
    "median" : "Median", 
    "stdv" : "Std. Dev.", 
    "var" : "Variance", 
    "skew" : "Skewness", 
    "kurt" : "Kurtosis",
    "ks" : "KS",
    "q5" : "5% Quantile",
    "q25" : "25% Quantile",
    "q75" : "75% Quantile",
    "q95" : "95% Quantile",
}

BUILTIN_TRACKS = {"_refstats", "_layerstats", "_readstats"}

CMP_GROUPS = {"dtwcmp"}

class RefstatsSplit:
    def __init__(self, stats, track_count):
        self.layer = [s for s in stats if s in LAYER_REFSTATS]
        self.compare = [s for s in stats if s in COMPARE_REFSTATS]

        self.layer_agg = [_REFSTAT_AGGS[s] for s in self.layer]

        if len(self.layer) + len(self.compare) != len(stats):
            bad_stats = [s for s in stats if s not in ALL_REFSTATS]
            raise ValueError("Unknown stats: %s (options: %s)" % (", ".join(bad_stats), ", ".join(ALL_REFSTATS)))

        if len(self.compare) > 0 and track_count != 2:
            raise ValueError("\"%s\" stats can only be computed using exactly two tracks" % "\", \"".join(self.compare))

from .io import M6anet, TSV, Eventalign, Tombo, BAM, ModelTrainer, INPUTS, OUTPUTS, INPUT_PARAMS, OUTPUT_PARAMS

class Tracks:
    def __init__(self, *args, **kwargs):
        self.read_index = kwargs.get("read_index", None)
        if self.read_index is not None:
            del kwargs["read_index"]

        if len(args) > 0 and isinstance(args[0], Tracks):
            self._init_slice(*args, **kwargs)
        else:
            self._init_new(*args, **kwargs)

    def _init_new(self, *args, model=None, **kwargs):


        self.conf, self.prms = config._init_group("tracks", copy_conf=True, *args, **kwargs)
        self.prms.refstats_layers = list(parse_layers(self.prms.refstats_layers))

        self.alignments = None
        self.layers = None

        self.new_alignment = False
        self.new_layers = set()

        self.alns = list()
        self.models = dict()
        self._tracks = dict()
        self.new_alignment = False
        self.new_layers = set()
        
        
        if self.read_index is None:
            self.read_index = ReadIndex(conf=self.conf) #(self.conf.read_index.paths, self.prms.read_filter, self.conf.read_index.read_index, self.conf.read_index.recursive)

        
        pm = self.conf.pore_model
        if not pm.has_workflow():
            flowcell, kit = self.read_index.default_model
            if flowcell is not None and kit is not None:
                pm.flowcell = flowcell
                pm.kit = kit

        if model is not None:
            self.set_model(model)
            pm.name = model.PRMS.name
        else:
            if len(pm.name) == 0 and pm.has_workflow():
                pm.name = PoreModel.PRESET_MAP.loc[pm.get_workflow(), "preset_model"]
                sys.stderr.write(f"Auto-detected flowcell='{pm.flowcell}' kit='{pm.kit}'\n")
            self.model = None


        self._init_io()

        if self.prms.basecaller_profile is None:# and pm.has_workflow():
            bp = self.model.name[:self.model.name.rfind("_")]
            if not bp in BASECALLER_PROFILES:
                wf = pm.get_workflow()
                if wf in PoreModel.PRESET_MAP.index:
                    bp = PoreModel.PRESET_MAP.loc[pm.get_workflow(), "basecaller_profile"]
                else:
                    bp = None

            self.prms.basecaller_profile = bp


        self.set_layers(self.prms.layers)


        self._aln_track_ids = [t.id for t in self.alns]

        if self.prms.ref is not None:
            self.index = load_index(self.model, self.prms.ref)
        else:
            self.index = None

        #    raise RuntimeError("Failed to load reference index")
        self.refstats = None

        #self.coords = self._ref_bounds_to_coords(self.prms.ref_bounds)
        self.coords = self.prms.ref_bounds

        if self.coords is not None:
            if not self.coords.has_bounds:
                self.coords.start = 0
                self.coords.end = self.index.get_ref_len(self.coords.name)

            if len(self._aln_track_ids) > 0 and self.prms.load_mat:
                self.load()

    def _init_slice(self, parent, coords, reads=None):
        self.conf = parent.conf 
        self.prms = parent.prms

        self.set_layers(self.prms.layers)
        self.prms.refstats_layers = list(parse_layers(self.prms.refstats_layers))

        self.index = parent.index
        self.read_index = parent.read_index
        self._aln_track_ids = parent._aln_track_ids
        self.refstats = None
        self.new_alignment = parent.new_alignment
        self.new_layers = parent.new_layers

        self.inputs = parent.inputs
        self.output = parent.output
        self.output_track = parent.output_track

        self.coords = coords
        self._tracks = dict()
        self.alns = list()
        self.model = parent.model

        self._add_tracks(parent._tracks)

        if parent.alignments is not None:
            if coords != parent.coords or reads is not None:
                mask = np.minimum(parent.alignments["ref_start"], coords.start) < np.maximum(parent.alignments["ref_end"], coords.end)
                if reads is not None:
                    mask &= parent.alignments["read_id"].isin(reads)
                self.alignments = parent.alignments[mask] #pd.concat(track_alns, axis=0, names=["track", "id"])
                ids = parent.layers.index
                self.layers = parent.layers.loc[ids]
            else:
                self.alignments = parent.alignments
                self.layers = parent.layers
        else:
            self.alignments = None
            self.layers = None
    

    def _add_tracks(self, tracks):
        for name,track in tracks.items():
            self._add_track(name, track)

    def _add_track(self, name, track):
        self._tracks[name] = track
        
        if name in BUILTIN_TRACKS:
            setattr(self, name[1:], track)
        elif isinstance(track, AlnTrack):
            self.alns.append(track)
        else:
            raise ValueError("Unrecognized track type: " + str(track))

    @property
    def track_count(self):
        return len(self.alns)

    @property
    def all_empty(self):
        for t in self.alns:
            if not t.empty:
                return False
        return True

    @property
    def any_empty(self):
        for t in self.alns:
            if t.empty:
                return True
        return False

    @property
    def track_names(self):
        return list(self._tracks.keys())

    def _ref_bounds_to_coords(self, ref_bounds):
        if ref_bounds is not None:
            if isinstance(ref_bounds, str):
                ref_bounds = str_to_coord(ref_bounds)
            elif isinstance(ref_bounds, tuple):
                ref_bounds = RefCoord(*ref_bounds)
            return self.index.get_coord_space(ref_bounds, self.conf.is_rna)
        return None

    def set_layers(self, layers):
        self.prms.layers = layers

        self.db_layers = list()
        self.fn_layers = list()
        self.cmp_layers = list()
        for group, layer in parse_layers(layers):
            if group in CMP_GROUPS:
                self.cmp_layers.append((group, layer))
            else:
                self.db_layers.append((group, layer))

    def __len__(self):
        return len(self._tracks)

    def keys(self):
        return self._tracks.keys()

    def __getitem__(self, i):
        return self._tracks[i]

    def _init_io(self):

        self.inputs = list()
        self.bam_in = None

        track_count = 0

        #if names is not None:

        for name,Cls in INPUTS.items():
            fnames = getattr(self.prms.io, name, None)
            if fnames is None: continue
            if isinstance(fnames, str):
                fnames = [fnames]
            assert isinstance(fnames, list)
                
            for filename in fnames:
                io = Cls(filename, False, self, track_count)

                p = io.conf.read_index
                if not self.prms.io.buffered and self.conf.read_index.load_signal:
                    if p.paths is not None:
                        self.read_index.load_paths(p.paths)
                    self.read_index.load_index_file(p.read_index)

                self.conf.load_config(io.conf)
                self.inputs.append(io)

                if isinstance(io, BAM) and self.bam_in is None:
                    self.bam_in = io

                track_count += len(io.aln_tracks)
                for track in io.aln_tracks:
                    self._add_track(track.name, track)

        in_prms = [getattr(self.prms.io, p) is not None for p in INPUT_PARAMS]
        out_prms = [getattr(self.prms.io, p) is not None for p in OUTPUT_PARAMS]

        if np.sum(out_prms) > 1:
            raise ValueError("No more than one output can be specified")

        if self.model is None:
            for track in self.alns:
                if track.model is not None:
                    self.model = track.model
                    break
        if self.model is None:
            raise ValueError("Failed to detect pore model. Please set '--pore-model', or '--flowcell' and '--kit'")

        if np.any(out_prms):
            out_format = OUTPUT_PARAMS[out_prms][0]
            filename = getattr(self.prms.io, out_format)
            new_track = True
            if out_format == "m6anet_out":
                self.output = M6anet(filename, True, self, track_count)
            elif out_format == "tsv_out":
                self.output = TSV(filename, True, self, track_count)
            elif out_format == "eventalign_out":
                self.output = Eventalign(filename, True, self, track_count)
            elif out_format == "bam_out":
                self.output = BAM(filename, True, self, track_count)
            elif out_format == "model_dir":
                self.output = ModelTrainer(filename, True, self, track_count)

            track_count += len(self.output.aln_tracks)
            for track in self.output.aln_tracks:
                self._add_track(track.name, track)

            self.output_track = self.output.aln_tracks[0].name
        else:
            self.output = None
            self.output_track = self.inputs[0].aln_tracks[0].name


    def set_model(self, model):
        self.conf.load_group("pore_model", model.PRMS)
        self.model = model
        for track in self.alns:
            track.model = model

    def get_read_fast5(self, read_id):
        if self.read_index is not None and self.read_index.read_files is not None:
            return self.read_index.read_files.loc[read_id]
        for io in self.inputs:
            if getattr(io, "read_id_in", None) == read_id:
                return io.fast5_in
        raise RuntimeError("Could not determine fast5 filename, may need to provide fast5 index (-x)")

    def aln_layers(self, layer_filter=None):
        ret = self.layers.columns
        if layer_filter is None:
            return ret
        return ret.intersection(layer_filter)
    
    def _track_or_default(self, track_name):
        if track_name is None:
            return self._tracks[self.output_track]
        elif track_name in self._tracks:
            return self._tracks[track_name]
        raise ValueError(f"Unknown track: {track_name}")

    def collapse_events(self, dtw, read=None):
        dtw["cuml_mean"] = dtw["length"] * dtw["current"]

        grp = dtw.groupby(level=0)

        lengths = grp["length"].sum()


        dtw = pd.DataFrame({
            "start"  : grp["start"].min().astype("uint32"),
            "length" : lengths.astype("uint32"),
            "current"   : grp["cuml_mean"].sum() / lengths,
            "kmer" : grp["kmer"].first(),
            "events" : grp["start"].count(),
        })

        if read is not None:
            dtw["current_sd"] = [
                np.std(
                    read.get_norm_signal(dtw.loc[i, "start"], dtw.loc[i, "start"]+dtw.loc[i, "length"])
                ) for i in dtw.index
            ]

        skip_counts = dtw["start"].value_counts().loc[dtw["start"]].to_numpy()
        dtw["events"] /= skip_counts

        return dtw.sort_index()

    def write_dtw_events(self, events=None, track_name=None, aln_id=None, read=None):
        if np.any(events.index.duplicated()):
            events = self.collapse_events(events, read=read)

        self.add_layers("dtw", events, track_name, aln_id, False, read)
        
    def add_layers(self, group, layers, track_name=None, aln_id=None, overwrite=False, read=None):
        track = self._track_or_default(track_name)

        if layers.index.names[0] == "mpos":
            layers = layers.set_index(self.index.mpos_to_pos(layers.index))
        elif layers.index.names[0] == "seq.pac":
            layers = layers.set_index(self.index.pac_to_pos(layers.index))

        layers.index.names = ("seq.pos",)

        track.add_layer_group(group, layers, aln_id, overwrite)

        self.new_layers.add(group)


    def init_alignment(self, track_name, aln_id, read, bam, *coord_args):
        track = self._track_or_default(track_name)

        fwd = not bam.is_reverse

        if self.index is not None:
            coords = RefCoord(bam.reference_name, *coord_args, fwd)
            seq = self.index.query(coords)

        else:
            kmers = self.model.str_to_kmers(bam.query_sequence)
            if not fwd:
                kmers = self.model.kmer_comp(kmers)
            if fwd == self.model.PRMS.reverse:
                kmers = self.model.kmer_rev(kmers)[::-1]
            seq = self.model[kmers]
            seq.set_coord(RefCoord(read.id, self.model.shift, len(kmers)+self.model.shift))

        return Alignment(aln_id, read, seq, bam, track_name)

    def write_alignment(self, aln):
        if self.output is not None:
            out = self.output
        else:
            out = self.inputs[0]

        out.write_alignment(aln)

    def set_read(self, read):
        self.output.read = read

    @property
    def input_count(self):
        return len(self.alns)

    def _verify_read(self):
        if len(self._aln_track_ids) == 0:
            raise RuntimeError("No input tracks have been loaded")

    def slice(self, ref_start=None, ref_end=None, reads=None, order=["fwd","ref_start"], tracks=None, full_overlap=False, shared_reads=False):
        if self.coords is None:
            raise IndexError("Cannot slice empty Tracks")

        hasbounds = (ref_start is not None, ref_end is not None)
        if np.all(hasbounds):
            coords = self.coords.intersection(RefCoord(self.coords.name, ref_start, ref_end))
        elif np.any(hasbounds):
            raise IndexError(f"Invalid bounds {ref_start}-{ref_end}")
        else:
            coords = self.coords

        stranded = True
        fwd = None

        if full_overlap:
            reads = self.get_full_overlap(reads)

        if shared_reads:
            if reads is None:
                reads = self.get_shared_reads()
            else:
                reads = self.get_shared_reads().intersection(reads)
            order = "read_id"

        return Tracks(self, coords, reads)

    def get_shared_reads(self):
        all_ids = self.alignments["read_id"]
        read_ids = pd.Index(all_ids.loc[self.alns[0].id])
        for track in self.alns[1:]:
            read_ids = read_ids.intersection(all_ids.loc[track.name])
        return read_ids

    def get_all_reads(self):
        return self.alignments["read_id"].drop_duplicates().to_numpy()

    def get_full_overlap(self, read_ids=None):
        rmin = self.coords.refs.min()
        rmax = self.coords.refs.max()
        for track in self.alns:
            alns = track.alignments[(track.alignments["ref_start"] <= rmin) & (track.alignments["ref_end"] >= rmax)]
            if read_ids is None:
                read_ids = pd.Index(alns["read_id"])
            else:
                read_ids = read_ids.intersection(alns["read_id"])
        return read_ids
            
    def load(self, ref_bounds=None, full_overlap=None, read_filter=None, load_mat=None):
        self._verify_read()
        self.new_alignment = True
        self.new_layers = set()

        if ref_bounds is not None:
            self.coords = self._ref_bounds_to_coords(ref_bounds)

        if read_filter is None:
            read_filter = self.prms.read_filter

        if self.coords is None:
            raise ValueError("Must set ref bounds")

        if full_overlap is None:
            full_overlap = self.prms.full_overlap

        track_alns = dict()
        track_layers = dict()

        for io in self.inputs:
            alns, layers = io.query(self.conf.tracks.layers, self.coords, ["aln.id","seq.pos","seq.fwd"], full_overlap=full_overlap, read_id=read_filter)


            track_alns.update(alns)
            track_layers.update(layers)

        self.layers = pd.concat(track_layers, axis=0, names=["aln.track", "aln.id", "seq.pos","seq.fwd"])
        self.alignments = pd.concat(track_alns, axis=0, names=["track", "id"]).sort_index()
    
        self.calc_refstats()

        if load_mat is None:
            load_mat = self.prms.load_mat
        if load_mat:
            self.init_mat()

        return self.alns

    def init_mat(self):
        df = self.layers.reset_index()

        self.mat = df.pivot(index=["aln.track","aln.id"], columns=["seq.pos"]) 
        self.mat = self.mat.rename_axis(("group","layer","seq.pos"), axis=1) 
        self.mat = self.mat.reindex(pd.RangeIndex(self.coords.start, self.coords.end), axis=1, level=2) 
        self.mat = self.mat.sort_index(axis=1)

        order = self.alignments.sort_values(["fwd", "ref_start"]).index

        self.mat = self.mat.reindex(order, copy=False)

        self.width = len(self.coords)
        self.height = len(self.alignments)

    def calc_refstats(self, cov=True):
        if self.prms.refstats is None or len(self.prms.refstats) == 0 or len(self.prms.refstats_layers) == 0 or self.alignments is None:
            self.refstats = None
            return None

        stats = RefstatsSplit(self.prms.refstats, len(self.alns))

        groups = self.layers[self.prms.refstats_layers].groupby(level=["aln.track", "seq.pos", "seq.fwd"])

        refstats = list()

        if cov:
            g,l = self.prms.refstats_layers[0]
            covs = groups.size()
            covs = covs.rename("cov") \
                    .reset_index() \
                    .pivot(index=["seq.pos","seq.fwd"], columns="aln.track") \
                    .reorder_levels([1,0],axis=1)
                    #}, axis=1
                #)
            covs.columns = pd.MultiIndex.from_product(covs.columns.levels +  [[""],[""]])
            refstats.append(covs)

        if len(stats.layer_agg) > 0:
            layerstats = groups.agg(stats.layer_agg) \
                             .reset_index() \
                             .pivot(index=["seq.pos","seq.fwd"], columns="aln.track") \
                             .reorder_levels([3,0,1,2],axis=1)
            refstats.append(layerstats)

        if len(stats.compare) > 0:
            if len(self.layers.index.unique(0)) != 2:
                def ks(df):
                    idx = df.index.levels[0]
                    d = dict()
                    for col in df.columns:
                        d[col+("ks",)] = pd.Series([pd.NA, pd,NA], index=["stat","pval"], name=df.name)
                    return pd.concat(d, axis=0)
            
            else:
                def ks(df):
                    tids = df.index.levels[0]
                    d = dict()
                    for col in df.columns:
                        a = df.loc[tids[0],col].dropna()
                        b = df.loc[tids[1],col].dropna()
                        if len(a) > 1 and len(b) > 1:
                            vals = scipy.stats.stats.ks_2samp(a,b,mode="asymp")
                            d[col+("ks",)] = pd.Series(vals, index=["stat","pval"], name=df.name)
                        else:
                            d[col+("ks",)] = pd.Series([np.nan, np.nan], index=["stat","pval"], name=df.name)
                    df = pd.concat(d, axis=0)
                    return df

            cmp_groups = self.layers[self.prms.refstats_layers].groupby(level=["seq.pos","seq.fwd"])
            cmpstats = cmp_groups.apply(ks).reorder_levels([2,0,1,3],axis=1)
            refstats.append(cmpstats)

        refstats = pd.concat(refstats, axis=1).sort_index()


        self.refstats = refstats.dropna()

        self._tracks["_refstats"] = self.refstats

        return refstats

    def iter_refs(self, ref_bounds=None):
        layer_iters = [
            io.iter_refs(
                self.db_layers, 
                chunksize=self.prms.io.ref_chunksize)
            for io in self.inputs]
        
        chunks = [(pd.DataFrame([]), pd.DataFrame([])) for i in layer_iters]
        chunk_hasnext = np.ones(len(chunks), bool)

        while np.any(chunk_hasnext):
            pac_start = np.inf
            pac_end = np.inf

            all_pacs = pd.Index([])

            strand = -1
            all_fwd = True
            all_rev = True

            all_done = True
            for i in range(len(chunks)):
                if not chunk_hasnext[i]: continue 
                if len(chunks[i][1]) == 0:
                    chunks[i] = next(layer_iters[i], (None, None))
                    if chunks[i][0] is None:
                        chunk_hasnext[i] = False
                        continue
                alns,layers = chunks[i]
                pac_end = min(pac_end, layers.index.get_level_values("seq.pac").max())
                all_done = False

            if all_done: break

            pos_end = self.index.pac_to_pos(pac_end)

            ret_layers = dict()
            ret_alns = dict()

            for i,(alns,layers) in enumerate(chunks):
                if layers is None: 
                    continue
                ret_layers[i] = layers.loc[:pac_end]
                ret_alns[i] = alns.loc[ret_layers[i].index.unique("aln.id")]
                leftover_layers = layers.loc[pac_end+1:] #drop(index=ret_layers[i].index)
                leftover_alns = alns.loc[leftover_layers.index.unique("aln.id")]
                chunks[i] = (leftover_alns, leftover_layers)

            alns = pd.concat(ret_alns, names=["track","id"])
            layers = pd.concat(ret_layers, names=["aln.track","seq.pos","seq.fwd","aln.id"])
            pos = self.index.pac_to_pos(layers.index.levels[1])
            layers.index = layers.index.set_levels(pos, level=1)

            coords = RefCoord(alns.iloc[0]["ref_name"], pos.min(), pos.max()+1)

            ret = self._tables_to_tracks(coords, alns, layers)

            yield ret

    def iter_reads(self, read_filter=None, ref_bounds=None, full_overlap=False, max_reads=None, ignore_bam=False, return_tracks=False):

        if read_filter is None and self.read_index is not None:
            read_filter = self.read_index.read_filter
        
        if ref_bounds is not None and not isinstance(ref_bounds, RefCoord):
            ref_bounds = RefCoord(ref_bounds)

        gen = self.iter_reads_db(read_filter, ref_bounds, full_overlap, max_reads, ignore_bam, return_tracks)

        return gen
            
    def iter_reads_slice(self, reads=None, ref_bounds=None):
        all_reads = self.get_all_reads()
        if reads is not None:
            all_reads = all_reads.intersection(reads)

        if ref_bounds is None or not ref_bounds.has_bounds():
            ref_start = ref_end = None
        else:
            ref_start = ref_bounds.start
            ref_end = ref_bounds.end
        
        for read_id in all_reads:
            yield read_id, self.slice(ref_start, ref_end, [read_id])

    def iter_reads_db(self, reads, ref_bounds, full_overlap, max_reads, ignore_bam, return_tracks):
        if ref_bounds is not None:
            self.conf.tracks.ref_bounds = ref_bounds
        if reads is None:
            reads = self.prms.read_filter
        if max_reads is None:
            max_reads = self.prms.max_reads

        main_io = None #list()
        cmp_iters = list()
        cmp_ios = list()

        for io in self.inputs:
            if not isinstance(io, BAM):
                if main_io is not None:
                    raise RuntimeError("Tracks.iter_reads() only supports multiple inputs for BAM files")
            elif ignore_bam:
                continue
            if main_io is None:
                main_io = io
            else:
                cmp_iters.append(io.iter_sam())
                cmp_ios.append(io)

        if len(cmp_iters) == 0:
            for aln in main_io.iter_alns():
                if return_tracks:
                    yield aln.read_id, self._alns_to_tracks(aln)
                else:
                    yield aln.read_id, aln
            return

        cmp_starts = [None for _ in cmp_iters]
        sam_cache = [defaultdict(list) for _ in cmp_iters]

        for aln in main_io.iter_alns():
            aln_start = aln.sam.reference_start
            aln_end = aln.sam.reference_end
            aln_coord = (aln.sam.reference_id, aln.sam.reference_end, aln.sam.is_reverse)
            sams = [list() for _ in sam_cache]
            alns_out = [aln]
            for i in range(len(sam_cache)):
                while cmp_starts[i] is None or cmp_starts[i] < aln_coord:
                    sam = next(cmp_iters[i], None)
                    if sam is None:
                        break
                    sam_cache[i][sam.query_name].append(sam)
                    cmp_starts[i] = (sam.reference_id, sam.reference_start, sam.is_reverse)

                best = None
                for sam in sam_cache[i][aln.read_id]:
                    ovl = sam.get_overlap(aln_start, aln_end)
                    if ovl > 0 and (best is None or best[0] < ovl):
                        best = (ovl, sam)
                
                if best is None:
                    alns_out.append(None)
                else:
                    alns_out.append(cmp_ios[i].sam_to_aln(best[1]))

            if return_tracks:
                yield aln.read_id, self._alns_to_tracks(alns_out)
            else:
                yield aln.read_id, alns_out

    def _alns_to_tracks(self, alns):
        if isinstance(alns, Alignment):
            alns = [alns]

        ref_name = alns[0].seq.name
        aln_rows = list()
        layer_rows = list()
        min_pos = np.inf
        max_pos = -np.inf

        if len(self.cmp_layers) > 0:
            if len(alns) != 2:
                raise RuntimeError("Can only compare two alignments of the same read")
            alns[0].calc_dtwcmp(alns[1].instance)

        for a in alns:
            if a is None: continue
            assert(a.seq.name == ref_name)
            min_pos = min(min_pos, a.seq.coord.get_start())
            max_pos = max(max_pos, a.seq.coord.get_end())
            aln_rows.append(a.attrs())
            layer_rows.append(a.to_pandas(self.prms.layers, index=["aln.track", "aln.id", "seq.pos", "seq.fwd"], bounds=self.coords))

        aln_df = pd.DataFrame(aln_rows, columns=aln_rows[0]._fields).set_index(["track", "id"]).sort_index()
        layer_df = pd.concat(layer_rows).sort_index()

        coords = RefCoord(ref_name, min_pos, max_pos+1)
        return self._tables_to_tracks(coords, aln_df, layer_df)

        
    def _tables_to_tracks(self, coords, alignments, layers):
        tracks = dict()
        self.new_alignment = False
        self.new_layers = set()

        layer_alns = layers.index.get_level_values("aln.id")

        if self.prms.shared_refs_only or self.prms.min_coverage > 1:
            track_covs = layers.index.droplevel("aln.id")
            track_covs = track_covs.value_counts()

            mask = track_covs >= self.prms.min_coverage
            if not np.any(mask):
                idx = layers.index[:0]
            elif self.prms.shared_refs_only:
                track_counts = pd.MultiIndex.from_tuples(track_covs[mask].index) \
                                   .droplevel(0) \
                                   .value_counts()

                shared = track_counts.index[track_counts == len(self.alns)]
                if len(shared) == 0: 
                    idx = layers.index[:0]
                else:
                    idx = pd.MultiIndex.from_tuples(shared)
            else:
                idx = track_covs[mask].index.unique()

            l = layers.index.drop
            layers = layers[layers.index.droplevel(["aln.track","aln.id"]).isin(idx)] # .loc[(slice(None),idx.get_level_values(0),idx.get_level_values(1),slice(None))]
            layer_alns = layers.index.droplevel(["seq.fwd","seq.pos"])
            alignments = alignments.loc[layer_alns.unique()]

        for parent in self.alns:
            track = AlnTrack(parent, coords)#, track_alns, track_layers)

            tracks[parent.name] = track

        tracks = Tracks(self, coords)
        tracks.layers = layers 
        tracks.alignments = alignments

        tracks.calc_refstats()

        return tracks

    def close(self):
        for io in self.inputs:
            io.close()
        if self.output is not None and (len(self.inputs) == 0 or self.output != self.inputs[0]):
            self.output.close()
