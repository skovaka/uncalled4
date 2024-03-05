from . import config
from .pore_model import PoreModel
from .ref_index import load_index, RefCoord, str_to_coord
from .read_index import ReadIndex, Fast5Reader, Slow5Reader

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
    #}, "moves" : {
    #    "start" : _Layer("Int32", "BC Sample Start", True),
    #    "length" : _Layer("Int32", "BC Sample Length", True),
    #    "end" : _Layer("Int32", "Sample End", False),
    #    "middle" : _Layer(np.float32, "Sample Middle", False),
    #    "indel" : _Layer("Int32", "Basecalled Alignment Indel", False),
    "dtwcmp" : {
        "aln_a" : _Layer("Int32", "Compare alignment V", True),
        "aln_b" : _Layer("Int32", "Compare alignment A", True),
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
        if layer in LAYER_META.index.get_level_values(0):
            group = layer
            _,layers = DEFAULT_LAYERS.get_loc_level(layer)
            for layer in layers:
                yield (group, layer)
            return
        
        matches = layer == LAYER_META.index.get_level_values(1)
        if np.sum(matches) > 1:
            raise ValueError(f"Ambiguous layer: {layer}. Please specify layer group.")
        elif np.sum(matches) == 1:
            for group,layer in LAYER_META.index[matches]:
                yield (group,layer)
            return

        group = "dtw"
    else:
        raise ValueError("Invalid layer: \"{layer}\"")

    if not (group, layer) in LAYER_META.index:
        raise ValueError(f"Invalid layer \"{group}.{layer}\"")

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

class Sequence:
    LAYERS = {"pos", "mpos", "pac", "name", "fwd", "strand", "kmer", "current", "base"}
    CONST_LAYERS = {"name", "fwd", "strand"}
    DEFAULT_LAYERS = ["pos", "kmer"]

    def __init__(self, seq, offset):
        self.instance = seq
        self.offset = offset
        self.index = self.instance.mpos

    @property
    def name(self):
        return self.coord.name

    @property
    def is_flipped(self):
        return self.index.start < 0

    @property
    def mpos(self):
        return self.index.expand().to_numpy()

    @property
    def pos(self):
        if self.is_flipped:
            return -self.mpos-1
        return self.mpos

    @property
    def pac(self):
        return self.offset + self.pos

    @property
    def strand(self):
        return "+" if self.fwd else "-"

    @property
    def base(self):
        return self.model.kmer_base(self.kmer, self.model.PRMS.shift)

    @property
    def fwd(self):
        return self.is_fwd

    def __len__(self):
        return len(self.instance)

    def _iter_layers(self, names):
        ret = list()

    def to_pandas(self, layers=None):
        if layers is None:
            layers = ["seq.pos", "kmer"]

        cols = dict()
        for name in layers:
            val = getattr(self, name)
            if name in self.CONST_LAYERS:
                val = np.full(len(self), val)
            cols[name] = val
        cols["index"] = self.mpos
        return pd.DataFrame(cols).set_index("index")

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise AttributeError(f"Sequence has no attribute '{name}'")
        return self.instance.__getattribute__(name)

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
        return 1000 * self.length / self.seq.model.PRMS.sample_rate
    
    def _get_series(self, name):
        vals = getattr(self, name, None)
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
        
        #idx = self.index.expand().to_numpy()
        #if index == "ref" and idx[0] < 0:
        #    idx = -idx-1
        #elif index != "mpos":
        #    raise ValueError(f"Unknown index column: {index}")
        #df[index] = idx

        return df.set_index("index")

    def __len__(self):
        return len(self.instance)

    def __getattr__(self, name):
        ret = getattr(self.instance, name, None)
        if ret is None:
            raise AttributeError(f"AlnDF has no attribute '{name}'")
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
        
        #idx = self.index.expand().to_numpy()
        #if index == "ref" and idx[0] < 0:
        #    idx = -idx-1
        #elif index != "mpos":
        #    raise ValueError(f"Unknown index column: {index}")
        #df[index] = idx
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
        #self.id = aln_id
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

        #Super = getattr(_uncalled4, f"_AlignmentK{seq.K}", None)
        #if Super is None:
        #    raise ValueError(f"Invalid k-mer length {seq.K}")
        self.instance = Super(aln_id, read_id, seq.instance)

        self.dtw = AlnDF(seq, self.instance._dtw)
        self.moves = AlnDF(seq, self.instance._moves)
        self.mvcmp = CmpDF(seq, self.instance._mvcmp)
        self.dtwcmp = CmpDF(seq, self.instance._dtwcmp)

    def __getattr__(self, name):
        if not hasattr(self.instance, name):
            raise AttributeError(f"Alignment has no attribute '{name}'")
        return self.instance.__getattribute__(name)
        
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
                #if idx is None:
                #    idx = vals[name].index
                #else:
                #    idx = vals[name].index.intersection(idx)

        df = pd.concat(vals, axis=1, names=["group", "layer"]).reset_index(drop=True)

        #df["aln","id"] = self.id
        
        if index is not None:
            df.set_index(list(index), inplace=True)
            if join_index:
                df.index.names = [".".join(idx) for idx in index]

        return df#.set_index(index)

    Attrs = namedtuple("Attrs", [
        "track", "id", "read_id", "ref_name", "ref_start", "ref_end", 
        "fwd", "sample_start", "sample_end", "coord"
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
            self.seq.fwd, samp_start, samp_end, self.seq.coord
        )

class AlnTrack:
    def __init__(self, *args, **kwargs):
        if isinstance(args[0], AlnTrack):
            self._init_slice(*args, **kwargs)
        else:
            self._init_new(*args, **kwargs)

        self.mat = None
        #self.coords = None
        #self.alignments = pd.DataFrame()
        #self.layers = pd.DataFrame()

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
        #self.set_data(coords, alignments, layers, order)

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


        #self.layers = self.layers.sort_index()

        self.has_fwd = np.any(self.alignments['fwd'])
        self.has_rev = not np.all(self.alignments['fwd'])

    @property
    def all_fwd(self):
        return np.all(self.alignments["fwd"])

    @property
    def all_rev(self):
        return not np.any(self.alignments["fwd"])

    #def slice(self, ref_start=0, ref_end=np.inf, aln_ids=None):
    def slice(self, coords=slice(None), aln_ids=None, reads=None, order=["fwd","ref_start"]):
        #ref_start = max(self.layer_refs.min(), ref_start)
        #ref_end = min(self.layer_refs.max()+1, ref_end)
        #coords = self.coords.ref_slice(ref_start, ref_end)
        if self.empty:
            return AlnTrack(self, coords, self.alignments, self.layers)
        layer_refs = self.layers.index.get_level_values("seq.pos")


        #layers = self.layers.loc[layer_refs.isin(coords.refs)]
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

        return AlnTrack(self, coords, alignments, layers, 
                        order=order) #TODO should handle order in tracks
        

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
        pd.set_option('display.max_rows', 50)


        if self.layers is None or overwrite:
            self.layers = df #pd.DataFrame({
            #    layer : df[layer].astype(LAYER_META.loc[layer,"dtype"]) 
            #    for layer in df}, columns=df.columns)
        else:
            #TODO don't always parse every layer twice
            #self.layers = pd.concat([self.layers, self._parse_layers(df)], axis=1)
            self.layers = pd.concat([self.layers, df], axis=1)
            #for layer in df:
            #    self.layers[layer] = df[layer]#.astype(LAYER_META.loc[layer,"dtype"])

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
                        #raise ValueError(f"Layer not found: {group}.{layer}")
                    vals = fn(self)
                    self.layers[group,layer] = vals

    def load_mat(self):
        df = self.layers.copy()
        #df["aln.id"] = df.index.get_level_values("aln.id")
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

        #cmp_df["aln_b"] = alns_b[0]
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
    "cov" : len,
    "mean" : np.mean, 
    "median" : np.median, 
    "q5" : (lambda x: np.quantile(x, 0.05)),
    "q95" : (lambda x: np.quantile(x, 0.95)),
    "q25" : (lambda x: np.quantile(x, 0.25)),
    "q75" : (lambda x: np.quantile(x, 0.75)),
    "stdv" : np.std, 
    "var"  : np.var,
    "skew" : scipy.stats.skew,
    "kurt" : scipy.stats.kurtosis,
    "min"  : np.min, 
    "max"  : np.min,
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

CMP_GROUPS = {"cmp", "mvcmp"}

class RefstatsSplit:
    def __init__(self, stats, track_count):
        if not "cov" in stats:
            stats = ["cov"] + stats
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
        else:
            if len(pm.name) == 0 and pm.has_workflow():
                pm.name = PoreModel.PRESET_MAP.loc[pm.get_workflow(), "preset_model"]
                sys.stderr.write(f"Auto-detected flowcell='{pm.flowcell}' kit='{pm.kit}'\n")
            self.model = None


        self._init_io()

        if self.prms.ref_index is None:
            raise RuntimeError("Failed to load reference index")

        if self.model is None:
            raise RuntimeError("Failed to detect model, please specify --flowcell and --kit")

        self.set_layers(self.prms.layers)


        self._aln_track_ids = [t.id for t in self.alns]

        self.index = load_index(self.model, self.prms.ref_index)
        self.refstats = None

        #self.coords = self._ref_bounds_to_coords(self.prms.ref_bounds)
        self.coords = self.prms.ref_bounds

        if self.coords is not None:
            if not self.coords.has_bounds:
                self.coords.start = 0
                self.coords.end = self.index.get_ref_len(self.coords.name)

            if len(self._aln_track_ids) > 0:
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
                self.layers = parent.layers.reset_index(level="seq.pos").loc[self.alignments.index].set_index("seq.pos",append=True)
            else:
                self.alignments = parent.alignments
                self.layers = parent.layers
        else:
            self.alignments = None
            self.layers = None
    
        #if self.alignments is not None:
        #    self.init_mat()

    def _add_tracks(self, tracks):
        for name,track in tracks.items():
            self._add_track(name, track)

    def _add_track(self, name, track):
        #if name in self._tracks:
        #    raise KeyError(f"Duplicate track name: {name}")
        self._tracks[name] = track
        
        if name in BUILTIN_TRACKS:
            setattr(self, name[1:], track)
        elif isinstance(track, AlnTrack):
        #elif isinstance(track, tuple):
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
            #for track in self.output.tracks:
            #    self.output_tracks[track.name] = track
        else:
            self.output = None
            self.output_track = self.inputs[0].aln_tracks[0].name

            #TODO each IO should have reference to its AlnTracks?
            #probably construct here in Tracks? or maybe in IO!
            #iterate through each input, it will populate its own AlnTracks

            #elif track.model is not None and self.model.K != track.model.K:
                #raise ValueError("Cannot load tracks with multiple k-mer lengths (found K={self.model.K} and K={track.model.K}")

        #pm = self.conf.pore_model
        #has_flowkit = not (len(pm.flowcell) == 0 or len(pm.kit) == 0)
        #if not has_flowkit:
        #    flowcell, kit = self.read_index.default_model
        #    if not (kit is None or flowcell is None):
        #        if len(pm.flowcell) == 0:
        #            pm.flowcell = flowcell
        #        if len(pm.kit) == 0:
        #            pm.kit = kit
        #        has_flowkit = True
        ##PoreModel.PRESET_MAP.loc[(pm.flowcell, pm.kit), "preset_model"]

        #if self.model is None and has_flowkit:
        #    if len(self.conf.pore_model.name) == 0:
        #        self.model = PoreModel(PoreModel.PRESET_MAP.loc[(pm.flowcell, pm.kit), "preset_model"])
        #    else:
        #        self.model = PoreModel(params=self.conf.pore_model)
        #    for track in self.alns:
        #        track.model = self.model

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
        #ret = pd.Index([])
        #for track in self.alns:
        #    layers = track.layers.columns
        #    if layer_filter is not None:
        #        layers = layers.intersection(layer_filter)
        #    ret = ret.union(layers)
        #return ret
    
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
            #dtw["stdv"] = np.std(read.signal[dtw["start"]:dtw["start"]+dtw["length"]], axis=1)
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


    def init_alignment(self, track_name, aln_id, read, ref_id, coords, sam=None):
        #TODO make model name paramter, initialize with track specific model
        track = self._track_or_default(track_name)
        #seq = self.index.instance.get_kmers(model, ref_id, coords, self.conf.is_rna)
        try:
            seq = self.index.query(coords)
        except RuntimeError:
            raise RuntimeError(f"Invalid coordinate for {read.id}: {coords}")
        #except:
        #    sys.exit(1)
        seq = Sequence(seq, self.index.get_pac_offset(ref_id))
        #except:
        #    raise RuntimeError(f"Read {read.id} fiale")
        return Alignment(aln_id, read, seq, sam, track_name)

    def write_alignment(self, aln):
        if self.output is not None:
            out = self.output
        else:
            out = self.inputs[0]

        #sd = aln.dtw.current_sd.to_numpy()
        #mask = (sd >= 0) & (sd <= self.prms.max_sd)
        #aln.dtw.mask(mask)

        out.write_alignment(aln)

    def set_read(self, read):
        self.output.read = read

    @property
    def input_count(self):
        return len(self.alns)

    def _verify_read(self):
        if len(self._aln_track_ids) == 0:
            raise RuntimeError("No input tracks have been loaded")

    #TODO read_ids, track_alns, max_cov(?)
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
                
        #if tracks is None:
        #    track_names = set(self.track_names)
        #else:
        #    track_names = set(tracks)
        #tracks = dict()
        #for name,track in self._tracks.items():
        #    if name not in track_names: continue
        #    if isinstance(track, pd.DataFrame):
        #        tracks[name] = track.loc[coords.refs]

        #    elif isinstance(track, AlnTrack):
        #        tracks[name] = track #.slice(coords, reads=reads, order=order)

        return Tracks(self, coords, reads)

    def get_shared_reads(self):
        all_ids = self.alignments["read_id"]
        read_ids = pd.Index(all_ids.loc[self.alns[0].id])
        for track in self.alns[1:]:
            read_ids = read_ids.intersection(all_ids.loc[track.name])
        return read_ids

    def get_all_reads(self):
        #all_ids = self.alignments["read_id"]
        #read_ids = pd.Index(all_ids.loc[self.alns[0].name])
        #for track in self.alns[1:]:
        #    if track.name in all_ids.index:
        #        read_ids = read_ids.union(all_ids.loc[track.name])
        #return read_ids
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
            
    def load(self, ref_bounds=None, full_overlap=None, read_filter=None, load_mat=False):
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

        if load_mat is None:
            load_mat = self.prms.load_mat

        track_alns = dict()
        track_layers = dict()

        for io in self.inputs:
            alns, layers = io.query(self.conf.tracks.layers, self.coords, ["aln.id","seq.pos"], full_overlap=full_overlap, read_id=read_filter)


            track_alns.update(alns)
            track_layers.update(layers)

        #TODO eliminate need for AlnTrack
        #just use self.layers, self.aln_layers
        #will need to initialize Alignment from these dataframes for Dotplot
        self.layers = pd.concat(track_layers, axis=0, names=["aln.track", "aln.id", "seq.pos"])
        self.alignments = pd.concat(track_alns, axis=0, names=["track", "id"]).sort_index()#.sort_values(["fwd","ref_start"])
    
        #self.init_mat()
            
        #self.load_compare(alignments.index.droplevel(0).to_numpy())
        self.calc_refstats()

        for track in self.alns:
            track.calc_layers(self.fn_layers)
            if load_mat:
                track.load_mat()

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

    def load_compare(self, aln_ids=None):
        if len(self.cmp_layers) == 0:
            return

        #TODO handle multiple inputs properly
        io = self.inputs[0]

        self.cmp = io.query_compare(self.cmp_layers, self._aln_track_ids, self.coords, aln_ids)


        if self.cmp is None:
            return

        groups = self.cmp.index.get_level_values("group_b").unique()
        if "moves" in groups:
            movess = self.cmp.loc[(slice(None), slice(None), slice(None), "moves"),:]
        else:
            movess = None
        if "dtw" in groups:
            dtws = self.cmp.loc[(slice(None), slice(None), slice(None), "dtw"),:]
        else:
            dtws = None

        for track in self.alns:
            def _add_group(group, df):
                df = df.reset_index(["aln_b", "group_b"])
                df = df[df.index.get_level_values("seq.pos").isin(track.layer_pacs)]
                df.rename(index=track.coords.pac_to_pos, level=0, inplace=True)
                df.index.names = ["seq.pos", "aln.id"]
                df = pd.concat({group : df.reindex(track.layers.index)}, axis=1)
                track.layers = pd.concat([track.layers, df], axis=1).dropna(axis=1,how="all")

            #try:
            if movess is not None:
                _add_group("mvcmp", movess)
            if dtws is not None:
                _add_group("cmp", dtws)
            #except:
            #    sys.stderr.write("Failed to write compare group\n")
            #    sys.stderr.write(str(track.alignments))

    def calc_compare(self, group_b, single_track, save):
        if len(self.alns) > 0:
            alns = self.alns
        else:
            raise ValueError("Must input at least one track")
        
        if self.output_track is not None:
            track_a = self._tracks[self.output_track]
        elif single_track or len(self.alns) == 1:
            track_a = self.alns[0]
        else:
            track_a = self.alns[0]

        if single_track or len(self.alns) == 1:
            track_b = track_a
        else:
            for track_b in self.alns:
                if track_b != track_a: break

        cols = track_a.layers.columns.get_level_values("group").unique()
        if (group_b == "dtw" and "cmp" in cols) or (group_b == "moves" and "mvcmp" in cols):
            sys.stderr.write(f"Read already has compare group. Skipping\n")
            return None

        if group_b == "dtw":
            if track_a == track_b:
                raise ValueError("Must input exactly two tracks to compare dtw alignments")

            df = track_a.cmp(track_b, True, True)

        elif group_b == "moves":
            df = track_a.mvcmp(track_b, True, True)

        df = df.dropna(how="all")
        self.add_layers("cmp", df)

    def calc_refstats(self, cov=False):
        if self.prms.refstats is None or len(self.prms.refstats) == 0 or len(self.prms.refstats_layers) == 0 or self.alignments is None:
            self.refstats = None
            return None

        stats = RefstatsSplit(self.prms.refstats, len(self.alns))

        refstats = dict()

        groups = self.layers[self.prms.refstats_layers].groupby(level=["aln.track", "seq.pos", "seq.fwd"])

        refstats = groups.agg(stats.layer_agg).reset_index().pivot(index=["seq.pos","seq.fwd"], columns="aln.track")
        #rename = ({
        #    old[-1] : new
        #    for old,new in zip(refstats.columns, stats.layer)
        #})
        #refstats.rename(columns=rename, inplace=True)
        #if cov:
        #    refstats[track.name].insert(0, "cov", groups.size())

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
                        vals = scipy.stats.stats.ks_2samp(df.loc[tids[0],col],df.loc[tids[1],col],mode="asymp")
                        d[col+("ks",)] = pd.Series(vals, index=["stat","pval"], name=df.name)
                    df = pd.concat(d, axis=0)
                    return df

            cmp_groups = self.layers[self.prms.refstats_layers].groupby(level=["seq.pos","seq.fwd"])
            df = cmp_groups.apply(ks)
            refstats = pd.concat([refstats, df], axis=1)

        self.refstats = refstats.dropna()

        #self.refstats["seq.kmer"] = self.index.get

        self._tracks["_refstats"] = self.refstats

        return refstats

    def iter_refs(self, ref_bounds=None):
        layer_iters = [
            io.iter_refs(
                self.db_layers, 
                #coords=coords, 
                chunksize=self.prms.io.ref_chunksize)
            for io in self.inputs]
        
        #for alns,layers in layer_iters[0]:
        #    layers = pd.concat({1 : layers}, names=["track.name","seq.fwd","seq.pac","aln.id"])
        #    alns = pd.concat({1 : alns}, names=["track.name","aln.id"])
        #    refs = layers.index.get_level_values("seq.pac")
        #    coords = RefCoord(alns.iloc[0]["ref_name"], refs.min(), refs.max()+1)
        #    yield self._tables_to_tracks(coords, alns, layers)
        #return

        #chunks = [next(i) for i in layer_iters]
        chunks = [(pd.DataFrame([]), pd.DataFrame([])) for i in layer_iters]
        chunk_hasnext = np.ones(len(chunks), bool)

        #while np.any([len(chunk[0]) > 0 for chunk in chunks]):
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

            #np.argmin([l.index.get_level_values("seq.pac").max() for l,c in chunks])

            ##TODO simplify! assume inputs segment by ref coord
            #for i in range(len(chunks)):
            #    chunk_alignments, chunk_layers = chunks[i]
            #    last_chunk = False

            #    pacs = chunk_layers.index.get_level_values("seq.pac").unique()
            #    while len(pacs) <= 1:
            #        alns, layers = next(layer_iters[i], (None, None))
            #        if alns is None:
            #            chunk_hasnext[i] = False
            #            break

            #        chunk_alignments = pd.concat([chunks[i][0], alns])
            #        chunk_alignments = chunk_alignments[~chunk_alignments.index.duplicated()]

            #        chunk_layers = pd.concat([chunk_layers, layers])
            #        chunks[i] = (chunk_alignments,chunk_layers)
            #        pacs = chunk_layers.index.get_level_values("seq.pac").unique()

            #    #idx = 
            #    pacs = chunk_layers.index.unique("seq.pac")

            #    all_pacs = all_pacs.union(pacs)
            #    pac_start = min(pac_start, pacs.min())

            #    #TODO this part is still key. only output up until the end of the least full chunk
            #    pac_end = min(pac_end, pacs.max())+1


            #layers = pd.concat({(i+1) : l for i,(a,l) in enumerate(chunks)}, names=["track.name","seq.fwd","seq.pac","aln.id"])
            #alns = pd.concat({(i+1) : a.loc[l.index.droplevel(["seq.fwd","seq.pac"]).unique()] for i,(a,l) in enumerate(chunks)}, names=["track.name","aln.id"])
            #aln_ids = alns.index

            #refs = layers.index.get_level_values("seq.pac")
            #coords = RefCoord(alns.iloc[0]["ref_name"], refs.min(), refs.max()+1)

            #ret = self._tables_to_tracks(coords, alns, layers)

            ##leftover collection
            #for i in range(len(chunks)):
            #    chunk_alns, chunk_layers = chunks[i]

            #    chunk_alns = chunk_alns.iloc[:0]
            #    chunk_layers = chunk_layers.iloc[:0]
            #    chunks[i] = (chunk_alns,chunk_layers)

            yield ret

    def iter_reads(self, read_filter=None, ref_bounds=None, full_overlap=False, max_reads=None, ignore_bam=False, return_tracks=False):

        if read_filter is None and self.read_index is not None:
            read_filter = self.read_index.read_filter
        
        if ref_bounds is not None and not isinstance(ref_bounds, RefCoord):
            ref_bounds = RefCoord(ref_bounds)

        if (self.coords is None or
            (read_filter is not None and 
             len(self.get_all_reads().intersection(read_filter)) < len(read_filter)) or
            (ref_bounds is not None and not self.coords.contains(ref_bounds))):
            gen = self.iter_reads_db(read_filter, ref_bounds, full_overlap, max_reads, ignore_bam, return_tracks)
        else:
            gen = self.iter_reads_slice(read_filter, ref_bounds)

        return gen
        #for read_id,chunk in gen:
        #    yield read_id,chunk
            
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
            self._set_ref_bounds(ref_bounds)
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
                #coords = aln.seq.coord#.name, pos.min(), pos.max()+1)
                #yield aln.read_id, aln
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
                    #alns_out.append(main_io.sam_to_aln(best[1],  load_moves=False))
                    alns_out.append(cmp_ios[i].sam_to_aln(best[1]))

            if return_tracks:
                yield aln.read_id, self._alns_to_tracks(alns_out)
            else:
                yield aln.read_id, alns_out

        #TODO handle multiple inputs properly
        #for io in self.inputs:
        #    if ignore_bam and io == self.bam_in:
        #        continue

        #    aln_iter = io.iter_alns(
        #        self.db_layers, 
        #        self._aln_track_ids,
        #        self.coords,
        #        read_id=reads,
        #        full_overlap=full_overlap,
        #        ref_index=self.index)

        #    for aln in aln_iter:
        #        yield aln

    def _alns_to_tracks(self, alns):
        if isinstance(alns, Alignment):
            alns = [alns]

        ref_name = alns[0].seq.name
        #aln_rows = defaultdict(list)
        #layer_rows = defaultdict(list)
        aln_rows = list()
        layer_rows = list()
        min_pos = np.inf
        max_pos = -np.inf

        for a in alns:
            if a is None: continue
            assert(a.seq.name == ref_name)
            min_pos = min(min_pos, a.seq.coord.get_start())
            max_pos = max(max_pos, a.seq.coord.get_end())
            aln_rows.append(a.attrs())
            layer_rows.append(a.to_pandas(self.prms.layers, index=["aln.track", "aln.id", "seq.pos", "seq.fwd"]), self.coords)

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
            #track_covs = alignments.loc[layer_alns, ["track.name"]] 
            #track_covs = track_covs.set_index(layers.index)
            #track_covs = track_covs.reset_index("aln.id", drop=True) 
            #track_covs = track_covs.set_index("track.name", append=True) 
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

        aln_groups = alignments.index.unique("aln.track")
        for parent in self.alns:
            if parent.id in aln_groups:
                track_alns = alignments.loc[parent.id]
                track_layers = layers.loc[parent.id].droplevel("seq.fwd")
            else:
                track_alns = alignments.iloc[:0] 
                track_layers = layers.iloc[:0]   

            track = AlnTrack(parent, coords)#, track_alns, track_layers)

            #if not track.empty:
            #track.calc_layers(self.fn_layers)

            tracks[parent.name] = track

        tracks = Tracks(self, coords)
        tracks.layers = layers #, axis=0, names=["track.name", "aln.id", "seq.pos"])
        tracks.alignments = alignments #, axis=0, names=["track", "id"]).sort_index()#.sort_values(["fwd","ref_start"])

        #if not tracks.all_empty:
        #    tracks.load_compare(alignments.index.to_numpy())
        tracks.calc_refstats()

        return tracks

    def close(self):
        for io in self.inputs:
            io.close()
        if self.output is not None and (len(self.inputs) == 0 or self.output != self.inputs[0]):
            self.output.close()
