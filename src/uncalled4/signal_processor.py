"""Stores datastructures to compute and store normalized events"""
import sys, os, argparse, re
from collections import defaultdict
import pandas as pd
import numpy as np

from .config import ParamGroup, Config
from . import EventDetector, Normalizer
from .pore_model import PoreModel
import _uncalled4

class ProcessedRead(_uncalled4._ProcessedRead):
    def __init__(self, read):
        if isinstance(read, pd.DataFrame):
            _uncalled4._ProcessedRead.__init__(self)
            if not "stdv" in read.columns:
                read["stdv"] = 0
            read = read[["mean","stdv","start","length"]].to_records(index=False)
            self.set_events(read)
        else:
            _uncalled4._ProcessedRead.__init__(self, read)

    def set_events(self, df):
        if isinstance(df, pd.DataFrame):
            if not "stdv" in df.columns:
                df["stdv"] = 0
            df = df[["mean","stdv","start","length"]].to_records(index=False)
        _uncalled4._ProcessedRead.set_events(self, df)

    def event_bounds(self, samp_start, samp_end):
        start = np.searchsorted(self.events["start"], samp_start)
        end = np.searchsorted(self.events["start"], samp_end)
        if start > 0: start -= 1
        return start, end

    def sample_range(self, start, end):
        evt_st,evt_en = self.event_bounds(start, end)
        return self.to_df()[evt_st:evt_en]

    def get_norm_signal(self, samp_min=0, samp_max=None):
        n = self.norm[0]
        if samp_max is None:
            samp_max = len(self.signal)
        norm = self.signal.to_numpy()[samp_min:samp_max]
        return norm

    def to_df(self):
        return pd.DataFrame(self.events)

class SignalProcessor:
    def __init__(self, model, conf):
        if isinstance(model.instance, _uncalled4.PoreModelU16):
            self.InstanceClass = _uncalled4.SignalProcessorU16
        elif isinstance(model.instance, _uncalled4.PoreModelU32):
            self.InstanceClass = _uncalled4.SignalProcessorU32
        else:
            raise ValueError(f"Unknown PoreModel type: {model.instance}")

        self.model = model

        if isinstance(model, PoreModel):
            model = model.instance

        self.instance = self.InstanceClass(model, conf.event_detector, conf.normalizer)

    def __getattr__(self, name):
        return self.instance.__getattribute__(name)

    def process(self, read, normalize=False):
        return ProcessedRead(self.instance.process(read, normalize))

    def process_signal(self, signal, normalize=False):
        read = _uncalled4.ReadBuffer("",0,0,0,signal)
        return ProcessedRead(self.instance.process(read, normalize))
