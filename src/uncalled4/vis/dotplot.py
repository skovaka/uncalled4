"""Plot signal-to-reference alignment dotplots"""

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import time
import sys

from .sigplot import Sigplot

from .. import config
from ..ref_index import str_to_coord
from ..tracks import Tracks, LAYER_META, parse_layers
from ..argparse import Opt, comma_split

class Dotplot:

    REQ_LAYERS = [
        "dtw.start", "dtw.length", "dtw.middle_sec", 
        "dtw.current", "dtw.current_sd", "dtw.dwell", "seq.kmer",  
        "moves.middle_sec", "seq.current",
        "dtw.start_sec", "dtw.length_sec"
    ]

    def __init__(self, *args, **kwargs):
        self.conf, self.prms = config._init_group("dotplot", *args, **kwargs)

        req_layers = self.REQ_LAYERS.copy()

        self.conf.tracks.layers = req_layers + self.prms.layers

        if isinstance(self.prms.tracks, str) or self.prms.tracks is None:
            self.tracks = Tracks(conf=self.conf)
        elif isinstance(self.prms.tracks, Tracks):
            self.tracks = self.prms.tracks
        else:
            raise ValueError("Dotplot tracks parameter must be string or Tracks instance")

        self.layers = list(parse_layers(self.prms.layers))

        self.tracks.set_layers(req_layers + self.layers)

        self.conf.load_config(self.tracks.conf)

        self.Scatter = go.Scatter if self.conf.vis.svg else go.Scattergl

        self.fig_config = {
                "toImageButtonOptions" : {
                    "format" : "svg" if self.conf.vis.svg else "png", 
                    "width" : None, "height" : None },
                "scrollZoom" : True, 
                "displayModeBar" : True
        }

    def iter_plots(self):
        t0 = time.time()
        #for read_id, aln in self.tracks.iter_reads(return_tracks=True):
        for read_id, tracks in self.tracks.iter_reads(return_tracks=True):
            yield read_id, self._plot(read_id, tracks)
            t0 = time.time()

    def plot(self, read_id):
        chunk = self.tracks.slice(reads=[read_id])
        return self._plot(read_id, chunk)

    def show(self, read_id):
        self.plot(read_id).show(config=self.fig_config)

    def _plot(self, read_id, tracks):
        cmp_stats = []

        column_widths=[6]+[1]*(len(self.layers)+len(cmp_stats))

        legend = set()

        fig = make_subplots(
            rows=2, cols=len(column_widths), 
            row_heights=[1,3],
            column_widths=column_widths,
            vertical_spacing=0.01,
            horizontal_spacing=0.01,
            shared_xaxes=True,
            shared_yaxes=True)

        tracks_filter,colors_filter = zip(*[(t,c) for t,c in zip(tracks.alns,self.conf.vis.track_colors) if t.name != self.prms.moves_track])

        Sigplot(tracks, conf=self.conf).plot(fig)

        hover_layers = [("dtw", "middle_sec"),("seq","kmer"),("dtw","current"),("dtw","dwell")] + self.layers
        hover_data = dict()

        coords = tracks.coords

        flipped = True
        fwd = True
        moves_plotted = False

        active_tracks = tracks.alignments.index.unique("track")

        for i,track in enumerate(tracks.alns):
            if track.name not in active_tracks: continue

            track_hover = list()

            has_moves = "moves" in tracks.layers.columns.get_level_values("group")
            only_moves = self.prms.moves_track == track.name
            model = track.model

            first_aln = True
            for aln_id, aln in tracks.alignments.loc[track.name].iterrows():
                layers = tracks.layers \
                              .loc[(track.name,aln_id)] \
                              .reset_index() \
                              .set_index("seq.pos") \
                              .sort_values(("dtw","start_sec"))
                              #droplevel("seq.fwd") #,slice(None)), slice(None)] \
                              #.droplevel("aln.id")

                end = layers.iloc[[-1]].copy()
                end["dtw","start_sec"] += end["dtw","length_sec"]
                layers = pd.concat([ layers,  end ])
                
                fwd = fwd and aln["fwd"]
                flipped = flipped and aln["fwd"] == self.conf.is_rna

                if self.prms.show_bands and "band" in layers:
                    bands = layers["band"].dropna(how="any")
                    fig.add_trace(self.Scatter(
                        x=bands["sample_start"], 
                        y=bands.index,
                        line={"color" : "orange"},
                        fillcolor="orange",
                        opacity=0.2,
                        fill="tonexty",
                        name="DTW Band",
                        legendgroup="band",
                    ), row=2, col=1)
                    fig.add_trace(self.Scatter(
                        x=bands["sample_end"], 
                        y=bands["ref_end"],
                        line={"color" : "orange"},
                        fillcolor="orange",
                        opacity=0.2,
                        fill="tonexty",
                        legendgroup="band",
                        showlegend=False
                    ), row=2, col=1)
                
                if has_moves and not moves_plotted:
                    moves_plotted = True
                    self._plot_moves(fig, legend, layers)

                if not only_moves: 
                    #end = layers.iloc[[-1]].copy()
                    #end["dtw","start_sec"] += end["dtw","length_sec"]
                    #df = pd.concat([ layers,  end ])
                    fig.add_trace(self.Scatter(
                        x=layers["dtw","start_sec"], y=layers.index,
                        name=track.name,
                        legendgroup=track.name,
                        mode="lines",
                        line={
                            "color":self.conf.vis.track_colors[i], 
                            "width":2, "shape" : "hv" }, #if not flipped else "vh" },
                        hoverinfo="skip",
                        showlegend=first_aln
                    ), row=2, col=1)
                if only_moves: continue

                track_hover.append(layers[hover_layers])

                first_aln = False

                for j,layer in enumerate(self.layers):
                    if not layer in layers.columns: continue
                    df = layers[layer].sort_index()
                    if layer[0] != "cmp":
                        fig.add_trace(self.Scatter(
                            x=df, y=df.index+0.5,
                            name=track.name, 
                            mode="lines",
                            line={
                                "color" : self.conf.vis.track_colors[i], 
                                "width":2, "shape" : "hv"},
                            legendgroup=track.name, showlegend=False,
                        ), row=2, col=j+2)



                    elif len(df.dropna()) > 0:
                        color= "rgba(255,0,0,1)"
                        fig.add_trace(go.Bar(
                            x=df.fillna(1.0), 
                            y=track.layer_refs, #TODO try vhv
                            base=0,
                            name="DTW Compare",
                            orientation="h",
                            width=1,
                            marker={"color":color,"line":{"color":color,"width":0.5}},
                            legendgroup="cmp",
                        ), row=2, col=j+2)


            if len(track_hover) > 0:
                hover_data[track.name] = pd.concat(track_hover)#.reset_index()
                hover_data[track.name] = track_hover[0]#.reset_index()

        if len(hover_data) > 0:
            hover_data = pd.concat(hover_data, axis=1)
            hover_coords = hover_data.xs("middle_sec", axis=1, level=2).mean(axis=1)

            hover_kmers = model.kmer_to_str(
                hover_data.xs("kmer", 1, 2)
                          .fillna(method="pad", axis=1)
                          .iloc[:,-1])

            customdata = hover_data.drop(["kmer","middle_sec"], axis=1, level=2).to_numpy()

            hover_rows = [
                "<b>" + coords.name + ":%{y:,d} [%{text}]</b>"
            ]
            labels = [LAYER_META.loc[(g,l),"label"] for g,l in hover_layers[2:]]

            for i,label in enumerate(labels):
                s = label
                fields = list()
                for j in range(len(tracks_filter)):
                    k = len(labels) * j + i
                    fields.append(
                        '<span style="color:%s;float:right"><b>%%{customdata[%d]:.2f}</b></span>' % 
                        (self.conf.vis.track_colors[j], k))
                hover_rows.append(s + ": " + ", ".join(fields))


            fig.add_trace(self.Scatter(
                x=hover_coords, y=hover_data.index,
                mode="markers", marker={"size":0,"color":"rgba(0,0,0,0)"},
                name="",
                customdata=customdata,
                hovertemplate="<br>".join(hover_rows),
                hoverlabel={"bgcolor":"rgba(255,255,255,1)"},
                text=hover_kmers,
                showlegend=False
            ), row=2,col=1)

        if flipped:
            fig.update_yaxes(autorange="reversed", row=2, col=1)
            fig.update_yaxes(autorange="reversed", row=2, col=2)

        if self.prms.select_ref is not None:
            fig.add_hline(y=self.prms.select_ref, line_color="red", row=2, col=1, opacity=0.5)
            i = hover_data.index.get_loc(self.prms.select_ref)

        strand = "+" if fwd else "-"
        fig.update_yaxes(row=2, col=1,
            title_text=tracks.coords.name + f" ({strand})")

        for i,(group,layer) in enumerate(self.layers):
            fig.update_xaxes(row=2, col=i+2,
                title_text=LAYER_META.loc[(group,layer),"label"])

        axis_kw = dict(
            showspikes=True,
            spikemode="across",
            spikecolor="darkgray",
            spikethickness=1)

        fig.update_xaxes(**axis_kw)
        fig.update_yaxes(**axis_kw)
        fig.update_yaxes(
            tickformat=",d", tickangle=-90, nticks=3
        )

        fig.update_xaxes(
            title_text="Time (s)", 
            row=2, col=1)

        fig.update_layout(
            margin={"l":100,"r":50},#, "b":50},
            barmode="overlay",
            hoverdistance=20,
            dragmode="pan", 
            showlegend=self.prms.show_legend,
            legend={
                "y" : 1.05, "yanchor" : "bottom",
                "traceorder" : "normal",
                "orientation" : "h", "bgcolor" : "#e6edf6"})

        return fig

    def _plot_moves(self, fig, legend, layers):
        fig.add_trace(self.Scatter(
            x=layers["moves","middle_sec"], y=layers.index,
            name="Basecalled Alignment",
            mode="markers", marker={"size":5,"color":"orange"},
            legendgroup="moves",
            hoverinfo="skip",
            showlegend="moves_starts" not in legend
        ), row=2, col=1)
        legend.add("moves_starts")



def dotplot(conf):
    """Plot signal-to-reference alignment dotplots"""

    dotplots = Dotplot(conf=conf)
    save = conf.out_prefix is not None
    for read_id, fig in dotplots.iter_plots():
        if save:
            fig.write_html(
                conf.out_prefix + read_id + ".html", 
                config=dotplots.fig_config)
        else:
            fig.show(config=dotplots.fig_config)
            sys.stderr.flush()
            sys.stdout.write("Press enter to plot next read, or type \"exit\"\n")
            sys.stdout.flush()
            choice = sys.stdin.readline().strip().lower()
            if choice == "exit": break
