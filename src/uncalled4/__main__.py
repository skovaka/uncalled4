from .argparse import ArgParser, Opt, MutexOpts, CONFIG_PARAM, comma_split, ref_coords
from .ref_index import RefCoord, str_to_coord

import sys
import textwrap
import importlib
import os

import pandas as pd

import cProfile

def parse_read_ids(reads):
    if reads is None:
        return []

    if isinstance(reads, str):
        if os.path.exists(reads):
            with open(reads) as reads_in:
                return [line.split()[0] for line in reads_in]
        else:
            return reads.split(",")

    return list(reads)

CONFIG_OPT = Opt(("-C", "--config"), type=str, default=None, required=False, help="Configuration file in TOML format", dest = CONFIG_PARAM)

DTW_OPTS = (
    #Opt("--ref", "tracks"), 
    MutexOpts("mode", [
        Opt("--ref", "tracks"),
        Opt("--self", "tracks", action="store_true"),
    ]),
    Opt("--reads", "read_index", "paths", nargs="+", type=str),
    Opt("--bam-in", "tracks.io", nargs="?", const="-", required=True, help="BAM file used to guide signal alignment. Must contain 'mv' and 'ts' basecaller tags."),
    
    Opt(("-p", "--processes"), "tracks.io"),

    Opt("--flowcell", "pore_model"),
    Opt("--kit", "pore_model"),
    Opt("--basecaller-profile", "tracks"),
    Opt("--rna", fn="set_r94_rna", help="RNA alignment (required for custom pore models)"),

    Opt("--ordered-out", "tracks.io", action="store_true"),
    Opt(("-f", "--overwrite"), "tracks.io", action="store_true"),
    Opt(("--kmer-shift"), "pore_model", "shift", default=None),
    Opt("--bam-chunksize", "tracks.io"),
    Opt(("--out-name"), "tracks.io"),

    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt(("-l", "--read-filter"), "read_index"),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-n", "--max-reads"), "tracks"),

    Opt("--count-events", "tracks", action="store_true"),
    Opt("--del-max", "dtw"),
    Opt("--ins-max", "dtw"),
    Opt("--unmask-splice", "dtw", action="store_true"),
    Opt("--method", "dtw"),
    Opt(("-c", "--cost-fn"), "dtw", choices=["abs_diff","z_score","norm_pdf"]),
    Opt(("-b", "--band-width"), "dtw"),
    Opt(("-s", "--band-shift"), "dtw"),
    Opt("--mvcmp-mask", "tracks"),
    Opt("--max-norm-dist", "tracks"),
    Opt("--max-sd", "tracks"),
    Opt("--min-aln-length", "tracks"),
    Opt(("-N", "--norm-mode"), "normalizer", "mode", choices=["ref_mom", "model_mom"]),
    Opt("--zero-ts", "tracks", action="store_true"),
    CONFIG_OPT,
)

CONVERT_OPTS = (
    Opt("--ref", "tracks", nargs="?"),
    Opt(("-p", "--processes"), "tracks.io"),
    MutexOpts("input", [
        Opt("--eventalign-in", "tracks.io", type=comma_split, nargs="?", const="-"),
        Opt("--tombo-in", "tracks.io", type=comma_split, action="extend"),
    ]),
    Opt("--bam-in", "tracks.io", type=comma_split, action="extend"),
    Opt(("-t", "--tracks"), "tracks.io", "input_names", type=comma_split),
    Opt(("-m", "--pore-model"), "pore_model", "name"),
    Opt("--kmer-shift", "pore_model", "shift"),
    Opt("--bam-chunksize", "tracks.io"),
    Opt("--max-sd", "tracks"),
    Opt("--basecaller-profile", "tracks"),

    MutexOpts("output", [
        Opt("--eventalign-out", "tracks.io", nargs="?", const="-"),
        Opt("--tsv-out", "tracks.io", nargs="?", const="-"),
        Opt("--bam-out", "tracks.io", nargs="?", const="-"),
        Opt("--m6anet-out", "tracks.io"),
    ]),
    Opt(("--out-name", "-o"), "tracks.io"),

    Opt("--tsv-cols", "tracks.io", type=comma_split, default="dtw"),
    Opt("--tsv-na", "tracks.io", nargs="?", const="-"),
    Opt("--tsv-noref", "tracks.io", action="store_true"),

    Opt("--eventalign-flags", "tracks.io", type=comma_split),
    Opt("--mask-skips", "tracks", nargs="?", const="all"),

    Opt("--flowcell", "pore_model"),
    Opt("--kit", "pore_model"),
    Opt("--reads", "read_index", "paths", nargs="+", type=str),
    Opt(("-l", "--read-filter"), "read_index"),
    Opt(("-x", "--read-index"), "read_index", required=False),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),
    #Opt(("-R", "--ref-bounds"), "tracks", type=str_to_coord),
    Opt(("-R", "--region"), "tracks", "ref_bounds", type=RefCoord),
    Opt(("-f", "--overwrite"), "tracks.io", action="store_true"),
    Opt(("-a", "--append"), "tracks.io", action="store_true"),
    CONFIG_OPT,
)

COMPARE_OPTS = (
    Opt("--bam-in", "tracks.io", nargs="+"), #, required=True
    Opt(("-t", "--tracks"), "tracks.io", "input_names", type=comma_split),
    Opt(("-l", "--read-filter"), "tracks", type=parse_read_ids),
    Opt(("-R", "--region"), "tracks", "ref_bounds", type=RefCoord),
    Opt(("-m", "--moves"), action="store_true", help="Compare against basecalled alignment. If two tracks input will look for \"moves\" group in second track, otherwise will look in the first track."),
    Opt("--tsv-cols", "tracks.io", type=comma_split, default="dtwcmp,mvcmp"),
    Opt("--tsv-na", "tracks.io", nargs="?", const="-"),
    Opt("--tsv-noref", "tracks.io", action="store_true"),
    Opt(("--tsv-out", "-o"), "tracks.io", nargs="?", default="-"),
    CONFIG_OPT,
)

ALL_REFSTATS = {"min", "max", "mean", "median", "stdv", "var", "skew", "kurt", "q25", "q75", "q5", "q95", "KS"}
REFSTATS_OPTS = (
    Opt("--layers", "tracks", nargs="+",
        help="Comma-separated list of layers over which to compute summary statistics"),# {%s}" % ",".join(LAYERS.keys())),
    Opt("--stats", dest="refstats", nargs="+",
        help="Comma-separated list of summary statistics to compute. Some statisitcs (ks) can only be used if exactly two tracks are provided {%s}" % ",".join(ALL_REFSTATS)),
    Opt("--bam-in", "tracks.io", nargs="+"), #, required=True
    Opt(("-t", "--tracks"), "tracks.io", "input_names", type=comma_split),
    Opt(("-R", "--region"), "tracks", "ref_bounds", type=RefCoord),
    Opt("--min-coverage", "tracks"),
    Opt("--bed-filter", "tracks.io"),
    Opt(("--ref-chunksize"), "tracks.io"),
    Opt(("--aln-chunksize"), "tracks.io"),
    Opt(("-c", "--cov"), action="store_true", help="Output track coverage for each reference position"),
    Opt("--ref", "tracks"), 
    Opt(("-m", "--pore-model"), "pore_model", "name", default=None),
    Opt(("-p", "--processes"), "tracks.io"),
    Opt(("-o", "--outfile"), type=str),
)

ALIGN_OPTS =  DTW_OPTS + (
    MutexOpts("output", [
        Opt(("-o", "--bam-out"), "tracks.io", nargs="?", const="-"),
        Opt("--tsv-out", "tracks.io", nargs="?", const="-"),
        Opt("--eventalign-out", "tracks.io", nargs="?", const="-"),
    ]),
    Opt(("-m", "--pore-model"), "pore_model", "name", default=None),
    Opt("--bam-f5c", "tracks.io", action="store_true"),
    Opt("--tsv-cols", "tracks.io", type=comma_split, default="dtw"),
    Opt("--tsv-na", "tracks.io", nargs="?", const="-"),
    Opt("--tsv-noref", "tracks.io", action="store_true"),
    Opt("--eventalign-flags", "tracks.io", type=comma_split),
    Opt(("--norm-iterations"), "dtw"),
    Opt("--mask-skips", "tracks", nargs="?", const="keep_best"),
    Opt("--skip-cost", "dtw"),
    Opt("--stay-cost", "dtw"),
    Opt("--move-cost", "dtw"),
)

TRAIN_OPTS = (
    Opt(("-i", "--train-iterations"), "train", "iterations"), 
    Opt(("-m", "--init-model"), "train"),
    Opt("--init-mode", "train"),
    Opt("--moves-avg", "train", type=comma_split),
    Opt(("-k", "--kmer-len"), "train"),
    Opt("--kmer-samples", "train"), 
    Opt("--buffer-size", "train"), 
    Opt(("-d", "--max-moves-dist"), "train"), 
    Opt(("--train-mean"), "train", action="store_true"), 
    Opt("--out-dir", "tracks.io", "model_dir"),
    Opt(("-a", "--append"), "train", action="store_true"),
    Opt("--skip-dtw", "train", action="store_true"),
    Opt("--mask-skips", "tracks", nargs="?", default="keep_best"),
    Opt(("--norm-iterations"), "dtw", default=1),
    Opt("--skip-cost", "dtw", default=4),
    Opt("--stay-cost", "dtw"),
    Opt("--move-cost", "dtw"),
) + DTW_OPTS


READSTATS_OPTS = (
    Opt("--bam-in", "tracks.io", nargs="+"), 
    Opt("--layers", "readstats", nargs="+"),
    Opt("--stats", "readstats", nargs="+"),
    Opt(("-R", "--region"), "tracks", "ref_bounds", type=RefCoord),
    Opt(("-s", "--summary-stats"), "readstats", type=comma_split),
    CONFIG_OPT,
)

REFPLOT_OPTS = (
    Opt("--bam-in", "tracks.io", nargs="+"), #, required=True
    Opt(("-R", "--region"), "tracks", "ref_bounds", type=RefCoord),
    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt(("-l", "--read-filter"), "read_index", type=parse_read_ids),
    Opt(("-L", "--layer"), "refplot"),
    CONFIG_OPT,
)

DOTPLOT_OPTS = (
    Opt("--bam-in", "tracks.io", nargs="+"), #, required=True

    Opt(("-o", "--out-prefix"), type=str, default=None, help="If included will output images with specified prefix, otherwise will display interactive plot."),

    Opt("--ref", "tracks"), 
    Opt("--names", "tracks.io", "input_names", type=comma_split), 
    Opt("--reads", "read_index", "paths", nargs="+", type=str),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),

    Opt(("-f", "--out-format"), default="svg", help="Image output format. Only has an effect with -o option.", choices={"pdf", "svg", "png"}),
    Opt(("-R", "--region"), "tracks", "ref_bounds", type=RefCoord),
    Opt(("-l", "--read-filter"), "read_index", type=parse_read_ids),
    Opt(("-L", "--layers"), "dotplot", "layers", type=comma_split),
    Opt(("-p", "--pore-model"), "pore_model", "name", default=None),
    Opt(("--multi-background"), "sigplot", action="store_true"),
    Opt(("--no-model"), "sigplot", action="store_true"),
    Opt(("--svg"), "vis", action="store_true"),
    CONFIG_OPT,
)

TRACKPLOT_OPTS = (
    Opt(("-R", "--region"), "tracks", "ref_bounds", type=RefCoord, required=True),
    Opt("--bam-in", "tracks.io", nargs="+", required=True),

    Opt("--ref", "tracks"), 
    Opt("--read-paths", "read_index", "paths", nargs="+", type=str),
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),
    Opt("--pore-model", "pore_model", "name"),
    Opt(("--svg"), "vis", action="store_true"),

    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt(("-l", "--read-filter"), "read_index", type=parse_read_ids),
    Opt(("-H", "--panel-heights"), "trackplot", nargs="+", type=int),
    Opt(("--shared-refs-only"), "tracks", action="store_true"),
    Opt(("--shared-reads-only"), "tracks", action="store_true"),
    Opt(("--share-reads"), "trackplot", action="store_true"),
    Opt(("--hover-read"), "trackplot", action="store_true"),
    Opt(("-o", "--outfile"), "trackplot"),
    CONFIG_OPT,
)

BROWSER_OPTS = (
    Opt(("-R", "--region"), "tracks", "ref_bounds", type=RefCoord, required=True, help="Reference coordinates to visualize (chr:start-end)"),
    Opt("--bam-in", "tracks.io", nargs="+"), #, required=True
    Opt(("--shared-reads-only"), "tracks", action="store_true"),

    Opt("--reads", "read_index", "paths", nargs="+", type=str),
    Opt("--ref", "tracks"), 
    Opt(("-x", "--read-index"), "read_index"),
    Opt(("-r", "--recursive"), "read_index", action="store_true"),
    Opt("--rna", fn="set_r94_rna", help="Should be set for direct RNA data"),

    Opt(("-l", "--read_filter"), "tracks", type=parse_read_ids),
    Opt(("-f", "--full-overlap"), "tracks", action="store_true"),
    Opt("--pore-model", "pore_model", "name"),
    Opt("--names", "tracks.io", "input_names", type=comma_split), 
    Opt(("-p", "--port"), help="Browser port", default=8000),
    Opt(("-o", "--outfile"), "trackplot"),
    CONFIG_OPT,
)

def panel_opt(name):
    return (lambda arg: (name, arg))

TRACKPLOT_PANEL_OPTS = (
    Opt("--bases", dest="panels",
        metavar="LAYER", action="append_const", const=("bases",True), #type=panel_opt("bases"),
        help="Display a ref-by-read matrix of specified alignment layer"), 

    Opt("--mat", dest="panels",
        metavar="LAYER", action="append", type=panel_opt("mat"),
        help="Display a ref-by-read matrix of specified alignment layer"), 

    Opt("--box", dest="panels", #"trackplot", "panels", 
        metavar="LAYER", action="append", type=panel_opt("box"),
        help="Display a boxplot of specified layer"), 

    Opt("--line", dest="panels", #"trackplot", "panels", 
        metavar="LAYER.STAT", action="append", type=panel_opt("line"),
        help="Display a line plot of specifed layer summary statistic"), 

    Opt("--scatter", dest="panels", #"trackplot", "panels", 
        metavar="LAYER.STAT", action="append", type=panel_opt("scatter"),
        help="Display a line plot of specifed layer summary statistic"), 
)


CMDS = {
    "align" : ("align", 
        "Perform DTW alignment guided by basecalled alignments", ALIGN_OPTS), 
    "convert" : ("io", "Convert between signal alignment file formats", CONVERT_OPTS),
    "train" : ("train", 
        "Iteratively train a new k-mer pore model", TRAIN_OPTS), 
    "refstats" : ("stats.refstats", "Calculate per-reference-coordinate statistics", REFSTATS_OPTS),
    "readstats" : ("stats.readstats", "Compute per-read summary statistics", READSTATS_OPTS),
    "compare" : ("stats.layerstats", "Compute distance between alignments of the same reads", COMPARE_OPTS),
    "dotplot" : ("vis.dotplot", "Plot signal-to-reference alignment dotplots", DOTPLOT_OPTS),
    "refplot" : ("vis.refplot", "Plot alignment tracks and per-reference statistics", REFPLOT_OPTS),
    "trackplot" : ("vis.trackplot", "Plot alignment tracks and per-reference statistics", TRACKPLOT_OPTS+TRACKPLOT_PANEL_OPTS),
    "browser" : ("vis.browser", "Interactive signal alignment genome browser", BROWSER_OPTS),
}

_help_lines = [
    "Utility for Nanopore Current ALignment to Large Expanses of DNA", "",
    "subcommand options:",
    "Signal Alignment:",
    "\talign      Perform DTW alignment guided by basecalled alignments",
    "\ttrain      Train new k-mer pore models",
    "\tconvert    Convert between signal alignment file formats",
    "",
    "Analysis:",
    "\trefstats   Calculate per-reference-coordinate statistics",
    "\treadstats  Compute per-read summary statistics",
    "\tcompare    Compare multiple alignments of the same read",
    #"\tlayerstats Compute, compare, and query alignment layers", "",
    "",
    "Visualization:",
    "\tdotplot    Plot signal-to-reference alignment dotplots",
    "\ttrackplot  Plot alignment tracks and per-reference statistics",
    "\tbrowser    Interactive signal alignment genome browser",
]

HELP = "\n".join([
    textwrap.fill(line,
        width=75,
        drop_whitespace=False, 
        replace_whitespace=False, 
        tabsize=2, 
        subsequent_indent=" "*13)
    for line in _help_lines])


def main():
    parser = ArgParser(CMDS, HELP)

    module, cmd, conf = parser.parse_args()


    if module is not None:
        m = importlib.import_module(f".{module}", "uncalled4")
        fn = getattr(m, cmd)
        if conf.cprof is not None:
            cProfile.runctx(f"fn(conf=conf)",
             {"fn" : fn, "conf": conf},{},
             conf.cprof)
        else:
            ret = fn(conf=conf)
            if isinstance(ret, pd.DataFrame):
                ret.to_csv(sys.stdout, sep="\t")

    else:
        parser.print_help()

if __name__ == "__main__":
    main()
