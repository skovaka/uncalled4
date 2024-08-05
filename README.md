# Uncalled4

**U**tility for **N**anopore **C**urrent **Al**ignment to **L**arge **E**xpanses of **D**NA (version 4)

A toolkit for nanopore signal alignment, analysis, and visualization

<img src="logo.png" width="250">

Performs accurate basecaller-guided nanopore signal alignment, similar to [nanopolish eventalign](https://github.com/jts/nanopolish) or [tombo resquiggle](https://github.com/nanoporetech/tombo), to map nanopore signal segments to reference nucleotides. Also supports conversion of any signal alignments to an efficient BAM format, allowing for interactive visualizations, modification detection, and other signal analyses.

For [real-time targeted sequencing](https://www.nature.com/articles/s41587-020-0731) via rapid signal mapping, see [UNCALLED](https://github.com/skovaka/UNCALLED)

### Read our pre-print on BioRxiv!

[**Uncalled4 improves nanopore DNA and RNA modification detection via fast and accurate signal alignment** \
S. Kovaka, P.W. Hook, K.M. Jenike, V. Shivakumar, L.B. Morina, R. Razaghi, W. Timp, M.C. Schatz \
BioRxiv (2024)
](https://www.biorxiv.org/content/10.1101/2024.03.05.583511v1)

## Table of Contents

- [Installation](#installation)
- [Overview](#overview)
  - [`align`: Perform DTW alignment guided by basecalled alignments](#align)
  - [`train`: Train new k-mer pore models](#train)
  - [`convert`: Import DTW alignments produced by Nanopolish or Tombo](#convert)
- [Visualization](#visualization)
  - [`dotplot`: Plot signal-to-reference alignment dotplots](#dotplot)
  - [`trackplot`: Plot alignment tracks and per-reference statistics](#trackplot)
  - [`browser`: Interactive signal alignment genome browser](#browser)
- [Analysis](#analysis)
  - [`refstats`: Calculate per-reference-coordinate statistics](#refstats)
  - [`readstats`: Calculate per-read summary statistics](#readstats)
  - [`compare`: Compare signal alignments](#compare)
- [Pore Models](#pore-models)
- [Alignment Layers](#alignment-layers)
- [Future Work](#future-work)
- [Release Notes](#release-notes)

## Installation

```
pip install uncalled4
```

OR

```
git clone https://github.com/skovaka/uncalled4.git
cd uncalled4
pip install .
```

Requires python >= 3.8 and GCC>=4.9.2

Uncalled4 is currently only tested Linux, primarily on Ubuntu 18.04 and 20.04

## Overview

The `example/` directory contains scripts to download example data and test each command:
```
cd uncalled4/example
./download.sh
./run.sh
```

Uncalled4 signal alignment is guided by BAM alignments output by the basecaller with "move" tags (Dorado `--emit-moves --emit-sam` or Guppy `--moves_out`), which include `mv:`, `ts:`, and related BAM tags. This BAM file should contain reference-aligned reads (Dorado `--reference` or guppy `--align_ref`), unmapped will be ignored. To (re)align reads while preserving these tags, you can use `samtools fastq -T "mv,ts,pi,sp,ns" <file>.bam > file.fastq` (v1.16 or newer) to convert the tags to FASTQ format, then align using `minimap2 -y ...` to propagate the tags to the SAM file. **Guppy only includes BAM tag on primary alignments by default**. See [scripts/bamprep.py](scripts/) to copy tags to supplemental and secondary aligments.

Uncalled4 primarily stores signal alignments in BAM alignment tags, including per-reference signal coordinates and current summary statistics required for most signal analyses. A sorted and indexed (via `samtools index`) Uncalled4 BAM file is required for most visualization and analysis commands. [Nanopolish](https://github.com/jts/nanopolish), [f5c](https://github.com/hasindu2008/f5c), or [Tombo](https://github.com/nanoporetech/tombo) alignments can be converted to Uncalled4 format via the `convert` command.

Signal alignment requires a pore model to map k-mers to expected current. Uncalled4 will attempt to automatically detect the appropriate pore model from the input data, but may require you to specify a preset pore model or custom pore model. This can be specified using the `--pore-model` flag or `--flowcell` and `--kit` flags. `uncalled4 train` trains new pore models, either starting from a initialization pore model, or from scratch using basecaller moves.

Preliminary support for RNA004 alignment is implemented on the `dev` branch. Signal alignment quality appears to be better than RNA002, but the exact implementation may change before migrating to the main branch.

### `align`

Perform DTW alignment guided by basecalled alignments

```
uncalled4 align [-h] [--ref REF | --self] [--reads READS [READS ...]] --bam-in [BAM_IN] [-p PROCESSES] [--flowcell FLOWCELL]
                [--kit KIT] [--basecaller-profile BASECALLER_PROFILE] [--rna] [--ordered-out] [-f] [--kmer-shift KMER_SHIFT]
                [--bam-chunksize BAM_CHUNKSIZE] [--out-name OUT_NAME] [-r] [-l READ_FILTER] [-x READ_INDEX] [-n MAX_READS]
                [--count-events] [--del-max DEL_MAX] [--ins-max INS_MAX] [--unmask-splice] [--method METHOD] 
                [-c {abs_diff,z_score,norm_pdf}] [-b BAND_WIDTH] [-s BAND_SHIFT] [--mvcmp-mask MVCMP_MASK]
                [--max-norm-dist MAX_NORM_DIST] [--max-sd MAX_SD] [--min-aln-length MIN_ALN_LENGTH] [-N {ref_mom,model_mom}]
                [--zero-ts] [-C CONFIG] [-o [BAM_OUT] | --tsv-out [TSV_OUT] | --eventalign-out [EVENTALIGN_OUT]] [-m PORE_MODEL]
                [--bam-f5c] [--tsv-cols TSV_COLS] [--tsv-na [TSV_NA]] [--tsv-noref] [--eventalign-flags EVENTALIGN_FLAGS]       
                [--norm-iterations NORM_ITERATIONS] [--mask-skips [MASK_SKIPS]] [--skip-cost SKIP_COST] [--stay-cost STAY_COST]
                [--move-cost MOVE_COST]                                                                                         

required arguments:
  --ref [fasta] | --self     Reference to align to, or perform signal-to-read self-alignment
  --reads [path]             Paths to fast5, slow5, or pod5 files, or to directories containing those files (optionally recursive)
  --bam-in [bam_in]          BAM input file (or "-"/no argument for stdin) (default: None)

output arguments (one required):
  -o, --bam-out     BAM output file
  --tsv-out         TSV output file (or "-"/no argument for stdout) 
  --eventalign-out  Eventalign (nanopolish) output file (or "-"/no argument for stdout)

selected optional arguments:
  -p, --processess      Number of parallel processes (default: 1)
  --flowcell FLOWCELL   Flowcell used for sequencing (e.g. FLO-MIN106)
  --kit KIT             Kit used for sequencing (e.g. SQK-LSK109)
  -m, --pore-model       Custom pore model
  --rna                 Required for custom RNA pore models
  --tsv-cols TSV_COLS   TSV file output alignment layers (comma-separated, can also include "read_id" (default: dtw)
  --tsv-na [TSV_NA]     Missing value representation for TSV output (default: *)
  --tsv-noref           Will NOT output reference coordinates to TSV if True (default: False)
  --eventalign-flags    Eventalign optional flags (comma-separated list of ""print-read-names", "signal-index", "samples")
  -h, --help            Print full command line usage

```

`--bam-in` must be a BAM file produced by Dorado using `--emit-moves --emit-sam` flags or the Guppy `--moves_out` flag.

The `--self` option is an alternative to `--ref`, and performs signal-to-read alignment to the basecalled read, which must be defined in the `SEQ` field in the input BAM file.


### convert

Convert between signal alignment file formats

```
uncalled4 convert [-h] [--bam-in BAM_IN | --eventalign-in [EVENTALIGN_IN] | --tombo-in TOMBO_IN]   
                  [--eventalign-out [EVENTALIGN_OUT] | --tsv-out [TSV_OUT] | 
                  --m6anet-out [M6ANET_OUT]] [--tsv-cols TSV_COLS] 
                  [--eventalign-flags EVENTALIGN_FLAGS]
                  [--mask-skips [MASK_SKIPS]] [--ref FASTA] [--reads READS [READS ...]]          
                  [-l READ_FILTER] [-x READ_INDEX] [-r] [--rna] [-R REF_BOUNDS] [-f] [-a]
```

Generally only one `--*-in` and one `--*-out` option should be specified, with the exception of `--bam-out` where a template bam file should be specified via `--bam-in`. 

[nanopolish](https://github.com/jts/nanopolish) or [f5c](https://github.com/hasindu2008/f5c) `eventalign` should be run with the `--signal-index` and `--scale-events` options, and can be converted with `uncalled4 convert --eventalign-in <eventalign.txt> --bam-in <mm2.bam>`, where `<mm2.bam>` is the exact BAM file used to guide the eventalign command.

`--m6anet-out` efficently implements [m6anet](https://m6anet.readthedocs.io/en/latest/) `dataprep` for sorted Uncalled4 BAM files. This should be used with an [m6Anet model trained on Uncalled4 alignments](https://github.com/skovaka/uncalled4_supplemental_data/tree/main/m6anet_model) 

### `train`

Iteratively train a new k-mer pore model

Accepts most of the same paramters as `uncalled4 align`, in addition to number of iterations. First iteration must use some starting pore model, while subsequent iterations use  the pore model from the previous iteration.

```
uncalled4 train [-h] [-i TRAIN_ITERATIONS] [-m INIT_MODEL] [--init-mode INIT_MODE] [--moves-avg MOVES_AVG] [-k KMER_LEN]
                [--kmer-samples KMER_SAMPLES] [--buffer-size BUFFER_SIZE] [-d MAX_MOVES_DIST] [--train-mean] [--out-dir OUT_DIR] [-a]
                [--skip-dtw] [--mask-skips [MASK_SKIPS]] [--norm-iterations NORM_ITERATIONS] [--skip-cost SKIP_COST]
                [--stay-cost STAY_COST] [--move-cost MOVE_COST] [--ref REF | --self] [--reads READS [READS ...]] --bam-in [BAM_IN]
                [-p PROCESSES] [--flowcell FLOWCELL] [--kit KIT] [--basecaller-profile BASECALLER_PROFILE] [--rna] [--ordered-out]
                [-f] [--kmer-shift KMER_SHIFT] [--bam-chunksize BAM_CHUNKSIZE] [--out-name OUT_NAME] [-r] [-l READ_FILTER]
                [-x READ_INDEX] [-n MAX_READS] [--count-events] [--del-max DEL_MAX] [--ins-max INS_MAX] [--unmask-splice]
                [--method METHOD] [-c {abs_diff,z_score,norm_pdf}] [-b BAND_WIDTH] [-s BAND_SHIFT] [--mvcmp-mask MVCMP_MASK]
                [--max-norm-dist MAX_NORM_DIST] [--max-sd MAX_SD] [--min-aln-length MIN_ALN_LENGTH] [-N {ref_mom,model_mom}]
                [--zero-ts] [-C CONFIG]                                                                                              
```
                                                                                          
## Visualization

All visualizations are generated using Plotly.

### `dotplot`

Plot signal-to-reference alignment dotplots

```
uncalled4 dotplot [-h] [--bam-in BAM_IN [BAM_IN ...]] [-o OUT_PREFIX] [--ref REF] [--names NAMES] [--reads READS [READS ...]]
                  [-x READ_INDEX] [-r] [--rna] [-f {svg,png,pdf}] [-R REGION] [-l READ_FILTER] [-L LAYERS] [-p PORE_MODEL]
                  [--multi-background] [--no-model] [--svg] [-C CONFIG]

options:
  -h, --help            show this help message and exit
  --bam-in BAM_IN [BAM_IN ...]
                        BAM input file (default: None)
  -o OUT_PREFIX, --out-prefix OUT_PREFIX
                        If included will output images with specified prefix, otherwise will display interactive plot. (default: None)
  --ref REF             Reference FASTA file, must match --bam-in reference (default: None)
  --names NAMES         Names of tracks to read from input(s) (default: None)
  --reads READS [READS ...]
                        Paths to FAST5, SLOW5, or POD5 files, or to directories containing those files (optionally recursive) (default:
                        None)
  -x READ_INDEX, --read-index READ_INDEX
                        File containing a mapping of read IDs to filenames (default: None)
  -r, --recursive       Recursively search 'paths' for FAST5, SLOW5, or POD5 files (default: False)
  --rna                 Should be set for direct RNA data (default: None)
  -f {svg,png,pdf}, --out-format {svg,png,pdf}
                        Image output format. Only has an effect with -o option. (default: svg)
  -R REGION, --region REGION
                        Only load reads which overlap these coordinates (default: None)
  -l READ_FILTER, --read-filter READ_FILTER
                        List of read IDs to load, or file containing one read ID per line (default: None)
  -L LAYERS, --layers LAYERS
  -p PORE_MODEL, --pore-model PORE_MODEL
                        Model preset name or TSV filename (default: None)
  --multi-background    Will plot multiple stacked background colors for multiple tracks if True (default: False)
  --no-model            Will not plot the expected reference signal if True (default: False)
  --svg                 Make SVG-friendly figures (default: False)
  -C CONFIG, --config CONFIG
                        Configuration file in TOML format (default: None)

```

### `trackplot`

Plot alignment tracks and per-reference statistics

Trackplots are defined by a series of panels displaying different layers. A `mat` panel display a heatmap of layer values for each ref/read coordinate on each track. A `box` panel displays boxplots of layer summary statistics for each track. `line` and `scatter` panels display [`refstats`](#refstats) summary statistics, specified by `<layer>.<statistic>` (e.g. `current.mean`, `model_diff.median`).

```
uncalled4 trackplot [-h] -R REGION --bam-in BAM_IN [BAM_IN ...] [--ref REF] [--read-paths READ_PATHS [READ_PATHS ...]]
                    [-x READ_INDEX] [-r] [--rna] [--pore-model PORE_MODEL] [--svg] [-f] [-l READ_FILTER]
                    [-H PANEL_HEIGHTS [PANEL_HEIGHTS ...]] [--shared-refs-only] [--shared-reads-only] [--share-reads] [--hover-read]
                    [-o OUTFILE] [-C CONFIG] [--bases] [--mat LAYER] [--box LAYER] [--line LAYER.STAT] [--scatter LAYER.STAT]

options:
  -h, --help            show this help message and exit
  -R REGION, --region REGION
                        Only load reads which overlap these coordinates (default: None)
  --bam-in BAM_IN [BAM_IN ...]
                        BAM input file (default: None)
  --ref REF             Reference FASTA file, must match --bam-in reference (default: None)
  --read-paths READ_PATHS [READ_PATHS ...]
                        Paths to FAST5, SLOW5, or POD5 files, or to directories containing those files (optionally recursive) (default:
                        None)
  -x READ_INDEX, --read-index READ_INDEX
                        File containing a mapping of read IDs to filenames (default: None)
  -r, --recursive       Recursively search 'paths' for FAST5, SLOW5, or POD5 files (default: False)
  --rna                 Should be set for direct RNA data (default: None)
  --pore-model PORE_MODEL
                        Model preset name or TSV filename (default: )
  --svg                 Make SVG-friendly figures (default: False)
  -f, --full-overlap    If true will only include reads which fully cover reference bounds (default: False)
  -l READ_FILTER, --read-filter READ_FILTER
                        List of read IDs to load, or file containing one read ID per line (default: None)
  -H PANEL_HEIGHTS [PANEL_HEIGHTS ...], --panel-heights PANEL_HEIGHTS [PANEL_HEIGHTS ...]
                        Relative height of each panel (default: None)
  --shared-refs-only    If true will only contain reference positions where all tracks have sufficient coverage (see min_coverage) (default:
                        False)
  --shared-reads-only   If true will only contain reads shared between all tracks (default: False)
  --share-reads         If True will only display reads shared by all alignment tracks with shared y-axis (default: False)
  --hover-read          If True will display read_id in mat hover (default: False)
  -o OUTFILE, --outfile OUTFILE
                        Output file (default: None)
  -C CONFIG, --config CONFIG
                        Configuration file in TOML format (default: None)
  --bases               Display a ref-by-read matrix of specified alignment layer (default: None)
  --mat LAYER           Display a ref-by-read matrix of specified alignment layer (default: None)
  --box LAYER           Display a boxplot of specified layer (default: None)
  --line LAYER.STAT     Display a line plot of specifed layer summary statistic (default: None)
  --scatter LAYER.STAT  Display a line plot of specifed layer summary statistic (default: None)

```

### `browser`

Interactive signal alignment genome browser

Integrates trackplot and dotplot for interactive alignment browsing

```
uncalled4 browser [-h] -R REGION [--bam-in BAM_IN [BAM_IN ...]] [--shared-reads-only] [--reads READS [READS ...]] [--ref REF]
                  [-x READ_INDEX] [-r] [--rna] [-l READ_FILTER] [-f] [--pore-model PORE_MODEL] [--names NAMES] [-p PORT] [-o OUTFILE]
                  [-C CONFIG]

options:
  -h, --help            show this help message and exit
  -R REGION, --region REGION
                        Reference coordinates to visualize (chr:start-end) (default: None)
  --bam-in BAM_IN [BAM_IN ...]
                        BAM input file (default: None)
  --shared-reads-only   If true will only contain reads shared between all tracks (default: False)
  --reads READS [READS ...]
                        Paths to FAST5, SLOW5, or POD5 files, or to directories containing those files (optionally recursive) (default:
                        None)
  --ref REF             Reference FASTA file, must match --bam-in reference (default: None)
  -x READ_INDEX, --read-index READ_INDEX
                        File containing a mapping of read IDs to filenames (default: None)
  -r, --recursive       Recursively search 'paths' for FAST5, SLOW5, or POD5 files (default: False)
  --rna                 Should be set for direct RNA data (default: None)
  -l READ_FILTER, --read_filter READ_FILTER
                        Only load reads which overlap these coordinates (default: None)
  -f, --full-overlap    If true will only include reads which fully cover reference bounds (default: False)
  --pore-model PORE_MODEL
                        Model preset name or TSV filename (default: )
  --names NAMES         Names of tracks to read from input(s) (default: None)
  -p PORT, --port PORT  Browser port (default: 8000)
  -o OUTFILE, --outfile OUTFILE
                        Output file (default: None)
  -C CONFIG, --config CONFIG
                        Configuration file in TOML format (default: None)
```

## Analysis

These functions compute statistics over reference and read coordinates. `refstats` computes summary and comparison statistics (e.g. Kolmogorov-Smirnov test) over reference coordinates, while `layerstats` maintains the read-by-reference dimensions of the DTW alignments.

(Coming soon: `readstats` to compute read-level statistics)

### `refstats`

Calculate per-reference-coordinate statistics

```
uncalled4 refstats [-h] [--bam-in BAM_IN [BAM_IN ...]] [--layers LAYERS [LAYERS ...]] [--stats REFSTATS [REFSTATS ...]] 
                   [-t TRACKS] [-R REGION] [--min-coverage MIN_COVERAGE] [--bed-filter BED_FILTER] [--ref-chunksize REF_CHUNKSIZE] 
                   [--aln-chunksize ALN_CHUNKSIZE] [-c] [--ref REF] [-m PORE_MODEL] [-p PROCESSES] [-o OUTFILE]

required arguments:
  --bam-in BAM_IN [BAM_IN ...]
                        BAM input file (default: None)
  --layers LAYERS [LAYERS ...]
                        Comma-separated list of layers over which to compute summary statistics (default: [])
  --stats REFSTATS [REFSTATS ...]
                        Comma-separated list of summary statistics to compute. Some statisitcs (ks) can only be used if exactly two tracks
                        are provided {q5,kurt,q75,KS,q95,max,var,skew,median,stdv,mean,q25,min} (default: None)

options:
  -h, --help            show this help message and exit
  -t TRACKS, --tracks TRACKS
                        Names of tracks to read from input(s) (default: None)
  -R REGION, --region REGION
                        Only load reads which overlap these coordinates (default: None)
  --min-coverage MIN_COVERAGE
                        Reference positions with less than this coverage will be excluded from each track (or all tracks if shared_refs_only
                        is true) (default: 1)
  --bed-filter BED_FILTER
                        Only parse regions in BED file (default: None)
  --ref-chunksize REF_CHUNKSIZE
                        Number of reference coordinates to query for iteration (default: 10000)
  --aln-chunksize ALN_CHUNKSIZE
                        Number of alignments to query for iteration (default: 500)
  -c, --cov             Output track coverage for each reference position (default: False)
  --ref REF             Reference FASTA file, must match --bam-in reference (default: None)
  -m PORE_MODEL, --pore-model PORE_MODEL
                        Model preset name or TSV filename (default: None)
  -p PROCESSES, --processes PROCESSES
                        Number of parallel processes (default: 1)
  -o OUTFILE, --outfile OUTFILE
```

### `readstats`

Compute per-read summary statistics

```
uncalled4 readstats [-h] [--bam-in BAM_IN [BAM_IN ...]] [--layers LAYERS [LAYERS ...]] [--stats STATS [STATS ...]] [-R REGION]
                    [-s SUMMARY_STATS] [-C CONFIG]

options:
  -h, --help            show this help message and exit
  --bam-in BAM_IN [BAM_IN ...]
                        BAM input file (default: None)
  --layers LAYERS [LAYERS ...]
                        Which layers to compute statistics (default: None)
  --stats STATS [STATS ...]
                        Summary statistics to compute (any builtin numpy function, e.g. mean, std, etc) (default: None)
  -R REGION, --region REGION
                        Only load reads which overlap these coordinates (default: None)
  -s SUMMARY_STATS, --summary-stats SUMMARY_STATS
                        Summary statistics to compute for "model_diff" command. (default: ['mean'])
  -C CONFIG, --config CONFIG
                        Configuration file in TOML format (default: None)
```

### `compare`

Compute distance between alignments of the same reads

```
uncalled4 compare [-h] [--bam-in BAM_IN [BAM_IN ...]] [-t TRACKS] [-l READ_FILTER] [-R REGION] [-m] [--tsv-cols TSV_COLS]
                  [--tsv-na [TSV_NA]] [--tsv-noref] [--tsv-out [TSV_OUT]] [-C CONFIG]

options:
  -h, --help            show this help message and exit
  --bam-in BAM_IN [BAM_IN ...]
                        BAM input file (default: None)
  -t TRACKS, --tracks TRACKS
                        Names of tracks to read from input(s) (default: None)
  -l READ_FILTER, --read-filter READ_FILTER
                        Only load reads which overlap these coordinates (default: None)
  -R REGION, --region REGION
                        Only load reads which overlap these coordinates (default: None)
  -m, --moves           Compare against basecalled alignment. If two tracks input will look for "moves" group in second track, otherwise
                        will look in the first track. (default: False)
  --tsv-cols TSV_COLS   TSV file output alignment layers (comma-separated, can also include "read_id" (default: dtwcmp,mvcmp)
  --tsv-na [TSV_NA]     Missing value representation for TSV output (default: *)
  --tsv-noref           Will NOT output reference coordinates to TSV if True (default: False)
  --tsv-out [TSV_OUT], -o [TSV_OUT]
                        TSV output file (or "-"/no argument for stdout) (default: -)
  -C CONFIG, --config CONFIG
                        Configuration file in TOML format (default: None)

```

## Pore Models

Uncalled4 pore models map k-mers to their expected current characteristics for a specific sequencing chemistry, minimally defining the expected mean current (`current.mean`) for each k-mer. Pore models trained by Uncalled4 include means and standard deviations for per-kmer current means, current standard deviations, and dwell time measured in raw sample length: `current.mean`/`current.stdv`, `current_sd.mean`/`current_sd.stdv`, `length.mean`/`length.stdv`. In addition to per-kmer statistics, each model has a defined k-mer `shift` used to define the central base that has the most effect on the current, the picoamp mean and standard deviation (`pa_mean` and `pa_stdv`) that can be used to scale the normalized current values to picoamps, and other model-specifc parameters like `sample_rate` and `bases_per_sec`. A `reverse` parameter is also included, which is set to `True` for RNA to indicated reversed sequencing direction.

Four pore model presets are provided by Uncalled4: `dna_r10.4.1_400bps_9mer`, `dna_r9.4.1_400bps_6mer`, `rna_r9.4.1_70bps_5mer`, and `rna004_130bps_9mer`. These are stored efficently in binary NumPy "npz" files, and can be converted to TSV format using the provided [model2tsv.py](scripts/model2tsv.py) script. 

Custom pore models can be provided in TSV format via the `--pore-model` command line argument, which should minimally include columns named `kmer` and `current.mean`. Column names for current levels are also aliased to support Nanopolish and other similar models, so `current.mean` can be `level_mean` or `mean`, `current_sd.mean` can be `sd_mean` or `stdv`, etc. K-mer offsets can also be defined for custom pore models using the `--kmer-shift` option. All pore models are automatically normalized such that the mean and standard deviation of `current.mean` is 0 and 1 respecively, which is required for BAM encoding, and the resulting BAM file will store scaling factors to convert to the original input values in `pa_mean` and `pa_stdv`.

Uncalled4 will attempt to automatically detect the sequencing chemistry based on POD5/FAST5/BLOW5 metadata, which is required even with custom pore models to determine the appropriate offset to use for basecaller moves metadata. If this cannot be automatically detected, the `--basecaller-profile` must also be provided, which is defined similar to the pore model presets but without a defined k-mer length: either `dna_r10.4.1_400bps`, `dna_r10.4.1_260bps`, `dna_r9.4.1_400bps`, `rna_r9.4.1_70bps`, or `rna004_130bps`. If you are using a sequenicng chemistry which does not have an appropirate preset, please submit an issue and we can implement one.

## Alignment Layers

Uncalled4 stores signal alignments as a set of **layers** associated with read and reference coordinates. Each layer is derived from the read signal (e.g. `current`), the reference sequence (e.g. `kmer`), or other layers (e.g. `model_diff`). Layers are organized into **layer groups**: `dtw` for signal alignments, `bcaln` for projected basecalled alignments, and `cmp` for alignment comparisons. Several subcommands detailed above take layer names as input, which should generaly be in the form `<group>.<layer>`. Below is a table of currently available layers:

| Group | Layer   | Description |
|-------|---------|-------------|
| dtw   | current | Normalized mean read signal current |
| dtw   | current_sd | Normalized read signal current standard deviation |
| dtw   | start | Read signal sample start index |
| dtw   | length | Read signal sample length |
| dtw   | dwell   | Signal dwell time (ms/nt, proportional to length) |
| dtw   | model_diff | Difference between predicted (via a pore model) and observed current|
| dtw   | abs_diff | Absolute value of dtw.model_diff |
| dtw   | kmer | Binarized reference k-mer |
| dtw   | events | Number of raw signal events aligned ("stays" > 1, "skips" < 1) |
| dtw   | base | Binarized reference base |
| moves   | start | Estimated read signal sample start index |
| moves   | length | Estimated read signal sample length |
| moves   | indel | Number of inserted (>0) or deleted (<0) nucleotides |
| mvcmp   | mean_ref_dist | Mean reference distance between signal alignment and ref-moves |
| mvcmp   | jaccard | Raw sample jaccard distances between signal alignment and ref-moves |
| dtwcmp   | mean_ref_dist | Mean reference distance between two alignments (must first run [`layerstats compare`](#compare)) |
| dtwcmp   | jaccard | Raw sample jaccard distances between two alignments (must first run [`layerstats compare`](#compare)) |

All tracks must be written to the same database for multi-track visualization and analysis (e.g. comparing alignments, calculating KS statistics). You can merge multiple databases into a single file using [`uncalled db merge`](#db)

## Release Notes
- v4.1.0:  Major update.
  - Added RNA004 support
  - Added signal-to-read alignment via `align --self`
  - Changed r10.4.1 output coordinates to be centered on central base
  - Changed all positional arguments to flags
- v4.0.0:  Pre-print release
For earlier development history, see https://github.com/skovaka/UNCALLED/tree/dev4
