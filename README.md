# Uncalled4

**U**tility for **N**anopore **C**urrent **Al**ignment to **L**arge **E**xpanses of **D**NA (version 4)

A toolkit for nanopore signal alignment, analysis, and visualization

<img src="logo.png" width="250">

Performs accurate basecaller-guided nanopore signal alignment, similar to [nanopolish eventalign](https://github.com/jts/nanopolish) or [tombo resquiggle](https://github.com/nanoporetech/tombo), to map nanopore signal segments to reference nucleotides. Also supports conversion of any signal alignments to an efficient BAM format, allowing for interactive visualizations, modification detection, and other signal analyses.

For [real-time targeted sequencing](https://www.nature.com/articles/s41587-020-0731) via rapid signal mapping, see [UNCALLED](https://github.com/skovaka/UNCALLED)

## Table of Contents

- [Installation](#installation)
- [Overview](#overview)
  - [`align`: Perform DTW alignment guided by basecalled alignments](#align)
  - [`train`: Train new k-mer pore models](#train)
  - [`convert`: Import DTW alignments produced by Nanopolish or Tombo](#convert)
- [Visualization](#visualization)
  - [`trackplot`: Plot alignment tracks and per-reference statistics](#trackplot)
  - [`browser`: Interactive signal alignment genome browser](#browser)
- [Analysis](#analysis)
  - [`refstats`: Calculate per-reference-coordinate statistics](#refstats)
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

Requires python >= 3.8

Uncalled4 is currently only tested Linux

## Overview

The `example/` directory contains scripts to download example data and test each command:
```
cd uncalled4/example
./download.sh
./run.sh
```

Uncalled4 signal alignment is guided by BAM alignments output by the basecaller (Guppy or Dorado) using the `--moves_out` option, which includes the `mv:` and `ts:` BAM tags. These tags encode an approximate alignment between the basecalled sequence and signal, which Uncalled4 uses to improve the speed and accuracy of signal alignment. To re-align reads while preserving these tags, you can use `samtools fasta -T "mv,ts" <file>.bam > file.fastq` (v1.16 or newer) to convert the tags to FASTQ format, then align using `minimap2 -y ...` to propagate the tags to the SAM file.

Uncalled4 primarily stores signal alignments in BAM alignment tags, including per-reference signal coordinates and current summary statistics required for most signal analyses. A sorted and indexed (via `samtools index`) Uncalled4 BAM file is required for most visualization and analysis commands. [Nanopolish](https://github.com/jts/nanopolish), [f5c](https://github.com/hasindu2008/f5c), or [Tombo](https://github.com/nanoporetech/tombo) alignments can be converted to Uncalled4 format via the `convert` command.

Signal alignment requires a pore model to map k-mers to expected current. Uncalled4 will attempt to automatically detect the appropriate pore model from the input data, but may require you to specify a preset pore model or custom pore model. This can be specified using the `--pore-model` flag or `--flowcell` and `--kit` flags. `uncalled4 train` trains new pore models, either starting from a initialization pore model, or from scratch using basecaller moves.

### `align`

Perform DTW alignment guided by basecalled alignments

```
usage: uncalled4 align [-h] ref_index paths [paths ...] --bam-in [BAM_IN] [-p PROCESSES] [--flowcell FLOWCELL] 
					   [--kit KIT] [--rna] [--ordered-out] [-f] [--kmer-shift KMER_SHIFT] [--bam-chunksize BAM_CHUNKSIZE] 
					   [--out-name OUT_NAME] [-r] [-l READ_FILTER] [-x READ_INDEX] [-n MAX_READS] [--count-events] 
                       [--del-max DEL_MAX] [--ins-max INS_MAX] [--mask-skips [MASK_SKIPS]] [--unmask-splice] 
                       [-c {abs_diff,z_score,norm_pdf}] [--skip-cost SKIP_COST] [--stay-cost STAY_COST] [--move-cost MOVE_COST] 
                       [-b BAND_WIDTH] [-s BAND_SHIFT] [--mvcmp-mask MVCMP_MASK] [--max-norm-dist MAX_NORM_DIST] [--max-sd MAX_SD]
                       [--min-aln-length MIN_ALN_LENGTH] [-N {ref_mom,model_mom}] [-C CONFIG] [-o [BAM_OUT] | --tsv-out
                       [TSV_OUT] | --eventalign-out [EVENTALIGN_OUT]] [-m PORE_MODEL] [--tsv-cols TSV_COLS] [--tsv-na [TSV_NA]]
                       [--tsv-noref] [--eventalign-flags EVENTALIGN_FLAGS]

required arguments:
  ref_index             FASTA file indexed by faidx
  paths                 Paths to fast5, slow5, or pod5 files, or to directories containing those files (optionally recursive)
  --bam-in [BAM_IN]     BAM input file (or "-"/no argument for stdin) (default: None)

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

`--bam-in` must be a BAM file produced by Guppy or Dorado using the `--moves_out` option.


### convert

Convert between signal alignment file formats

```
usage: uncalled4 convert [-h] [--bam-in BAM_IN | --eventalign-in [EVENTALIGN_IN]
                         | --tombo-in TOMBO_IN]                                                   
                         [--eventalign-out [EVENTALIGN_OUT] | --tsv-out [TSV_OUT] | 
                         --m6anet-out [M6ANET_OUT]] [--tsv-cols TSV_COLS] 
                         [--eventalign-flags EVENTALIGN_FLAGS]
                         [--mask-skips [MASK_SKIPS]] [--reads READS [READS ...]]          
                         [-l READ_FILTER] [-x READ_INDEX] [-r] [--rna] [-R REF_BOUNDS] [-f] [-a]
                         ref_index       
```

Generally only one `--*-in` and one `--*-out` option should be specified, with the exception of `--bam-out` where a template bam file should be specified via `--bam-in`. 

[nanopolish](https://github.com/jts/nanopolish) or [f5c](https://github.com/hasindu2008/f5c) `eventalign` should be run with the `--signal-index` and `--scale-events` options, and can be converted with `uncalled4 convert --eventalign-in <eventalign.txt> --bam-in <mm2.bam>`, where `<mm2.bam>` is the exact BAM file used to guide the eventalign command.

`--m6anet-out` efficently implements [m6anet](https://m6anet.readthedocs.io/en/latest/) `dataprep` for sorted Uncalled4 BAM files. This should be used with an m6Anet model trained on Uncalled4 alignments (link coming soon) 

### `train`

Iteratively train a new k-mer pore model

Accepts most of the same paramters as `uncalled4 align`, in addition to number of iterations. First iteration must use some starting pore model, while subsequent iterations use  the pore model from the previous iteration.

```
usage: uncalled4 train [-h] [-i ITERATIONS] [--kmer-samples KMER_SAMPLES] [--buffer-size BUFFER_SIZE]          
                       [-d MAX_BCALN_DIST] [--use-median] [--out-dir OUT_DIR] [-a] [--skip-dtw] [-p PROCESSES] 
                       [--bam-chunksize BAM_CHUNKSIZE] [--guppy-in GUPPY_IN] --bam-in [BAM_IN]                 
                       [--out-name OUT_NAME] [-r] [-l READ_FILTER] [-x READ_INDEX] [-n MAX_READS]             
                       [--del-max DEL_MAX] [--mask-skips [MASK_SKIPS]] [--mask-indels MASK_INDELS] [-f]        
                       [-m PORE_MODEL] [--kmer-shift KMER_SHIFT] [--save-bands] [--full-overlap] [--rna]       
                       [-R REF_BOUNDS] [-c {abs_diff,z_score,norm_pdf}] [--skip-cost SKIP_COST]                
                       [--stay-cost STAY_COST] [--move-cost MOVE_COST] [-b BAND_WIDTH] [-s BAND_SHIFT]         
                       [-N {ref_mom,model_mom}] [--norm-median] [--norm-seg] [--bc-group BC_GROUP] [-C CONFIG] 
                       ref_index read_files [read_files ...]                                              
```
                                                                                          
## Visualization

All visualizations are generated using Plotly.

### `trackplot`

Plot alignment tracks and per-reference statistics

Trackplots are defined by a series of panels displaying different layers. A `mat` panel display a heatmap of layer values for each ref/read coordinate on each track. A `box` panel displays boxplots of layer summary statistics for each track. `line` and `scatter` panels display [`refstats`](#refstats) summary statistics, specified by `<layer>.<statistic>` (e.g. `current.mean`, `model_diff.median`).

```
usage: uncalled4 trackplot [-h] ref_bounds bam_in [bam_in ...] [--ref REF]         
                          [--reads READS [READS ...]] [-x READ_INDEX] [-r] [--rna]
                          [--pore-model PORE_MODEL] [-f] [-l READ_FILTER]                  
                          [-H PANEL_HEIGHTS [PANEL_HEIGHTS ...]] [--shared-refs-only]
                          [--shared-reads-only] [--share-reads] [--hover-read] [-o OUTFILE]     
                          [--mat LAYER] [--box LAYER] [--line LAYER.STAT] [--scatter LAYER.STAT]
                                    
```

### `browser`

Interactive signal alignment genome browser

This feature is in very early stages. Currently it features trackplot visualization where you can click on different locations to display read/reference information.

```
usage: uncalled4 browser [-h] ref_bounds bam_in [bam_in ...] [--ref REF]
                        [--reads READS [READS ...]] [-x READ_INDEX] [-r] [--rna]
                        [-l READ_FILTER] [-f] [--pore-model PORE_MODEL] [-p BROWSER_PORT]
                        [-o OUTFILE]
                        
```

## Analysis

These functions compute statistics over reference and read coordinates. `refstats` computes summary and comparison statistics (e.g. Kolmogorov-Smirnov test) over reference coordinates, while `layerstats` maintains the read-by-reference dimensions of the DTW alignments.

(Coming soon: `readstats` to compute read-level statistics)

### `refstats`

Calculate per-reference-coordinate statistics

```
usage: uncalled4 refstats [-h] [-t TRACKS] [-R REF_BOUNDS] [--min-coverage MIN_COVERAGE]
                          [--bed-filter BED_FILTER] [--ref-chunksize REF_CHUNKSIZE]
                          [--aln-chunksize ALN_CHUNKSIZE] [-c] [--ref-index REF_INDEX]
                          [-m PORE_MODEL] [-p PROCESSES] [-o OUTFILE]
                          layers refstats bam_in [bam_in ...]

Calculate per-reference-coordinate statistics

positional arguments:
  layers                Comma-separated list of layers over which to compute summary statistics
  refstats              Comma-separated list of summary statistics to compute. Some statisitcs (ks)
                        can only be used if exactly two tracks are provided
                        {var,q95,min,mean,q5,q25,q75,max,stdv,kurt,median,KS,skew}
  bam_in                BAM input file (or "-"/no argument for stdin)

optional arguments:
  -h, --help            show this help message and exit
  -t TRACKS, --tracks TRACKS
                        Names of tracks to read from input(s) (default: None)
  -R REF_BOUNDS, --ref-bounds REF_BOUNDS
                        Only load reads which overlap these coordinates (default: None)
  --min-coverage MIN_COVERAGE
                        Reference positions with less than this coverage will be excluded from each
                        track (or all tracks if shared_refs_only is true) (default: 1)
  --bed-filter BED_FILTER
                        Only parse regions in BED file (default: None)
  --ref-chunksize REF_CHUNKSIZE
                        Number of reference coordinates to query for iteration (default: 10000)
  --aln-chunksize ALN_CHUNKSIZE
                        Number of alignments to query for iteration (default: 500)
  -c, --cov             Output track coverage for each reference position (default: False)
  --ref-index REF_INDEX
                        BWA index prefix (default: None)
  -m PORE_MODEL, --pore-model PORE_MODEL
                        Model preset name or TSV filename (default: None)
  -p PROCESSES, --processes PROCESSES
                        Number of parallel processes (default: 1)
  -o OUTFILE, --outfile OUTFILE

```

## Alignment Layers

Uncalled4 stores signal alignments as a set of **layers** associated with read and reference coordinates. Each layer is derived from the read signal (e.g. `current`), the reference sequence (e.g. `kmer`), or other layers (e.g. `model_diff`). Layers are organized into **layer groups**: `dtw` for signal alignments, `bcaln` for projected basecalled alignments, and `cmp` for alignment comparisons. Several subcommands detailed above take layer names as input, which should generaly be in the form `<group>.<layer>`. For brevity, group can be excluded for `dtw` layers (e.g. you can simply input `current`, not `dtw.currnt`). Below is a table of currently available layers:

| Group | Layer   | Description |
|-------|---------|-------------|
| dtw   | current | Mean read signal current (pA) |
| dtw   | start | Read signal sample start index |
| dtw   | length | Read signal sample length |
| dtw   | dwell   | Signal dwell time (ms/nt, proportional to length) |
| dtw   | model_diff | Difference between predicted (via a pore model) and observed current|
| dtw   | abs_diff | Absolute value of dtw.model_diff |
| dtw   | kmer | Binarized reference k-mer |
| dtw   | events | Number of raw signal events aligned ("stays" > 1, "skips" < 1) |
| dtw   | base | Binarized reference base |
| bcaln   | start | Estimated read signal sample start index |
| bcaln   | length | Estimated read signal sample length |
| bcaln   | indel | Number of inserted (>0) or deleted (<0) nucleotides |
| cmp   | mean_ref_dist | Mean reference distance between two alignments (must first run [`layerstats compare`](#compare)) |
| cmp   | jaccard | Raw sample jaccard distances between two alignments (must first run [`layerstats compare`](#compare)) |

All tracks must be written to the same database for multi-track visualization and analysis (e.g. comparing alignments, calculating KS statistics). You can merge multiple databases into a single file using [`uncalled db merge`](#db)

## Release Notes
- v4.0.0:  Pre-print release
For development history, see https://github.com/skovaka/UNCALLED/tree/dev4
