# Uncalled4 Utility Scripts

## bamprep.py
Run on basecalled read alignments to label multimapping reads with "mv:" and "ts:" tags required for Uncalled4 alignment, which are otherwise only included on primary alignments. 

```
usage: bamprep.py [-h] [-e EXT] [-m] [-t THREADS] [-o OUTFILE] infiles [infiles ...]

Copy BAM tags from primary to supplemental and secondary alignents. Useful to label all reads
with basecaller move tags. Multiple files and directories can be included, and output
optionally merged at the end

positional arguments:
  infiles               Input BAM files or directories containing BAM files

optional arguments:
  -h, --help            show this help message and exit
  -e EXT, --ext EXT     Only files with this extension will be included in directory search
  -m, --merge           Merge output files at the end. Creates a temporary directory
  -t THREADS, --threads THREADS
                        Number of threads to use ONLY for merging
  -o OUTFILE, --outfile OUTFILE
                        Output directory, or file if --merge specified
```

Most efficently run by directly inputting basecaller paths and outputting a single merged BAM file: 

`python3 bamprep.py -m pass/ fail/ -o all_tagged.bam`

Currently large individual BAM files are slow and memory-intensive, relying on [pysam.IndexedReads](https://pysam.readthedocs.io/en/latest/api.html#pysam.IndexedReads) to find primary alignments. Using name-sorted BAMs could allow for more efficent processing of large files, but the provided implementation is most efficent for Guppy's output, which consists of many small coordinate-sorted files. 

## t2g.py

This script converts transcriptome to genome coordinates in a tabular format (e.g. TSV) guided by a GTF annotation. Currently, the file must contain a header labeling the column names.

```
usage: t2g.py [-h] [-t TRANS_COL] [-p POS_COL] [-d DELIM] [-o OUTFILE] gtf infile

Translates transcriptome to genome coordinates

positional arguments:
  gtf                   GTF annotation file
  infile                Tab-delimited file (or set by --delim) containing transcript
                        names and coordinates

optional arguments:
  -h, --help            show this help message and exit
  -t TRANS_COL, --trans-col TRANS_COL
                        Name of column with transcript IDs (matching GTF transcript_id)
  -p POS_COL, --pos-col POS_COL
                        Name of column with transcript coordinates
  -d DELIM, --delim DELIM
                        Field seperator
  -o OUTFILE, --outfile OUTFILE
                        Output file
```

This can be used to translate m6Anet output to genome coordinates by specifying ',' as the field delimiter: 

`python3 t2g.py -d , anno.gtf data.site_proba.csv > data.site_proba.t2g.csv`

In the Uncalled4 manuscript, we average the translated results to the gene level using Pandas groupby:

```
import pandas as pd
df_trans = pd.read_csv("data.site_proba.t2g.csv")
grp = df_trans.groupby(["gene_id","chr_loc"])
df_gene = pd.DataFrame({
	"probability_modified" : grp["probability_modified"].mean(),
	"mod_ratio" : grp["mod_ratio"].mean(),
	"n_reads" : grp["n_reads"].mean(),
	"kmer" : grp["kmer"].first(),
})
return df_gene
```
