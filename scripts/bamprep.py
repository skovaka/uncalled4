import pysam
import sys
from time import time
import os

import argparse

parser = argparse.ArgumentParser(description="Copy BAM tags from primary to supplemental and secondary alignents. Useful to label all reads with basecaller move tags. Multiple files and directories can be included, and output optionally merged at the end")
parser.add_argument("infiles", type=str, nargs="+", help="Input BAM files or directories containing BAM files")
#parser.add_argument("-t", "--tags", default=["mv","ts"], nargs="+")
parser.add_argument("-e", "--ext", default=".bam", help="Only files with this extension will be included in directory search")
parser.add_argument("-m", "--merge", action="store_true", help="Merge output files at the end. Creates a temporary directory")
parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use ONLY for merging")
parser.add_argument("-o", "--outfile", required=False, type=str, default=None, help="Output directory, or file if --merge specified")

def main(args):

    if args.merge:
        out_dir = args.outfile + "_tmp"
    else:
        out_dir = args.outfile

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    out_paths = list()

    for path in args.infiles:
        if os.path.isdir(path):
            files = [
                os.path.join(path,fname) 
                for fname in os.listdir(path) 
                if fname.endswith(args.ext) ]
        else: 
            files = [path]

        for in_path in files:
            out_path = os.path.join(out_dir, in_path.replace("/","."))
            copy_tags(in_path, out_path)
            out_paths.append(out_path)

    if args.merge:
        print("Merging...")
        pysam.merge("-@", f"{args.threads}", args.outfile, *out_paths)

def copy_tags(fname_in, fname_out):
    #try:
    iterfile = pysam.AlignmentFile(fname_in, "rb")
    idxfile = pysam.AlignmentFile(fname_in, "rb")
   # except:
   #     sys.stderr.write(f"Error: {fname_in} failed to parse\n")
   #     return
    outfile = pysam.AlignmentFile(fname_out, "wb", template=iterfile)

    idx = pysam.IndexedReads(idxfile, multiple_iterators=True)
    idx.build()

    for aln in iterfile:
        if not aln.has_tag("mv"):
            found = False
            for other in idx.find(aln.query_name):
                if other.has_tag("mv"):
                    for name,val,dtype in other.get_tags(True):
                        if not aln.has_tag(name):
                            aln.set_tag(name,val)
                    found = True
            if not found:
                sys.stderr.write(f"Warning: 'mv' tag not found for read {aln.query_name}\n")
        outfile.write(aln)

    iterfile.close()
    idxfile.close()
    outfile.close()

if __name__ == "__main__":
    args = parser.parse_args(sys.argv[1:])
    main(args)
