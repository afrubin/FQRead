from __future__ import print_function
from sys import stderr
import argparse
import os.path
from fqread import read_fastq, split_fastq_path, create_compressed_outfile


def trim_fastq(outdir, files, start, end, length, compression):
    """
    """
    if len(files) == 0:
        print("Error: no files provided", file=stderr)
        return

    elif length is not None:
        length_mode = True
        if start is not None and end is not None:
            if (end - start + 1) != length:
                print("Error: start/end/length do not agree", file=stderr)
                return
        elif end is not None: # start is None, so length should equal end
            if length != end:
                print("Error: length/end do not agree", file=stderr)
                return
        elif start is None:
            start = 1
    else:
        length_mode = False
        if start is None:
            start = 1

    for f in files:
        # open the output file
        _, base, ext, _ = split_fastq_path(f)
        outname = base + ".trim" + ext
        outname = os.path.join(outdir, outname)
        handle = create_compressed_outfile(outname, compression)

        # trim the reads and write them
        for read in read_fastq(f):
            if length_mode:
                read.trim_length(length, start)
            else:
                read.trim(start, end)
            print(read, file=handle)

        handle.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim reads in FASTQ files.")
    parser.add_argument("files", metavar="FQ", nargs="+",
                        help="FASTQ files to trim")
    parser.add_argument("-s", "--start", metavar="N", 
                        help="start of trim region", type=int)
    parser.add_argument("-e", "--end", metavar="N", 
                        help="end of trim region", type=int)
    parser.add_argument("-l", "--length", metavar="N", 
                        help="length of trim region", type=int)
    parser.add_argument("-o", "--output", default=".", metavar="DIR",
                        help="output directory")
    parser.add_argument("--bz2", dest="compression", action="store_const", 
                        const="bz2", default=None,
                        help="compress output with bzip2")
    parser.add_argument("--gz", dest="compression", action="store_const", 
                        const="gz", default=None,
                        help="compress output with gzip")

    args = parser.parse_args()

    trim_fastq(args.output, args.files, args.start, args.end, args.length, args.compression)
    