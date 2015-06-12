import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter

__author__ = "Ryan Smith (ryanpsmith@wustl.edu) with code by Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-1-23 12:18 $"

def bamremovenonpairs(bamfile, ids_file):

    # set input files
    in_bam = pysam.Samfile(bamfile, "rb")

    mapped_ids = set([line.strip() for line in open(ids_file)])

    # set output file
    out_bam = pysam.Samfile('-', 'wb', template=in_bam)

    # Write out matching alignments
    for al in in_bam:
    	if al.qname in mapped_ids:
    		out_bam.write(al)

    in_bam.close()
    out_bam.close()

# ============================================
# functions
# ============================================


def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamgroupreads.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Pulls anchored reads paired to transposon alignments")
    parser.add_argument('-a', required=True, help='BAM file to be filtered')
    parser.add_argument('-b', required=True, help='list of hits')

    # parse the arguments
    args = parser.parse_args()

    # send back the user input
    return args

# ============================================
# driver
# ============================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():
	args = get_args()
	bamremovenonpairs(args.a, args.b)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise
