#!/usr/bin/env python

# for tgi cluster:
#/gapp/x64linux/opt/pythonbrew/venvs/Python-2.7.6/gemini/bin/python
# for uva cluster:

import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter
import string
from string import *
from collections import defaultdict

__author__ = "Ryan Smith (ryanpsmith@wustl.edu) with code by Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-12-15 11:43 $"

def filter_realigned_pairs(bamfile, anchorfile, is_sam, out_bam, out_anchors):
    # set input file
    if bamfile == None:
        if is_sam:
            in_bam = pysam.Samfile("-", "r")
        else:
            in_bam = pysam.Samfile('-', 'rb')
    else:
        if is_sam:
            in_bam = pysam.Samfile(bamfile, 'r')
        else:
            in_bam = pysam.Samfile(bamfile, "rb")
        
    anchors_in = pysam.Samfile(anchorfile, 'r')
    out_bam = pysam.Samfile(out_bam, 'wb', template=in_bam)
    out_anchors = pysam.Samfile(out_anchors, 'wb', template=anchors_in)

    prev_al = False
    d = defaultdict(list)
    for al in in_bam:
        namelist = al.qname.split("_")
        if len(namelist) == 2:
            name = namelist[0]
            num = namelist[1]
            d[name].append(al)

    for al in anchors_in:
        name, num = al.qname.split("_")
        if name in d:
            if len(d[name]) > 1:
                del d[name]
                continue
            mei_al = d[name][0]
            namelist = mei_al.qname.split("_")
            if len(namelist) == 2:
                mei_name = namelist[0]
                mei_num = namelist[1]
            if num != mei_num and mismatch(mei_al) <= 0.1:
                al.qname = name
                mei_al.qname = name
                out_anchors.write(al)
                out_bam.write(mei_al)
                del d[name]

        

# ============================================
# functions
# ============================================
def mismatch(al):
    cigar = al.cigar
    #these flags are all mismatches (I,D,N,S,H,X)
    mismatches = set([1,2,3,4,5,8])
    NM = sum([sgmt[1] for sgmt in cigar if sgmt[0] in mismatches])
    length = len(al.seq)
    percent = NM/(float(length))
    return percent

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamgroupreads.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Group BAM file by read IDs without sorting")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    parser.add_argument('-a', required=True, help='Input pairs anchor file')
    parser.add_argument('-o', required=True, help='Output MEI BAM')
    parser.add_argument('-ao', required=True, help='Output anchors BAM file')
    # parse the arguments
    args = parser.parse_args()
    
    # bail if no BAM file
    if args.input is None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
    
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
    filter_realigned_pairs(args.input, args.a, args.S, args.o, args.ao)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

