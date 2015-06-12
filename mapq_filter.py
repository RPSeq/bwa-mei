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

__author__ = "Ryan Smith (ryanpsmith@wustl.edu) with code by Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-12-15 11:43 $"

def scan(bamfile, is_sam, ref_cutoff, mei_cutoff):
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
    
    out_bam = pysam.Samfile("-", 'wb', template=in_bam)
    d = {}
    for al in in_bam:
        if al.qname not in d:
            d[al.qname] = mei_pair(al, in_bam)
        else:
            d[al.qname].add_anchor(al, in_bam)
        
        if d[al.qname].is_complete():
            d[al.qname].mq_filter(ref_cutoff, mei_cutoff, out_bam)
            del d[al.qname]
            
# ============================================
# functions
# ============================================
class mei_pair:
    def __init__(self, al, bamfile):
        self.ref_anchor = None
        self.mei_anchor = None
        self.complete = False
        self.add_anchor(al, bamfile)
        
    def add_anchor(self, al, bamfile):
        if bamfile.getrname(al.rname).startswith("moblist"):
            self.mei_anchor = al
        else:
            self.ref_anchor = al
        return self.is_complete()
    
    def is_complete(self):
        if self.mei_anchor and self.ref_anchor:
            return True
        return False
    
    def mq_filter(self, ref_cutoff, mei_cutoff, out_bam):
        if (self.ref_anchor.mapq >= ref_cutoff and self.mei_anchor.mapq >= mei_cutoff):
            #self.mei_anchor.mapq = ref_cutoff+5
            out_bam.write(self.ref_anchor)
            out_bam.write(self.mei_anchor)
            return True
        return False
        

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamgroupreads.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Filter name-sorted BWA-MEI discordant or splitters .bam file (default discordants)")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    parser.add_argument('-a', required=True, metavar='REF_MAPQ', type=int, help='Mininum anchor mapq')
    parser.add_argument('-u', required=True, metavar='MEI_MAPQ', type=int, help='Mininum unanchor mapq')
    
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
    scan(args.input, args.S, args.a, args.u)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

