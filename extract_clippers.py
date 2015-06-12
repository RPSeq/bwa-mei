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

def extract_clippers(bamfile, is_sam, anchorfile, unanchorfile, clip_len):
    # set input file
    clip_len = int(clip_len)
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
    
    unanchors_out_fastq = open(unanchorfile, 'w')
    anchors_out_bam = pysam.Samfile(anchorfile, 'wb', template=in_bam)
    
    for al in in_bam:
        # add read name to dictionary if not already there
        if al.is_secondary:
            continue
        #arbitrary mapq cutoff
        if al.mapq > 5:
            #skip if there are secondaries
            #try:
            #    SA = al.opt('SA')
            #    continue
            #except:
            #    pass
            
            #NOTE cigar = tuple list format [(0,14),(4,20)], cigarstr = string format
            #cigar_dict = {0:'M',1:'I',2:'D',3:'N',4:'S',5:'H',6:'P',7:'=',8:'X'}
            cigar = al.cigar
            clip_L = False
            clip_R = False
            #if clipped at least clip_len at al start or end:
            if cigar[0][0] == 4 and cigar[0][1] >= clip_len:
                clip_L = True
            elif cigar[-1][0] == 4 and cigar[-1][1] >= clip_len:
                clip_R = True
            if (clip_L and clip_R) or (not clip_L and not clip_R):
                continue
            if (clip_L or clip_R):
                if clip_L:
                    sequence = al.seq[:al.qstart]
                    quals = al.qual[:al.qstart]
                elif clip_R:
                    sequence = al.seq[al.qend:]
                    quals = al.qual[al.qend:]
                if al.is_reverse:
                    sequence = reverse_complement(sequence)
                    quals = quals[::-1]
                    cigar = cigar[::-1]
                    al.cigar = None
                    al.cigar = cigar
                if al.is_read1:
                    al.qname = al.qname + "_1"
                elif al.is_read2:
                    al.qname = al.qname + "_2"
                unanchors_out_fastq.write("@"+al.qname+" OC:Z:"+al.cigarstring+"\n"+sequence+"\n+\n"+quals+"\n")
                anchors_out_bam.write(al)


# ============================================
# functions
# ============================================

def reverse_complement(sequence):
    complement = maketrans("ACTGactg", "TGACtgac")
    sequence = sequence[::-1].translate(complement)
    return sequence

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamgroupreads.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Group BAM file by read IDs without sorting")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    parser.add_argument('-a', '--anchors', required=True, help='Output anchors bamfile')
    parser.add_argument('-u', '--unanchors', required=True, help='Output unanchors fastq')
    parser.add_argument('-c', '--clip', required=True, help='Minimum clip length')
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
    extract_clippers(args.input, args.S, args.anchors, args.unanchors, args.clip)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

