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
__date__ = "$Date: 2015-01-27 15:52 $"

def fixpairflags(bamfile):
    # set input file
    in_bam = pysam.Samfile(bamfile, "rb")

    # set output file
    out_bam = pysam.Samfile("-", 'wh', template=in_bam)


    d = {}
    for al in in_bam:
        # add read name to dictionary if not already there
    	key = al.qname
        if key not in d:
            d.setdefault(key,Namegroup(al))
        # print matched read pairs
        else:
            d[key].add_alignment(al)
            if d[key].is_complete():
            	anchor = False
                unanchor = False
                for al in d[key].alignments:
                # check if secondary (don't use secondaries) and if anchored/unanchored end
                	if al.is_secondary:
                		continue

                	# THE FOLLOWING TWO ELIF BLOCKS SHOULD REPAIR THE FLAG AND MATE INFO FOR PAIRS
                	# if al.is_paired, this is the anchored end.
                	if al.is_paired:
                		anchor = al

                	# if al is not paired, this is the remapped end.
                	elif not al.is_paired:
                		unanchor = al

                	if anchor and unanchor:
                		unanchor.is_paired = True
                		if anchor.is_read1:
                			unanchor.is_read2 = True
                 		elif anchor.is_read2:
                 			unanchor.is_read1 = True
                 		unanchor.mate_is_reverse = anchor.is_reverse
                 		unanchor.rnext = anchor.rname
                 		unanchor.pnext = anchor.pos
                 		anchor.mate_is_reverse = unanchor.is_reverse
                 		anchor.rnext = unanchor.rname
                 		anchor.pnext = unanchor.pos
                 		out_bam.write(anchor)
                 		out_bam.write(unanchor)
                #this might miss some groups.
                # i use this script after heavy BAM filtering.
                # all secondary reads might not be present,
                # which is how .is_complete determines a complete group.
                del d[key]

    #SA tag won't necessarily work, so check remaining groups
    for group in d:
        for al in d[group].alignments:
            # check if secondary (don't use secondaries) and if anchored/unanchored end
            if al.is_secondary:
                continue

            # THE FOLLOWING TWO ELIF BLOCKS SHOULD REPAIR THE FLAG AND MATE INFO FOR PAIRS
            # if al.is_paired, this is the anchored end.
            if al.is_paired:
                anchor = al

            # if al is not paired, this is the remapped end.
            elif not al.is_paired:
                unanchor = al

            if anchor and unanchor:
                unanchor.is_paired = True
                if anchor.is_read1:
                    unanchor.is_read2 = True
                elif anchor.is_read2:
                    unanchor.is_read1 = True
                unanchor.mate_is_reverse = anchor.is_reverse
                unanchor.rnext = anchor.rname
                unanchor.pnext = anchor.pos
                anchor.mate_is_reverse = unanchor.is_reverse
                anchor.rnext = unanchor.rname
                anchor.pnext = unanchor.pos
                out_bam.write(anchor)
                out_bam.write(unanchor)

    #if len(d) != 0:
        #sys.stderr.write('Warning: %s unmatched name groups\n' % len(d))

# ============================================
# functions
# ============================================

# class that holds reads from a sequence fragment
class Namegroup():
    def __init__(self, al):
        self.alignments = list()
        self.name = al.qname
        self.sa = 0
        self.num_prim = 0
        self.add_alignment(al)

    def add_alignment(self, al):
        self.alignments.append(al)
        if not al.is_secondary:
            self.num_prim += 1
            try:
                self.sa += len(al.opt('SA').rstrip(';').split(';'))
                # print self.sa
            except KeyError:
                pass

    def is_complete(self):
        return self.num_prim == 2 and len(self.alignments) == self.sa + 2

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamgroupreads.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Group BAM file by read IDs without sorting")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
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
    fixpairflags(args.input)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

