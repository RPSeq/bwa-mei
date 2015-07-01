#!/usr/bin/env python

# for tgi cluster:
#/gapp/x64linux/opt/pythonbrew/venvs/Python-2.7.6/gemini/bin/python
# for uva cluster:

import pysam, sys, argparse, string
from argparse import RawTextHelpFormatter
#from string import *

__author__ = "Ryan Smith (ryanpsmith@wustl.edu) with code by Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.2 $"
__date__ = "$Date: 2015-02-09 23:44 $"

# ============================================
# functions
# ============================================
def repair_cigars(mei_anchor, ref_anchor, clips):

    ##get original cigar from OC tag
    #try:
    #    old_mei_cigarstr = dict(mei_anchor.tags)['OC']
    #except:
    #    exit("Error: Input BAM must have OC:Z:Cigar (original cigar) tag")
        
    #get cigar of ref alignment
    ref_cigar = ref_anchor.cigar
    #if ref al is reverse, cigar still corresponds to + strand. flip cigar to match query 5' and 3'
    if ref_anchor.is_reverse:
        ref_cigar = ref_cigar[::-1]

    #get new partial cigar from mei realignment
    mei_cigar = mei_anchor.cigar

    #NOTE pysam al.cigar = tuple list format [(0,14),(4,20)], al.cigarstring = string format
    #cigar_dict = {0:'M',1:'I',2:'D',3:'N',4:'S',5:'H',6:'P',7:'=',8:'X'}
    
    #if ref is clipped at 5' end
    if ref_cigar[0][0] == 4 and ref_cigar[0][1] >= clips:
        #get length of rest of ref (from next cigar tuple to end of cigar)
        length = sum([tupe[1] for tupe in ref_cigar[1:]])
        #add that length as H to mei cigar
        #if ref is reversed, add to end of MEI cigar
        if ref_anchor.is_reverse:
            mei_cigar.append((5,length))
        #if ref is not reversed, add to start of MEI cigar
        else:
            mei_cigar.insert(0, (5, length))
    
    #if ref is clipped at 3' end
    elif ref_cigar[-1][0] == 4 and ref_cigar[-1][1] >= clips:
        #get length of rest of ref (from first cigar tuple to second to last)
        length = sum([tupe[1] for tupe in ref_cigar[:-1]])
        #add that length as H to mei cigar
        #if ref is reverse, add to start of MEI cigar
        if ref_anchor.is_reverse:
            mei_cigar.insert(0, (5, length))
        #if ref if not reverse, add to end of MEI cigar
        else:
            mei_cigar.append((5,length))
    
    #if ref and mei anchor are same orientation, must flip mei_cigar because we also flipped the ref_cigar
    if ref_anchor.is_reverse == mei_anchor.is_reverse:
        mei_cigar = mei_cigar[::-1]

    #reset cigar to new value and update flags
    mei_anchor.cigar = None
    mei_anchor.cigar = mei_cigar
    mei_anchor.is_secondary = True
    mei_anchor.is_paired = True
    if ref_anchor.is_read1:
        mei_anchor.is_read1 = True
    elif ref_anchor.is_read2:
        mei_anchor.is_read2 = True
    return mei_anchor


def groupsplitters(bamfile, clips):
    # set input file
    in_bam = pysam.Samfile(bamfile, "rb")

    # check name sorted
    try:
        if in_bam.header['HD']['SO'] != 'queryname':
            exit("Error: File must be name sorted")
    except:
        exit("Error: File must be name sorted")

    # set output file
    out_bam = pysam.Samfile("-", 'wh', template=in_bam)
    
    #starting with second alignment...
    for al_1 in in_bam:

        al_2 = in_bam.next()

        #MEI anchor has OC tag
        if 'moblist' in in_bam.getrname(al_1.rname):
            al_1 = repair_cigars(al_1, al_2, clips)

        elif 'moblist' in in_bam.getrname(al_2.rname):
            al_2 = repair_cigars(al_2, al_1, clips)

        out_bam.write(al_1)
        out_bam.write(al_2)

    in_bam.close()
    out_bam.close()


def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamgroupreads.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Group BAM file by read IDs without sorting")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
    parser.add_argument('-o', '--output', metavar='OUT', required=False, help='Output BAM file')
    parser.add_argument('-c', '--clips', metavar='CLIPLEN', type=int, required=False, default=False, help='Input is clipped read BAM')
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
    groupsplitters(args.input, args.clips)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

