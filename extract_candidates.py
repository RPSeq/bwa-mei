#!/usr/bin/env python

# for tgi cluster:
#/gapp/x64linux/opt/pythonbrew/venvs/Python-2.7.6/gemini/bin/python
# for uva cluster:

import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter
from string import maketrans

__author__ = "Ryan Smith (ryanpsmith@wustl.edu) with code by Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-12-15 11:43 $"

#main loop function
def extract_clippers(bamfile, is_sam, bam_out, uncompressed_out, anchors, single_fq, pair_fq, clip_len, max_opp_clip=7):
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
    
    # set output file
    if uncompressed_out:
        anchors_out_bam = pysam.Samfile(anchors, 'wbu', template=in_bam)
    elif bam_out:
        anchors_out_bam = pysam.Samfile(anchors, 'wb', template=in_bam)
    else:
        anchors_out_bam = pysam.Samfile(anchors, 'wh', template=in_bam)
        
    pair_fq = open(pair_fq, 'w')
    single_fq = open(single_fq, 'w')
    
    mates = {}
    for al in in_bam:
        # add read name to dictionary if not already there
        #not interested in secondary alignments
        if al.is_secondary or al.is_duplicate:
            continue
        #if we flagged the mate, write the pair and remove from cache
        elif al.qname in mates:
            write_pairs(al, mates[al.qname], anchors_out_bam,single_fq, pair_fq)
            del mates[al.qname]
        #grab all non-proper pairs.
        elif not al.is_proper_pair:
            mates[al.qname] = al
        #might be some signal here too... one end unique, one end repeat, but still flagged as proper_pair
        elif (al.mapq >= 10 and dict(al.tags)['MQ'] < 10) or (al.mapq < 10 and dict(al.tags)['MQ'] >= 10):
            mates[al.qname] = al 
        
        #this case: proper pair, unique mapping. check for clips, if so, flag the read
        elif al.mapq > 0:
            cigar = al.cigar
            #NOTE cigar = tuple list format [(0,14),(4,20)], cigarstr = string format
            #cigar_dict = {0:'M',1:'I',2:'D',3:'N',4:'S',5:'H',6:'P',7:'=',8:'X'}
            #if clipped at least clip_len on L:
            if cigar[0][0] == 4 and cigar[0][1] >= clip_len:
                #if opposite is not clipped more than max opposite clip len, write for realignment
                if cigar[-1][0] != 4 or (cigar[-1][0] == 4 and cigar[-1][1] <= max_opp_clip):
                    write_clip(al, 'L', anchors_out_bam, single_fq)
            #elif clipped at least clip_len on R:
            elif cigar[-1][0] == 4 and cigar[-1][1] >= clip_len:
                #if opposite is not clipped more than max opposite clip len, write for realignment
                if cigar[0][0] != 4 or (cigar[0][0] == 4 and cigar[0][1] <= max_opp_clip):
                    write_clip(al, 'R', anchors_out_bam, single_fq)
                    
    if len(mates) > 0:
        sys.stderr.write("Warning: {0} unmatched mates".format(len(mates)))
            
                
# ============================================
# functions
# ============================================

def reverse_complement(sequence):
    complement = maketrans("ACTGactg", "TGACtgac")
    sequence = sequence[::-1].translate(complement)
    return sequence

def write_fastq(al, fq) :
    seq = al.seq
    quals = al.qual
    if al.is_reverse:
        seq = reverse_complement(seq)
        quals = quals[::-1]
    fq.write("@"+al.qname+" OC:Z:"+al.cigarstring+"\n"+seq+"\n+\n"+quals+"\n")
    
def write_pairs(al1, al2, anchors, single_fq, pair_fq):
    #both reads uniquely mapped
    if al1.mapq > 0 and al2.mapq > 0:
        #realign both
        write_fastq(al1, pair_fq)
        write_fastq(al2, pair_fq)
    elif al1.mapq == 0 and al2.mapq > 0:
        #realign al1
        write_fastq(al1, single_fq)
        anchors.write(al2)
    elif al1.mapq > 0 and al2.mapq == 0:
        #realign al2
        anchors.write(al1)
        write_fastq(al2, single_fq)
    return
    
def write_clip(al, side, anchors, single_fq):
    if side=="L":
        seq = al.seq[:al.qstart]
        quals = al.qual[:al.qstart]
    elif side=="R":
        seq = al.seq[al.qend:]
        quals = al.qual[al.qend:]
    if al.is_reverse:
        seq = reverse_complement(seq)
        quals = quals[::-1]
        
    if al.is_read1:
        al.qname += "_1"
    elif al.is_read2:
        al.qname += "_2"
    anchors.write(al)
    single_fq.write("@"+al.qname+" OC:Z:"+al.cigarstring+"\n"+seq+"\n+\n"+quals+"\n")
    return
    
def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamgroupreads.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Extract candidates for MEI re-alignment")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
    parser.add_argument('-a', '--anchors', required=True, help='Output anchors bamfile')
    parser.add_argument('-s', '--single', required=True, help='Output single-read unanchors fastq')
    parser.add_argument('-p', '--pair', required=True, help='Output paired-read unanchors fastq')
    parser.add_argument('-c', '--clip', required=True, type=int, help='Minimum clip length')
    parser.add_argument('-oc', '--opp_clip', required=False, type=int, help='Max opposite clip length')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    parser.add_argument('-b', required=False, action='store_true', help='Output BAM format')
    parser.add_argument('-u', required=False, action='store_true', help='Output uncompressed BAM format (implies -b)')
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
    extract_clippers(args.input, args.S, args.b, args.u, args.anchors, args.single, args.pair, args.clip, args.opp_clip)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

