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

# ====================
# SAM Class
# ====================
class sam_al(object):
    '''Class representing a SAM file alignment entry'''

    def __init__(self, sam, in_sam=False):
        #manual overloading based on arg types
        if type(sam) == pysam.AlignedRead and in_sam:
            self.read_pysam(sam, in_sam)
        elif type(sam) == str or type(sam) == list:
            self.read(sam)
        else:
            exit("Error creating sam_al.\nUsage:sam_al(samlist), \
                sam_al(samstr), sam_al(pysam.al, pysam.infile)")

    def read(self, sam):
        if type(sam)==str:
            sam = sam.strip("\n").split("\t")
        self.qname = sam[0]
        self.flag = int(sam[1])
        self.rname = sam[2]
        self.pos = int(sam[3])
        self.mapq = int(sam[4])
        self.cigar = sam[5]
        self.rnext = sam[6]
        self.pnext = sam[7]
        self.tlen = int(sam[8])
        self.seq = sam[9]
        self.qual = sam[10]
        self.tags = {}
        for i in range(10,len(sam)):
            tag, ttype, val = sam[i].split(":")
            tags[tag]=(val)
        return

    def read_pysam(self, al, in_sam):
        self.qname = al.qname
        self.flag = al.flag
        self.rname = in_sam.getrname(al.tid)
        self.pos = al.pos
        self.mapq = al.mapq
        self.cigarstring = al.cigarstring
        self.cigar = al.cigar
        self.rnext = in_sam.getrname(al.mrnm)
        self.pnext = al.pnext
        self.tlen = al.tlen
        self.seq = al.seq
        self.qual = al.qual
        self.tags = dict(al.tags)
        self.is_secondary = al.is_secondary
        self.is_duplicate = al.is_duplicate
        self.is_proper_pair = al.is_proper_pair
        self.is_read1 = al.is_read1
        self.is_read2 = al.is_read2
        self.qstart = al.qstart
        self.qend = al.qend
        self.is_reverse = al.is_reverse
        self.is_paired = al.is_paired
        return

    def sam_str(self, tag=False):
        name = self.qname
        if tag:
            if self.is_read1:
                name = self.qname+"_1"
            elif self.is_read2:
                name = self.qname+"_2"
        outlist = [name, str(self.flag), self.rname, 
                    str(self.pos), str(self.mapq), self.cigarstring, 
                    self.rnext, str(self.pnext), str(self.tlen), self.seq, self.qual]

        for tag, val in self.tags.viewitems():
            if type(val)==int:
                ttype = 'i'
            else:
                ttype = 'Z'
            outlist.append("{0}:{1}:{2}".format(tag, ttype, val))
        return "\t".join(outlist)+"\n"

    def fastq(self):
        seq = self.seq
        quals = self.qual
        if self.is_reverse:
            seq = reverse_complement(seq)
            quals = quals[::-1]
        fq = "@"+self.qname+" OC:Z:"+self.cigarstring+"\n"+seq+"\n+\n"+quals+"\n"
        return fq

    def write(self, output):
        line = self.sam_str()
        output.write(line)
        return

#main loop function
def extract_clippers(bamfile, is_sam, pair_anchors, pair_fastq, clip_anchors, clip_fastq, clip_len, max_opp_clip=7):
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

    header = "@HD\tVN:1.3\tSO:unsorted\n"
    header+="\n".join(in_bam.text.split("\n")[1:])
    
    pair_anchors_out = open(pair_anchors, 'w')
    pair_fastq_out = open(pair_fastq, 'w')

    clip_anchors_out = open(clip_anchors, 'w')
    clip_fastq_out = open(clip_fastq, 'w')

    pair_anchors_out.write(header)
    clip_anchors_out.write(header)

    mates = {}
    pair_anchor_batch = []
    clip_anchor_batch = []
    pair_fq_batch = []
    clip_fq_batch = []

    batchsize = 1000000
    for al in in_bam:
        # add read name to dictionary if not already there
        #not interested in secondary alignments
        if len(pair_anchor_batch) >= batchsize:
            pair_anchors_out.write("".join(pair_anchor_batch))
            del pair_anchor_batch[:]
            #pair_anchor_batch = []
        if len(pair_fq_batch) >= batchsize:
            pair_fastq_out.write("".join(pair_fq_batch))
            del pair_fq_batch[:]
            #pair_fq_batch = []
        if len(clip_fq_batch) >= batchsize:
            clip_fastq_out.write("".join(clip_fq_batch))
            del clip_fq_batch[:]
            #clip_fq_batch = []
        if len(clip_anchor_batch) >= batchsize:
            clip_anchors_out.write("".join(clip_anchor_batch))
            del clip_anchor_batch[:]
            #clip_anchor_batch = []
        rname = in_bam.getrname(al.rname)
        mrname = in_bam.getrname(al.mrnm)
        if al.is_secondary or al.is_duplicate or rname.startswith("GL") or mrname.startswith("GL"):
            continue

                #if mate mapq suggests unique, look for clips
        if al.mapq > 0:
            cigar = al.cigar
            if len(cigar) < 1:
                pass
            #NOTE cigar = tuple list format [(0,14),(4,20)], cigarstr = string format
            #cigar_dict = {0:'M',1:'I',2:'D',3:'N',4:'S',5:'H',6:'P',7:'=',8:'X'}
            #if clipped at least clip_len on L:
            elif cigar[0][0] == 4 and cigar[0][1] >= clip_len:
                #if opposite is not clipped more than max opposite clip len, write for realignment
                if cigar[-1][0] != 4 or (cigar[-1][0] == 4 and cigar[-1][1] <= max_opp_clip):
                    fq, anc = write_clip(al, 'L', in_bam)
                    clip_fq_batch.append(fq)
                    clip_anchor_batch.append(anc)

            #elif clipped at least clip_len on R:
            elif cigar[-1][0] == 4 and cigar[-1][1] >= clip_len:
                #if opposite is not clipped more than max opposite clip len, write for realignment
                if cigar[0][0] != 4 or (cigar[0][0] == 4 and cigar[0][1] <= max_opp_clip):
                    fq, anc = write_clip(al, 'R', in_bam)
                    clip_fq_batch.append(fq)
                    clip_anchor_batch.append(anc)

#THIS SHOULD BE CHANGED SO USES ZSCORE 3 DIFFERENCE FROM EXPECTED TEMPLATE LENGTH FOR EACH READ GROUP
#SHOULD ALSO JUST USE ONE IF (cond1 | cond2 | cond3):
        # if not al.is_proper_pair:
        #     fq, anc = write_pairs(al, in_bam)
        #     if fq and anc:
        #         pair_fq_batch.append(fq)
        #         pair_anchor_batch.append(anc)
        #     elif anc:
        #         pair_anchor_batch.append(anc)

        #mate mapped to diff chrom 
        if al.rname != al.mrnm:
            #print("rname")
            fq, anc = write_pairs(al, in_bam)
            if fq and anc:
                pair_fq_batch.append(fq)
                pair_anchor_batch.append(anc)
            elif anc:
                pair_anchor_batch.append(anc)

        #mate wrong orientation
        elif al.is_reverse == al.mate_is_reverse:
            #print("reverse")
            fq, anc = write_pairs(al, in_bam)
            if fq and anc:
                pair_fq_batch.append(fq)
                pair_anchor_batch.append(anc)
            elif anc:
                pair_anchor_batch.append(anc)

        #one mapped, one unmapped
        elif al.is_unmapped != al.mate_is_unmapped:
            #print("unmapped")
            fq, anc = write_pairs(al, in_bam)
            if fq and anc:
                pair_fq_batch.append(fq)
                pair_anchor_batch.append(anc)
            elif anc:
                pair_anchor_batch.append(anc)

    pair_anchors_out.write("".join(pair_anchor_batch))
    pair_fastq_out.write("".join(pair_fq_batch))
    clip_fastq_out.write("".join(clip_fq_batch))
    clip_anchors_out.write("".join(clip_anchor_batch))
                    
    if len(mates) > 0:
        sys.stderr.write("Warning: {0} unmatched mates\n".format(len(mates)))
            
# ============================================
# functions
# ============================================

def reverse_complement(sequence):
    complement = maketrans("ACTGactg", "TGACtgac")
    sequence = sequence[::-1].translate(complement)
    return sequence

def write_fastq(al) :
    seq = al.seq
    quals = al.qual
    name = al.qname
    if al.is_reverse:
        seq = reverse_complement(seq)
        quals = quals[::-1] 
    if al.is_read1:
        name += "_1" 
    elif al.is_read2:
        name += "_2"
    return "@"+name+" OC:Z:"+al.cigarstring+"\n"+seq+"\n+\n"+quals+"\n"

def write_pairs(al1, in_bam):
    #both reads uniquely mapped
    mate_mapq = dict(al1.tags)['MQ']
    if al1.mapq > 0 and mate_mapq > 0:
        #realign both
        fq = write_fastq(al1)
        anc = sam_al(al1, in_bam).sam_str(1)
        return fq, anc
    elif al1.mapq == 0 and mate_mapq > 0:
        #realign al1
        fq = write_fastq(al1)
        anc = sam_al(al1, in_bam).sam_str(1)
        return fq, anc
    elif al1.mapq > 0 and mate_mapq == 0:
        #dont realign this mate
        #fq = write_fastq(al2)
        anc = sam_al(al1, in_bam).sam_str(1)
        return False, anc
    else:
        return False, False
    
def write_clip(al, side, in_bam):
    name = al.qname
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
        name += "_1"
    elif al.is_read2:
        name += "_2"
    anc = sam_al(al, in_bam).sam_str(1)
    fq = "@"+name+" OC:Z:"+al.cigarstring+"\n"+seq+"\n+\n"+quals+"\n"
    return fq, anc

def zscore(val, mean, stdev):
    return abs((float(val)-float(mean))/float(stdev))

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
bamgroupreads.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Extract candidates for MEI re-alignment")
    parser.add_argument('-i', '--input', metavar='BAM', required=False, help='Input BAM file')
    parser.add_argument('-pa', '--pair_anchors', metavar='SAM', required=True, help='Output paired anchors SAM')
    parser.add_argument('-pf', '--pair_fastq', metavar='FASTQ', required=True, help='Output paired realign FASTQ')
    parser.add_argument('-ca', '--clip_anchors', metavar='SAM', required=True, help='Output clip anchors SAM')
    parser.add_argument('-cf', '--clip_fastq', metavar='FASTQ', required=True, help='Output clip realign FASTQ')
    parser.add_argument('-c', '--clip', metavar='LEN', required=True, type=int, help='Minimum clip length')
    parser.add_argument('-oc', '--opclip', metavar='LEN', required=False, type=int, help='Max opposite clip length')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    # parser.add_argument('-m', metavar='MEAN ILEN', required=True, type=float, help='Mean insert size')
    # parser.add_argument('-sd', metavar='STDEV ILEN', required=True, type=float, help='Insert size STDEV')
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
    #extract_clippers(bamfile, is_sam, pair_anchors, pair_fastq, clip_anchors, clip_fastq, clip_len, max_opp_clip=7)
    extract_clippers(args.input, args.S, args.pair_anchors, args.pair_fastq, args.clip_anchors, args.clip_fastq, args.clip, args.opclip)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

