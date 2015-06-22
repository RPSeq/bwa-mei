#!/usr/bin/env python

import pysam
import sys
import argparse
from argparse import RawTextHelpFormatter
from string import maketrans
from multiprocessing import Process, Pool, Manager, Queue, Lock

__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-6-21 11:43 $"


#main
def main():

    #get command line args
    args = get_args() 

    #get IO streams
    in_bam = get_sam_IO(args.input, args.S)
    #define the output streams

    #Open read and write queues
    readQueue = Queue()
    writeQueue = Queue()

    workers = args.t
    if workers==1:
        mates = dict()
    elif workers>1:
        manager = Manager()
        mates = manager.dict()
    else:
        exit("Must have at least 1 thread")


    #Create reader and writer processes
    reader = Process(target=bam_reader, args=((in_bam, 100, readQueue, workers)))
    writer = Process(target=bam_writer, args=((writeQueue, workers, args.anchors, args.single, args.pair)))

    #Start them first
    reader.start()
    writer.start()

    #Start pool of worker processes
    worker = Pool(workers, al_scanner, ((readQueue, writeQueue, mates, args.clip, args.opclip)))


    reader.join()

    worker.close()
    worker.join()

    writer.join()

    if len(mates) > 0:
         sys.stderr.write("Warning: {0} unmatched mates".format(len(mates)))

# ============================================
# functions
# ============================================

#opens all output streams and returns the file objects
def get_sam_IO(bamfile, is_sam):

    # set input BAM file
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
    
    #now we are manually writing SAM (possible uncompressed BAM too)

    #return output file objects for writing
    return in_bam
               
#main worker function
#REALLY need to stop using the mates dict, this slows down the workers badly
def al_scanner(readQueue, writeQueue, mates, clip_len, max_opp_clip):
    while True:
        als = readQueue.get()
        if als == 'DONE':
            writeQueue.put('DONE')
            return
        for al in als:
            # add read name to dictionary if not already there
            #not interested in secondary alignments
            if al.is_secondary or al.is_duplicate:
                continue
            #if we flagged the mate, write the pair and remove from cache
            elif al.qname in mates:
                writeQueue.put(('write_pairs', al, mates[al.qname]))
                #write_pairs(al, mates[al.qname], anchors_out,single_fq, pair_fq)
                del mates[al.qname]
            #grab all non-proper pairs.
            elif not al.is_paired and not al.is_proper_pair:
                mates[al.qname] = al
            #might be some signal here too... one end unique, one end repeat, but still flagged as proper_pair
            elif (al.mapq >= 10 and al.tags['MQ'] < 10) or (al.mapq < 10 and al.tags['MQ'] >= 10):
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
                        writeQueue.put(('write_clip', al, 'L'))
                        #write_clip(al, 'L', anchors_out, single_fq)
                #elif clipped at least clip_len on R:
                elif cigar[-1][0] == 4 and cigar[-1][1] >= clip_len:
                    #if opposite is not clipped more than max opposite clip len, write for realignment
                    if cigar[0][0] != 4 or (cigar[0][0] == 4 and cigar[0][1] <= max_opp_clip):
                        #write_clip(al, 'R', anchors_out, single_fq)
                        writeQueue.put(('write_clip', al, 'R'))


def bam_reader(in_bam, chunk, readQueue, workers):
    als = []
    i = 0
    for al in in_bam:
        if al.is_duplicate or al.is_secondary:
            continue
        al = sam_al(al, in_bam)
        als.append(al)
        i+=1
        if i==chunk:
            readQueue.put(als)
            als = []
            i = 0
    readQueue.put(als)
    #each worker must recieve a 'DONE' token in order to exit
    for i in range(workers):
        readQueue.put('DONE')
    return
    
#file output worker
def bam_writer(writeQueue, workers, anchors_out, single_fq, pair_fq):
    anchors = open(anchors_out, 'w+')
    single_fq = open(single_fq, 'w+')
    pair_fq = open(pair_fq, 'w+')
    finished = 0
    while finished < workers:
        writelist = writeQueue.get()
        if writelist == 'DONE':
            finished+=1
            continue
        if writelist[0] == 'write_clip':
            write_clip(writelist[1:], anchors, single_fq)
        elif writelist[0] == 'write_pairs':
            write_pairs(writelist[1:], anchors, single_fq, pair_fq)
    anchors.close()
    single_fq.close()
    pair_fq.close()
    return

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
    
def write_pairs(writelist, anchors, single_fq, pair_fq):
    al1, al2 = writelist
    #both reads uniquely mapped
    if al1.mapq > 0 and al2.mapq > 0:
        #realign both
        write_fastq(al1, pair_fq)
        write_fastq(al2, pair_fq)
    elif al1.mapq == 0 and al2.mapq > 0:
        #realign al1
        write_fastq(al1, single_fq)
        al2.write(anchors)
    elif al1.mapq > 0 and al2.mapq == 0:
        #realign al2
        al1.write(anchors)
        write_fastq(al2, single_fq)
    return
    
def write_clip(writelist, anchors, single_fq):
    al, side = writelist
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
    al.write(anchors)
    single_fq.write("@"+al.qname+" OC:Z:"+al.cigarstring+"\n"+seq+"\n+\n"+quals+"\n")
    return
    
def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
                                            extract_candidates.py\n\
                                            author: " + __author__ + "\n\
                                            version: " + __version__ + "\n\
                                            description: Extract candidates for MEI re-alignment")
    parser.add_argument('-i', '--input', metavar='BAM/SAM', required=False, help='Input file')
    parser.add_argument('-a', '--anchors', metavar='SAM', required=True, help='Output anchors samfile')
    parser.add_argument('-s', '--single', metavar='FQ', required=True, help='Output single-read unanchors fastq')
    parser.add_argument('-p', '--pair', metavar='FQ', required=True, help='Output paired-read unanchors fastq')
    parser.add_argument('-c', '--clip', metavar='LEN', required=True, type=int, help='Minimum clip length')
    parser.add_argument('-oc', '--opclip', metavar='LEN', required=False, type=int, help='Max opposite clip length')
    parser.add_argument('-S', required=False, action='store_true', help='Input is SAM format')
    parser.add_argument('-t', required=True,type=int, help='Number of processing threads')

    # parse the arguments
    args = parser.parse_args()
    
    # bail if no BAM file
    if args.input is None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
    
    # send back the user input
    return args


# ====================
# SAM Class
# ====================
class sam_al(object):
    def __init__(self, sam, in_sam=False):
        #case: input is pysam al
        if type(sam)==pysam.AlignedRead and type(in_sam)==pysam.Samfile:
                self.read_pysam(sam, in_sam)
        elif type(sam)==str:
            sam = sam.strip("\n").split("\t")
            if len(sam) < 11:
                exit("Error, sam_al input is not complete SAM entry")
        elif type(sam)!=list:
            exit("Error, sam_al() requires list, str, or pysam.Samfile objects as input")
        else:
            exit("Error: sam_al requires pysam.Samfile AND pysam.AlignedRead")
        return 

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
            tags[tag]=(ttype, val)
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

    def write(self, output):
        outlist = [self.qname, str(self.flag), self.rname, 
                    str(self.pos), str(self.mapq), self.cigarstring, 
                    self.rnext, str(self.pnext), str(self.tlen), self.seq, self.qual]

        for tag, val in self.tags.viewitems():
            if type(val)==int:
                ttype = 'i'
            else:
                ttype = 'Z'
            outlist.append("{0}:{1}:{2}".format(tag, ttype, val))
        output.write("\t".join(outlist)+"\n")

# ============================================
# driver
# ============================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise

