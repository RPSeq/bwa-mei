#!/usr/bin/env python

import argparse, sys
import math, time, re
from argparse import RawTextHelpFormatter
from collections import defaultdict

__author__ = "Ryan Smith (ryan.smith.p@gmail.com)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2015-06-16 10:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
vcfToBedpe\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Convert a VCF file to a BEDPE file")
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default=None, help='VCF input (default: stdin)')
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout, help='Output BEDPE to write (default: stdout)')
    
    # parse the arguments
    args = parser.parse_args()
    
    # if no input, check if part of pipe and if so, read stdin.
    if args.input == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input = sys.stdin
    
    # send back the user input
    return args

class Vcf(object):
    def __init__(self):
        self.file_format = 'VCFv4.2'
        # self.fasta = fasta
        self.reference = ''
        self.sample_list = []
        self.info_list = []
        self.format_list = []
        self.alt_list = []
        self.add_format('GT', 1, 'String', 'Genotype')
    
    def add_header(self, header):
        for line in header:
            if line.split('=')[0] == '##fileformat':
                self.file_format = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##reference':
                self.reference = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##INFO':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_info(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##ALT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_alt(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##FORMAT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_format(*[b.split('=')[1] for b in r.findall(a)])
            elif line[0] == '#' and line[1] != '#':
                self.sample_list = line.rstrip().split('\t')[9:]
    # return the VCF header
    def get_header(self):
        header = '\n'.join(['##fileformat=' + self.file_format,
                            '##fileDate=' + time.strftime('%Y%m%d'),
                            '##reference=' + self.reference] + \
                           [i.hstring for i in self.info_list] + \
                           [a.hstring for a in self.alt_list] + \
                           [f.hstring for f in self.format_list] + \
                           ['\t'.join([
                               '#CHROM',
                               'POS',
                               'ID',
                               'REF',
                               'ALT',
                               'QUAL',
                               'FILTER',
                               'INFO',
                               'FORMAT'] + \
                                      self.sample_list
                                  )])
        return header
    
    def add_info(self, id, number, type, desc):
        if id not in [i.id for i in self.info_list]:
            inf = Info(id, number, type, desc)
            self.info_list.append(inf)
    
    def add_alt(self, id, desc):
        if id not in [a.id for a in self.alt_list]:
            alt = Alt(id, desc)
            self.alt_list.append(alt)
    
    def add_format(self, id, number, type, desc):
        if id not in [f.id for f in self.format_list]:
            fmt = Format(id, number, type, desc)
            self.format_list.append(fmt)
    
    def add_sample(self, name):
        self.sample_list.append(name)

class Info(object):
    def __init__(self, id, number, type, desc):
        self.id = str(id)
        self.number = str(number)
        self.type = str(type)
        self.desc = str(desc)
        # strip the double quotes around the string if present
        if self.desc.startswith('"') and self.desc.endswith('"'):
            self.desc = self.desc[1:-1]
        self.hstring = '##INFO=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

class Alt(object):
    def __init__(self, id, desc):
        self.id = str(id)
        self.desc = str(desc)
        # strip the double quotes around the string if present
        if self.desc.startswith('"') and self.desc.endswith('"'):
            self.desc = self.desc[1:-1]
        self.hstring = '##ALT=<ID=' + self.id + ',Description=\"' + self.desc + '\">'

class Format(object):
    def __init__(self, id, number, type, desc):
        self.id = str(id)
        self.number = str(number)
        self.type = str(type)
        self.desc = str(desc)
        # strip the double quotes around the string if present
        if self.desc.startswith('"') and self.desc.endswith('"'):
            self.desc = self.desc[1:-1]
        self.hstring = '##FORMAT=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

class Variant(object):
    def __init__(self, var_list, vcf):
        self.chrom = var_list[0]
        self.pos = int(var_list[1])
        self.var_id = var_list[2]
        self.ref = var_list[3]
        self.alt = var_list[4]
        self.qual = var_list[5]
        self.filter = var_list[6]
        self.sample_list = vcf.sample_list
        self.info_list = vcf.info_list
        self.info = dict()
        self.format_list = vcf.format_list
        self.active_formats = list()
        self.gts = dict()
        # make a genotype for each sample at variant
        for i in xrange(len(self.sample_list)):
            s_gt = var_list[9+i].split(':')[0]
            s = self.sample_list[i]
            self.gts[s] = Genotype(self, s, s_gt)
        # import the existing fmt fields
        for i in xrange(len(self.sample_list)):
            s = self.sample_list[i]
            for j in zip(var_list[8].split(':'), var_list[9+i].split(':')):
                self.gts[s].set_format(j[0], j[1])
        
        self.info = dict()
        i_split = [a.split('=') for a in var_list[7].split(';')] # temp list of split info column
        for i in i_split:
            if len(i) == 1:
                i.append(True)
            self.info[i[0]] = i[1]
    
    def set_info(self, field, value):
        if field in [i.id for i in self.info_list]:
            self.info[field] = value
        else:
            sys.stderr.write('\nError: invalid INFO field, \"' + field + '\"\n')
            exit(1)
    
    def get_info(self, field):
        return self.info[field]
    
    def get_info_string(self):
        i_list = list()
        for info_field in self.info_list:
            if info_field.id in self.info.keys():
                if info_field.type == 'Flag':
                    i_list.append(info_field.id)
                else:
                    i_list.append('%s=%s' % (info_field.id, self.info[info_field.id]))
        return ';'.join(i_list)
    
    def get_format_string(self):
        f_list = list()
        for f in self.format_list:
            if f.id in self.active_formats:
                f_list.append(f.id)
        return ':'.join(f_list)
    
    def genotype(self, sample_name):
        if sample_name in self.sample_list:
            return self.gts[sample_name]
        else:
            sys.stderr.write('\nError: invalid sample name, \"' + sample_name + '\"\n')
    
    def get_var_string(self):
        s = '\t'.join(map(str,[
            self.chrom,
            self.pos,
            self.var_id,
            self.ref,
            self.alt,
            '%0.2f' % self.qual,
            self.filter,
            self.get_info_string(),
            self.get_format_string(),
            '\t'.join(self.genotype(s).get_gt_string() for s in self.sample_list)
        ]))
        return s

class Genotype(object):
    def __init__(self, variant, sample_name, gt):
        self.format = dict()
        self.variant = variant
        self.set_format('GT', gt)
    
    def set_format(self, field, value):
        if field in [i.id for i in self.variant.format_list]:
            self.format[field] = value
            if field not in self.variant.active_formats:
                self.variant.active_formats.append(field)
                # sort it to be in the same order as the format_list in header
                self.variant.active_formats.sort(key=lambda x: [f.id for f in self.variant.format_list].index(x))
        # else:
        #     sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
        #     exit(1)
    
    def get_format(self, field):
        return self.format[field]
    
    def get_gt_string(self):
        g_list = list()
        for f in self.variant.active_formats:
            if f in self.format:
                if type(self.format[f]) == float:
                    g_list.append('%0.2f' % self.format[f])
                else:
                    g_list.append(self.format[f])
            else:
                g_list.append('.')
        return ':'.join(map(str,g_list))
        
#class mei_variant(object):        
#    def __init__(call_list):
#        self.chrom =
#        self.pos =
#        self.mei =
#        self.evs =
        
class mei_call(object):
    #should also make a call group class to store these
    def __init__(self, meilist):
        self.chrom = meilist[0]
        self.startlist = meilist[1]
        self.endlist = meilist[2]
        self.alt = meilist[3]
        self.mei = meilist[4]
        self.evs = meilist[5]
        self.strands = get_strands(self.alt)

def get_strands(altstr):
    # 4 possible alt configurations:
    if altstr.startswith("N]"):
        ref_strand, mei_strand = "+","+"
    elif altstr.startswith("N["):
        ref_strand, mei_strand = "+","-"
    elif altstr.startswith("]"):
        ref_strand, mei_strand = "-","+"
    elif altstr.startswith("["):
        ref_strand, mei_strand = "-","-"
    return ref_strand,mei_strand


def merge_meis(var_cluster):
    for mei_group in var_cluster:
        ori_groups = defaultdict(list)
        cluster = var_cluster[mei_group]
        for var in cluster:
            ori_groups[var.alt[:2]].append(var)
        
        merged = []
        #this will iterate over groups of calls that are the same element
        #and in the same ref/mei orientation
        for ori, variants in ori_groups.viewitems():
            SU,SR,PE=0,0,0
            startlist = []
            poslist = []
            endlist = []
            for var in variants:
                SU += int(var.info['SU'])
                SR += int(var.info['SR'])
                PE += int(var.info['PE'])
                poslist.append(var.pos)
                
                #stuff to calculate CI span
                span = map(int, var.info['CIPOS95'].split(","))
                s1 = var.pos + span[0] - 1
                e1 = var.pos + span[1]
                startlist.append(s1)
                endlist.append(e1)
                
                sep = "]"
                if sep not in var.alt:
                    sep = "["
                if var.alt.startswith("N"):
                    new_alt = "N"+sep+mei_group+sep
                else:
                    new_alt = sep+mei_group+sep+"N"
            
            #need to make a new variant here
            chrom = variants[0].chrom
            callset = mei_call([chrom, startlist, endlist, new_alt, mei_group, {'SU':SU, 'SR':SR, 'PE':PE}])
            merged.append(callset)
       

        #five = []
        #three = []
        starts = []
        ends = []
        SU5,SU3,SR5,SR3,PE5,PE3=0,0,0,0,0,0
        
        
        for call in merged:
            starts.extend(call.startlist)
            ends.extend(call.endlist)
            if call.strands[0] == "+":
                #five.append(min(call.poslist))
                SU5+=call.evs['SU']
                SR5+=call.evs['SR']
                PE5+=call.evs['PE']
            elif call.strands[0] == "-":
                #three.append(max(call.poslist))
                SU3+=call.evs['SU']
                SR3+=call.evs['SR']
                PE3+=call.evs['PE']
        SR = SR5+SR3    
        PE = PE5+PE3
        passed = False
        #no filtering for now.
        if SR + PE >= 2:
            passed = True
        
            
        if passed:
            bstart = min(startlist)
            bend = max(endlist)
            infostr = "SU5="+str(SU5)+";SU3="+str(SU3)+";SR5="+str(SR5)+";SR3="+str(SR3)+";PE5="+str(PE5)+";PE3="+str(PE3)
            return [merged[0].chrom, str(bstart), str(bend), mei_group.split("_")[1], infostr]
        else:
            return False
            
                
        

# primary function
def vcfToBedpe(vcf_file, bedpe_out, mei_prefix="moblist", window=1):
    vcf = Vcf()
    in_header = True
    header = []
    prev_var = False
    var_cluster = defaultdict(list)
    
    for line in vcf_file:
        if in_header:
            if line[0] == '#':
                header.append(line)
                if line[1] != '#':
                    sample_list = line.rstrip().split('\t')[9:]
                continue
            else:
                bedpe_out.write("\t".join(["#CHROM","START","END","MEI","INFO"])+"\n")
                in_header = False
                vcf.add_header(header)
        
        v = line.rstrip().split('\t')
        var = Variant(v, vcf)
        
        if var.info['SVTYPE'] != 'BND':
            continue
        
        else:
            if mei_prefix in var.chrom:
                continue
            
            ref_break = var.pos
            sep = '['
            if sep not in var.alt:
                sep = ']'
            r = re.compile(r'\%s(.+?)\%s' % (sep, sep))
            mei_chrom, mei_break = r.findall(var.alt)[0].split(':')
            mei_break = int(mei_break)
            mei_group = mei_chrom.split(".")[0]
            
            if prev_var:
                #catch (some) unsorted VCF files
                if (var.chrom == prev_var.chrom) and (var.pos - prev_var.pos < 0):
                    exit("Error: Input .vcf must be coordinate sorted!")
                elif (var.chrom != prev_var.chrom) or (var.pos - prev_var.pos > window):
                    merged = merge_meis(var_cluster)
                    if merged:
                        bedpe_out.write("\t".join(merged)+"\n")
                    var_cluster.clear()
            
            prev_var = var
            var_cluster[mei_group].append(var)
    
    # close the files
    bedpe_out.close()
    vcf_file.close()
    
    return

# --------------------------------------
# wrapper function

def main():
    # parse the command line args
    args = get_args()
    
    # call primary function
    vcfToBedpe(args.input, args.output)

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
