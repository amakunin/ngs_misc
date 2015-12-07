#! /usr/bin/env python
###############################################################################
# based on https://github.com/Niknafs/NGSTools/blob/master/baseParser.py
# Returns substition counts for bam file relative to fasta reference    
# Usage:
# bam_rsub.py ref.fa sample.bam
# Warning: pysam mpileup produces coverage reduced by 1 compared to samtools    
# mpileup in some regions.
###############################################################################

import argparse
import pysam
import os
import sys
from collections import Counter, OrderedDict

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Returns substition counts for bam file relative to fasta reference. 
                    Requires: pysam package
                    """
                    )
    parser.add_argument("fasta_file", help="reference fasta file")
    parser.add_argument("bam_file", help="input bam file")
    parser.add_argument("-l", "--max_length", type=int, help="maximum read length filter")
    parser.add_argument("-e", "--max_edit_dist", type=int, help="maximum edit distance to reference for the read")
    return parser.parse_args()

class parseString(object):
    
    def __init__(self, ref, string):
        self.ref = ref.upper()
        self.string = string.upper()
        self.types = Counter()  
        self.process()
        
    def process(self):
        # remove end of read character
        self.string = self.string.replace('$','')
        while self.string != '':
            if self.string[0] == '^':
                # skip two characters when encountering '^' as it indicates
                # a read start mark and the read mapping quality
                self.string = self.string[2:]
            elif self.string[0] == '*':
                # unknown?
                self.types['un'] += 1
                # skip to next character
                self.string = self.string[1:]
            
            elif self.string[0] in list('.,'):
                if (len(self.string)== 1) or (self.string[1] not in ['+','-']):
                    # a reference base
                    self.types[self.ref+self.ref] += 1
                    self.string = self.string[1:]
                elif self.string[1] == '+':
                    insertionLength = int(self.string[2])
                    #insertionSeq = self.string[3:3+ insertionLength]
                    self.types['ins'] += 1 # .append(insertionSeq)
                    self.string = self.string[3+insertionLength:]
                elif self.string[1] == '-':
                    deletionLength = int(self.string[2])
                    #deletionSeq = self.string[3:3+deletionLength]
                    self.types['del'] += 1 #.append(deletionSeq)
                    self.string = self.string[3+deletionLength:]
                    
            elif self.string[0] in list('ACTG'):
                # one of the four bases
                self.types[self.ref + self.string[0]] += 1
                self.string = self.string[1:]
            else:
                # unrecognized character
                # or a read that reports a substitition followed by an insertion/deletion
                print self.string
                self.types['un'] += 1
                self.string = self.string[1:]
        return
    def __repr__(self):
        return self.types

def main():

    args = parse_command_line_arguments()

    # length filter
    samfile = pysam.AlignmentFile(args.bam_file, "rb")
    if args.max_length:
        with  pysam.AlignmentFile('tmp.bam', "wb", template=samfile) as tmpfile:
            for read in samfile.fetch(until_eof=True):
                if read.query_length <= args.max_length:
                    tmpfile.write(read)
        samfile.close()
        samfile = pysam.AlignmentFile('tmp.bam', "rb")
        os.remove('tmp.bam')

    # edit distance filter
    if args.max_edit_dist:
        with  pysam.AlignmentFile('tmp.bam', "wb", template=samfile) as tmpfile:
            for read in samfile.fetch(until_eof=True):
                edit_dist = read.get_tag('NM')
                if edit_dist <= args.max_edit_dist:
                    tmpfile.write(read)
        samfile.close()
        samfile = pysam.AlignmentFile('tmp.bam', "rb")
        os.remove('tmp.bam')

    # pileup generation
    with pysam.AlignmentFile('tmp.bam', "wb", template=samfile) as tmpfile:
        for read in samfile.fetch(until_eof=True):
            tmpfile.write(read)
        samfile.close()

    pileup = pysam.mpileup('-f',args.fasta_file,'tmp.bam')
    os.remove('tmp.bam')

    # pileup parsing 
    subst = {'AA':0,'TT':0,'CC':0,'GG':0,
             'AC':0,'AT':0,'AG':0,
             'CA':0,'CT':0,'CG':0,
             'TA':0,'TC':0,'TG':0,
             'GA':0,'GC':0,'GT':0,
             'ins':0,'del':0,'un':0}
    for line in pileup:
        toks = line.strip('\n').split('\t')
        ref = toks[2].upper()
        alt = parseString(ref, toks[4]).__repr__()
        for alt_type,count in alt.iteritems():
            if ref in list('ACTG'):
                try:
                    subst[alt_type] += count
                except:
                    print alt, toks[1]
                    sys.exit()
            else:
                subst['un'] += count

    # output
    keyorder = "AA\tTT\tCC\tGG\tAC\tAT\tAG\tCA\tCT\tCG\tTA\tTC\tTG\tGA\tGC\tGT\tins\tdel\tun"
    print  keyorder
    subst = OrderedDict(sorted(subst.items(), key=lambda i: keyorder.index(i[0])))
    print  '\t'.join([str(x) for x in subst.values()])
            

if __name__ == '__main__':
    main()
