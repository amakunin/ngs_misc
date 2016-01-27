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
                    Or filters bam file by read length and edit distance.
                    Requires: pysam package
                    """
                    )
    parser.add_argument("bam_file", help="input bam file")
    parser.add_argument("-l", "--max_length", type=int, help="maximum read length filter")
    action1 = parser.add_mutually_exclusive_group(required=False)
    action1.add_argument("-e", "--max_edit_dist", type=int, help="maximum edit distance to reference for the read")
    action1.add_argument("-d", "--max_derived", type=int, help="maximum percent of derived positions in read calculated from edit distance")
    #parser.add_argument("-r", "--remove_trimming", action='store_true', help="discard all trimming data")
    action2 = parser.add_mutually_exclusive_group(required=True)
    action2.add_argument("-f", "--filter_only", action='store_true', help="only filter bam")
    action2.add_argument("-r", "--ref_fasta", help="reference fasta file")
    

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
                #sys.stdout.write(self.string[0])
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
                self.types['un'] += 1
                #sys.stdout.write(self.string[0])
                self.string = self.string[1:]
        return
    def __repr__(self):
        return self.types

def main():

    args = parse_command_line_arguments()
     
    # defaults and naming 
    if not args.max_length:
        args.max_length = 1000
    if not args.max_edit_dist and not args.max_derived:
        args.max_edit_dist = 1000
    if args.max_edit_dist:
        outbam = '%slen%dNM%d.bam' % (args.bam_file[:-3], args.max_length, args.max_edit_dist)
    elif args.max_derived:
        outbam = '%slen%dder%d.bam' % (args.bam_file[:-3], args.max_length, args.max_derived)

    with pysam.AlignmentFile(args.bam_file, "rb") as samfile:
        print '%d reads before filtering' % samfile.count()
        i = 0
        with  pysam.AlignmentFile(outbam, "wb", template=samfile) as tmpfile:
            for read in samfile.fetch():
                edit_dist = read.get_tag('NM')
                read_length = read.query_length
                if read_length <= args.max_length:
                    if args.max_edit_dist:
                        if edit_dist <= args.max_edit_dist:
                            tmpfile.write(read)
                            i += 1
                    elif args.max_derived:
                        if 100*edit_dist/read_length < args.max_derived:
                            tmpfile.write(read)
                            i += 1
        print '%d reads after filtering' % i

    # pileup generation
    if not args.filter_only:
        pileup = pysam.mpileup('-f',args.ref_fasta,outbam)
        os.remove(outbam)

        # pileup parsing 
        subst = {'AA':0,'TT':0,'CC':0,'GG':0,
                 'AC':0,'AT':0,'AG':0,
                 'CA':0,'CT':0,'CG':0,
                 'TA':0,'TC':0,'TG':0,
                 'GA':0,'GC':0,'GT':0,
                 'ins':0,'del':0,'un':0,'un_ref':0}
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
                    subst['un_ref'] += count

        # output
        keyorder = "AA\tTT\tCC\tGG\tAC\tAT\tAG\tCA\tCT\tCG\tTA\tTC\tTG\tGA\tGC\tGT\tins\tdel\tun\tno_ref"
        print  keyorder
        subst = OrderedDict(sorted(subst.items(), key=lambda i: keyorder.index(i[0])))
        print  '\t'.join([str(x) for x in subst.values()])
            

if __name__ == '__main__':
    main()
