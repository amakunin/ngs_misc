#!/usr/bin/env python

import pysam
import sys
import os
import tempfile
import argparse
import subprocess

def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description=    
                    """
                    Removes contamination reads from bam file using bowtie2 alignment to alternate reference
                    Requirements: pysam python package, bowtie2 (in $PATH) and bedtools (in $PATH)

                    Workflow:
                    1) bedtools bamtofastq (assume unpaired reads)
                    2) bowtie2 align to alternate reference (use same alignment strategy as for input bam generation)
                    3) compare edit distances for all reads in original vs alternate alignment
                    """
                    )
    parser.add_argument("target_file",
                        help="input sam/bam alignment to be filtered")
    parser.add_argument("-b",
                        help="bowtie2 arguments including path to alternate reference index")
    parser.add_argument("-o",
                        help="output bam alignment. Default: <input_prefix>.decontam.bam")
    parser.add_argument("-q", '--min_quality', default=0,
                        help="Minimum alignment quality for filtered file. Default - 0.")
    parser.add_argument("-e", '--max_edit_distance', default=99,
                        help="Maximum alignment edit distance for filtered file. Default - 99.")
    return parser.parse_args()

def bam_sort(in_file, out_file, name_sort=False):
    samtools_release = pysam.version.__samtools_version__[0] 
    if samtools_release == '0':
        if name_sort:
            pysam.sort('-n', in_file, out_file[:-4])
        else:
            pysam.sort(in_file, out_file[:-4])
    elif samtools_release == '1':
        if name_sort:
            pysam.sort('-n', '-T', '/tmp/bam_nsort', '-o', out_file, in_file)
        else:
            pysam.sort('-T', '/tmp/bam_sort', '-o', out_file, in_file)
    else:
        raise Exception('Invalid samtools verion: %s. Valid versions 0.*.*, 1.*.*.' % (pysam.version.__samtools_version__[0]))

def exec_command(args):
    print ' '.join(args)
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
    (out, err) = p.communicate()
    if p.returncode != 0:
        raise Exception("%s error:\n%s" % (args[0], err))

args = parse_command_line_arguments()
filter_bam = args.o if args.o else args.target_file[:-4] + '.decontam.bam'
min_qual = int(args.min_quality)
max_edit = int(args.max_edit_distance)

#bam to fastq
fq = tempfile.NamedTemporaryFile(suffix = '.fastq', delete=False)
fq.close()
btf_command = ['bedtools','bamtofastq','-i',args.target_file,'-fq',fq.name]
exec_command(btf_command)

# bowtie2 alignment
alt_aln = tempfile.NamedTemporaryFile(suffix = '_alt.sam', delete=False)
alt_aln.close()
bt2_command = ['bowtie2','-U',fq.name,'-S',alt_aln.name] + args.b.split()
exec_command(bt2_command)


# name-sort
t = tempfile.NamedTemporaryFile(suffix = '_t.bam', delete=False)
c = tempfile.NamedTemporaryFile(suffix = '_c.bam', delete=False)
t.close()
c.close()
bam_sort(args.target_file, t.name, name_sort=True)
bam_sort(alt_aln.name, c.name, name_sort=True)

# filter to the intermediate output
f = tempfile.NamedTemporaryFile(suffix = '_f.bam', delete=False)
f.close()
with pysam.AlignmentFile(t.name) as tfile, pysam.AlignmentFile(c.name) as cfile:
    filter_file = pysam.AlignmentFile(f.name,'wb', template=tfile)
    #contam_file = pysam.AlignmentFile(contam_bam,'wb', template=cfile)
    #unmap_file = pysam.AlignmentFile(unmap_bam,'wb', template=tfile)
    i = 0
    j = 0
    for tread in tfile:
        i += 1
        cread = cfile.next()
        assert tread.query_name == cread.query_name
        if tread.mapping_quality >= min_qual and tread.opt('NM') <= max_edit:
            if not cread.has_tag('NM') or tread.opt('NM') <= cread.opt('NM'):
                j += 1
                filter_file.write(tread)
            
    filter_file.close()

# coordinate-sort output
bam_sort(f.name, filter_bam)
sys.stderr.write('Read %d alignments from input file %s.\nWritten %d alignments to filtered output file %s.\n'%(i, args.target_file, j, filter_bam))

#clean-up
os.unlink(fq.name)
os.unlink(alt_aln.name)
os.unlink(t.name)
os.unlink(c.name)
os.unlink(f.name)

