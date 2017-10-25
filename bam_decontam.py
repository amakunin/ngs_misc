#!/usr/bin/env python

import os
import sys
import pysam
import tempfile
import argparse
import subprocess


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description="""
                    Removes contamination reads from bam file \
                    using bowtie2 alignment to alternate reference
                    Requirements: pysam python package; bowtie2 \
                    and bedtools in $PATH

                    Workflow:
                    1) bedtools bamtofastq (assume unpaired reads)
                    2) bowtie2 align to  contamination reference
                    3) compare edit distances to target and contamination \
                    references
                    """)
    parser.add_argument("in_file",
                        help="input sam/bam alignment to be filtered")
    # parser.add_argument("target_ref",
    #                   help="target reference genome fasta or bowtie2 prefix")
    parser.add_argument("contam_ref",
                        help="contamination reference genome fasta or bowtie2 \
                        prefix")
    parser.add_argument("-b", default="--very-sensitive-local -p 1",
                        help="bowtie2 performance adjustment arguments")
    parser.add_argument("-o",
                        help="output bam name.\
                         Default: <input_prefix>.decontam.bam")
    parser.add_argument("-q", '--min_quality', default=0,
                        help="Minimum alignment quality for filtered file. \
                        Default - 0.")
    parser.add_argument("-e", '--max_edit_distance', default=99,
                        help="Maximum alignment edit distance for filtered \
                        file. Default - 99.")
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
        raise Exception('Invalid samtools verion.')


def exec_command(args):
    # sys.stderr.write(' '.join(args) + '\n')
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = p.communicate()
    if p.returncode != 0:
        raise RuntimeError("%s error:\n%s" % (args[0], err))


def index(fasta):
    try:
        exec_command('bowtie2-inspect', fasta)
    except RuntimeError:
        exec_command('bowtie2-build', fasta, fasta)


args = parse_command_line_arguments()
filter_bam = args.o if args.o else args.in_file[:-4] + '.decontam.bam'
min_qual = int(args.min_quality)
max_edit = int(args.max_edit_distance)

# bam to fastq
fq = tempfile.NamedTemporaryFile(suffix='.fastq', delete=False)
fq.close()
btf_command = ['bedtools', 'bamtofastq', '-i',
               args.in_file, '-fq', fq.name]
exec_command(btf_command)

# bowtie2 alignment
r_aln = tempfile.NamedTemporaryFile(suffix='.sam', delete=False)
r_aln.close()
bt2_aln = 'bowtie2 -x %s -U %s -S %s' % (args.contam_ref, fq.name, r_aln.name)
exec_command(bt2_aln.split() + args.b.split())

# name-sort
t_aln = tempfile.NamedTemporaryFile(suffix='_t.bam', delete=False)
c_aln = tempfile.NamedTemporaryFile(suffix='_c.bam', delete=False)
t_aln.close()
c_aln.close()
bam_sort(args.in_file, t_aln.name, name_sort=True)
bam_sort(r_aln.name, c_aln.name, name_sort=True)

# filter to the intermediate output
f_aln = tempfile.NamedTemporaryFile(suffix='_f.bam', delete=False)
f_aln.close()
with pysam.AlignmentFile(t_aln.name) as tfile, \
        pysam.AlignmentFile(c_aln.name) as cfile:
    filter_file = pysam.AlignmentFile(f_aln.name, 'wb', template=tfile)
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
bam_sort(f_aln.name, filter_bam)
sys.stderr.write("%d alignments in %s.\n"
                 "%d alignments removed as belonging to %s.\n"
                 "%d alignments remained after filtering in %s.\n" %
                 (i, args.in_file, i-j, args.contam_ref, j, filter_bam))

# clean-up
os.unlink(fq.name)
os.unlink(r_aln.name)
os.unlink(t_aln.name)
os.unlink(c_aln.name)
os.unlink(f_aln.name)
