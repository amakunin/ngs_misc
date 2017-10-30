#!/usr/bin/env python

import os
import sys
import pysam
import tempfile
import argparse
import subprocess


def parse_command_line_arguments():

    parser = argparse.ArgumentParser(description="""Remove contamination
                    reads from bam file
                    """)
    parser.add_argument("in_file",
                        help="input sam/bam alignment to be filtered")
    parser.add_argument("target_ref",
                        help="target reference genome fasta or bowtie2 prefix")
    parser.add_argument("contam_ref",
                        help="contamination reference genome fasta or prefix")
    parser.add_argument("-a", "--algorithm", default="bwa",
                        help="alignment algorithm: 'bwa' (bwa mem, Default) \
                            'bt2' (bowtie2 --very-sensitive), or 'bt2l' \
                            (bowtie2 --very-sensitive-local)")
    parser.add_argument("-o", "--output", help="output bam prefix.\
                            Default: <input_prefix>.decontam")
    parser.add_argument("-q", '--min_quality', type=int, default=0,
                        help="Minimum alignment quality for filtered file. \
                            Default: 0.")
    parser.add_argument("-e", '--max_edit_distance', type=int, default=99,
                        help="Maximum alignment edit distance for filtered \
                            file. Default: 99.")
    return parser.parse_args()


def exec_command(command, out_file='std', verbose=True):
    """Execute command using list of arguments and raise RuntimeError"""
    if verbose:
        sys.stderr.write(command)
        sys.stderr.write('\n' if out_file == 'std' else ' > %s\n' % out_file)
    p = subprocess.Popen(command.split(),
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (out, err) = p.communicate()
    if p.returncode != 0:
        raise RuntimeError("%s error:\n%s" % (args[0], err))
    if out_file != 'std':
        with open(out_file, 'w') as o:
            o.write(out)


def bam_sort(in_file, out_file, name_sort=False):
    """Sort alignment by coordinate (default) or by read name"""
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


def index(fasta, algorithm):
    """Generate reference index if it does not exist"""
    if algorithm == 'bwa':
        if not os.path.isfile(fasta + '.bwt'):
            exec_command('bwa index ' + fasta)
    elif algorithm == 'bt2' or algorithm == 'bt2l':
        if not os.path.isfile(fasta + '.1.bt2'):
            exec_command('bowtie2-build %s %s' % (fasta, fasta))


def align(fastq, ref, algorithm, aln):
    """Align fastq to reference using specified algorithm"""
    if algorithm == 'bwa':
        command = 'bwa mem %s %s' % (ref, fastq)
        exec_command(command, out_file=aln)
    else:
        profile = '--very-sensitive-local' if algorithm == 'bt2l' else \
                '--very-sensitive'
        command = 'bowtie2 %s -x %s -U %s -S %s' % (profile, ref,
                                                    fastq, aln)
        exec_command(command)


args = parse_command_line_arguments()
if args.output:
    base = args.output
    output = base + '.bam'
else:
    assert args.in_file.endswith('.bam')
    base = args.in_file[:-4]
    output = base + '.decontam.bam'
    base = output[:-13]

fnames = {
         'reads': base + '.fastq',
         'raw_target': base + '.tgt.sam',
         'raw_contam': base + '.cnt.sam',
         'target': base + '.tgt.bam',
         'contam': base + '.cnt.bam',
         'filter': base + '.flt.bam',
         }

# bam to fastq
btf = 'bedtools bamtofastq -i %s -fq %s' % (args.in_file, fnames['reads'])
exec_command(btf)
in_reads = sum(1 for line in open(fnames['reads']))/4

# index and align
for (aln, ref) in ((fnames['raw_target'], args.target_ref),
                   (fnames['raw_contam'], args.contam_ref)):
    index(ref, args.algorithm)
    align(fnames['reads'], ref, args.algorithm, aln)

# sort alignments by name
for (in_aln, out_aln) in ((fnames['raw_target'], fnames['target']),
                          (fnames['raw_contam'], fnames['contam'])):
    bam_sort(in_aln, out_aln, name_sort=True)

# filter to the intermediate output
with pysam.AlignmentFile(fnames['target']) as tfile, \
        pysam.AlignmentFile(fnames['contam']) as cfile:
    filter_file = pysam.AlignmentFile(fnames['filter'], 'wb', template=tfile)
    (i, j, k) = (0, 0, 0)
    for tread in tfile:
        i += 1
        cread = cfile.next()
        assert tread.query_name == cread.query_name
        # low-quality
        if tread.mapping_quality < args.min_quality:
            j += 1
        elif not tread.has_tag('NM'):
            j += 1
        elif tread.opt('NM') > args.max_edit_distance:
            j += 1
        # higher similarity to contamination
        elif cread.has_tag('NM') and cread.opt('NM') > tread.opt('NM'):
            k += 1
        else:
            filter_file.write(tread)

    filter_file.close()

# sort output by coordinate
bam_sort(fnames['filter'], output)
sys.stderr.write("%d alignments in the input file %s.\n"
                 "%d alignments aligned to reference %s.\n"
                 "%d alignments removed not passing quality filters.\n"
                 "%d alignments removed as belonging to %s.\n"
                 "%d alignments remained after filtering in %s.\n" %
                 (in_reads, args.in_file, i,  args.target_ref,
                  j, k, args.contam_ref, i-j-k, output))

# clean-up
for f in fnames:
    os.unlink(fnames[f])
