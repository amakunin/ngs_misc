#!/usr/bin/env python
# split a bam file by read group ID
# modified from script by Sean Davis <seandavi@gmail.com> (March 10, 2012)
# to comply with new pysam datatypes

import pysam
import argparse
import logging

logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser()
parser.add_argument('filename',help="The bam filename")

opts = parser.parse_args()

infile = pysam.Samfile(opts.filename,'rb')

header = infile.header.to_dict()

readgroups = header['RG']
# remove readgroups from header
del(header['RG'])
outfiles = {}
for rg in readgroups:
    tmphead = header
    tmphead['RG']=[rg]
    logging.info('Creating new BAM file: %s', (rg['ID']+'.bam', 'wb'))
    outfiles[rg['ID']] = pysam.Samfile(rg['ID']+'.bam', 'wb', header=tmphead)

j=0
for read in infile:
    j+=1
    idtag = [x[1] for x in read.tags if x[0]=='RG'][0]
    if((j % 100000)==0):
        logging.info('read and wrote %d records', (j))
    outfiles[idtag].write(read)

for f in outfiles.values():
    f.close()

infile.close()
