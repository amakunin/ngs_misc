#! /usr/bin/env python
###############################################################################
# Returns taxonomy counts (one per query) for tabular blast file. 
# Subject taxonomy data is taken from the last column of the first blast hit per query. 
# Usage:
# blast_tax_count.py blast_out.tabular
###############################################################################

import os
import sys
import operator

tax_count = dict()

with open(sys.argv[1]) as hit_file:
    prev_query = '0'
    for line in hit_file:
        ll = line.split('\t')
        curr_query = ll[0]
        if prev_query != curr_query:
            tax_id = ll[-1].strip('\n')
            if tax_id in tax_count.keys():
                tax_count[tax_id] += 1
            else:
                tax_count[tax_id] = 1
        prev_query = curr_query

sorted_tax_count = sorted(tax_count.items(), key=operator.itemgetter(1), reverse=True)

for tax in sorted_tax_count:
    print '\t'.join([tax[0],str(tax[1])])
