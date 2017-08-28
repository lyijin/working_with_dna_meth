#!/usr/bin/env python3

"""
> filter_miscalled_Cs.py <

Script to process the .cov files produced by bismark_methylation_extractor, and
removes methylated positions that potentially arise by chance.

A super-conservative estimate of 1% miscall rate is assumed (might be as low
as 0.1%, evidenced by the CHG and CHH call rate in S. pistillata libraries).

The script ONLY considers positions that have at least one methylated position.
There's no way to "miscall" a C as methylated... if it's not methylated. As
a result, the B-H correction only takes into consideration the universe of 
positions with at least one methylated read.

A binomial test is run to check the cumulative P value (P >= k), where k is the 
number of methylated sites, and P values are further B-H corrected.
"""
import argparse
import csv
import sys
import re

import scipy.stats

import correct_p_values
import natural_sort

parser = argparse.ArgumentParser(description="""
Script to process the .cov files produced by bismark_methylation_extractor, and
removes methylated positions that potentially arise by chance.""")

parser.add_argument('bismark_cov', metavar="cov_filename",
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin, help="Bismark .cov filename.")
parser.add_argument('-v', action='store_true',
                    help="verbose mode, prints progress to stderr.")
args = parser.parse_args()

# calculate binomial P values first
counter_rows = 0
p_values = {}
with args.bismark_cov as tsv_file:
    tsv_reader = csv.reader(tsv_file, delimiter='\t')
    for row in tsv_reader:
        if not row: continue 

        if args.v:
            counter_rows += 1
            if counter_rows % 100000 == 0:
                print ('{} rows processed...'.format(counter_rows), file=sys.stderr)
            
        # column 4: sum(methylated reads)
        # column 5: sum(non-methylated reads)
        meth = int(row[4])
        non_meth = int(row[5])
        
        # ignore lines with meth == 0
        if meth > 0:
            # scipy.stats.binom.sf calculates P(X > n) by default; the "-1"
            # is added to achieve P(X >= n).
            binom_p = scipy.stats.binom.sf(meth-1, meth + non_meth, 0.01)
            p_values[str(row)] = binom_p

# run B-H correction
if args.v:
    print ('correcting p values...', file=sys.stderr)

p_values = correct_p_values.correct_p_values(p_values)

# then print valid rows out
if args.v:
    print ('printing output and exiting.', file=sys.stderr)
    
for p in natural_sort.natural_sort(p_values):
    # print out rows that have corrected P values < 0.05
    if p_values[p] < 0.05:
        print ('\t'.join(eval(p)))