#!/usr/bin/env python3

"""
> merge_bismark_cov.py <

Given vanilla/augmented bismark cov files, combine the meth and unmeth reads
together, recompute the meth %, and leave everything else unchanged.
"""
import argparse
import csv
import sys

import natural_sort

parser = argparse.ArgumentParser(description="""
Given vanilla/augmented bismark cov files, combine the meth and unmeth reads
together, recompute the meth %, and leave everything else unchanged.""")

parser.add_argument('cov_files', metavar="cov_files",
                    type=argparse.FileType('r'), nargs='+', 
                    help="List of cov files from Bismark.")
parser.add_argument('-v', action='store_true',
                    help="verbose mode, prints progress to stderr.")
args = parser.parse_args()

# combined_data is in the format
#   combined_data[scaffold][1-based_coord] = [meth, unmeth]
combined_data = {}
counter_rows = 0
for c in args.cov_files:
    tsv_reader = csv.reader(c, delimiter='\t')
    for row in tsv_reader:
        if not row: continue
        
        if args.v:
            counter_rows += 1
            if counter_rows % 1000000 == 0:
                print ('{} rows processed...'.format(counter_rows), 
                       file=sys.stderr)
        
        scaf = row[0]
        if scaf not in combined_data: combined_data[scaf] = {}
        
        pos = int(row[1])
        if pos not in combined_data[scaf]:
            combined_data[scaf][pos] = [0, 0]
        
        meth = int(row[4])
        unmeth = int(row[5])
        combined_data[scaf][pos] = [combined_data[scaf][pos][0] + meth,
                                    combined_data[scaf][pos][1] + unmeth]

# printing
for scaf in natural_sort.natural_sort(combined_data):
    for pos in sorted(combined_data[scaf]):
        # recalculate meth % from the new meth + unmeth numbers
        meth, unmeth = combined_data[scaf][pos]
        meth_pct = round(meth / (meth + unmeth) * 100, 4)
        
        print (scaf, pos, pos, meth_pct, meth, unmeth, sep='\t')
