#!/usr/bin/env python3

"""
> tabulate_tsvs.py <

Script to merge similar tsvs together, by checking common keys in the first
column (by default, can be toggled with --key).

Headers are assumed to be NOT PRESENT by default!
"""
import argparse
import csv
import re
import sys

import natural_sort

parser = argparse.ArgumentParser(description="""
Script to merge similar tsvs together, by checking common keys in the first
column.
""")

parser.add_argument('tsv_files', metavar="tsv_filenames",
                    type=argparse.FileType('r'), nargs='+',
                    help="tab-separated filenames.")
parser.add_argument('--header', action='store_true',
                    help="tsv files have a header line (default: no headers!).")
parser.add_argument('--key', '-k', metavar="key_columns",
                    type=int, nargs='+', default=[0],
                    help="(0-based) columns as keys (default: 0).")
parser.add_argument('--col', '-c', metavar="retained_columns",
                    type=int, nargs='+',
                    help="(0-based) columns as values (default: all except -k)")
parser.add_argument('-v', action='store_true',
                    help="verbose mode, prints extra details to stderr.")
args = parser.parse_args()

giant_dict = {}
max_cols = len(args.col) if args.col else 0

if args.v:
    print ('Files used: {}'.format(', '.join([x.name for x in args.tsv_files])),
        file=sys.stderr)
        
    if args.header:
        print ('Headers are PRESENT.', file=sys.stderr)
    else:
        print ('Headers are NOT PRESENT.', file=sys.stderr)

# read data
for tsv_file in args.tsv_files:
    tsv_reader = csv.reader(tsv_file, delimiter='\t')
    
    # skip header row if args.header
    if args.header:
        header = next(tsv_reader)
    
    if args.v:
        print ('\rReading file #{}'.format(args.tsv_files.index(tsv_file) + 1),
            end='', file=sys.stderr)
    
    for row in tsv_reader:
        # skip empty rows
        if not row: continue
        
        row_key = '\t'.join([row[k] for k in args.key])
        
        if args.col:
            row_val = '\t'.join([row[c] for c in args.col])
        else:
            row_val = '\t'.join([row[c] for c in range(len(row))
                                 if c not in args.key])
            if len(row) > max_cols:
                max_cols = len(row_val.split('\t'))
        
        if row_key not in giant_dict:
            giant_dict[row_key] = {}

        giant_dict[row_key][tsv_file.name] = row_val

if args.v:
    print ('\nUnion of all files produces {} rows.'.format(len(giant_dict)),
        file=sys.stderr)

# write data
# header lines
print ('\t'.join([''] * len(args.key) +\
                 [x.name + '\t' * (max_cols - 1) for x in args.tsv_files]))
if args.header:
    header_key = '\t'.join([header[k] for k in args.key])
    if args.col:
        header_val = '\t'.join([header[c] for c in args.col])
    else:
        header_val = '\t'.join([header[c] for c in range(len(header))
                                if c not in args.key])
    print ('\t'.join([header_key] + [header_val] * len(args.tsv_files)))
    
for g in natural_sort.natural_sort(giant_dict):
    print ('\t'.join([g] + [giant_dict[g][x.name] if x.name in giant_dict[g]
                            else '\t' * (max_cols - 1) for x in args.tsv_files]))
