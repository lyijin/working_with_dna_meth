#!/usr/bin/env python3

"""
> natural_sort.py <

Does it says on the tin: can be imported into other scripts for the 
natural_sort function, or can be called from the command line to sort an input
list.

Natural sort: 15, 9, 101 --> 9, 15, 101
UNIX sort:    15, 9, 101 --> 101, 15, 9
"""
import argparse
import re
import sys

def natural_sort(input_list, reversed=False):
    tryint = lambda x: int(x) if x.isdigit() else x
    chunked_text = lambda x: [tryint(y) for y in re.split('([0-9]+)', x)]
    sorted_list = sorted(input_list, key=chunked_text, reverse=reversed)

    return sorted_list

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Script takes in an input file, natural sorts it, and spits it back out.""")
    parser.add_argument('unsorted', metavar='text_file',
                        type=argparse.FileType('r'), nargs='?',
                        default=sys.stdin, help='unsorted file.')
    parser.add_argument('-r', '--reverse', action='store_true', default=False,
                        help='reverse the sort order.')
    
    args = parser.parse_args()

    # read all lines (separators included, thus no need for the csv module) 
    # into a list:
    all_lines = []
    with args.unsorted as f:
        for line in f:
            all_lines.append(line.strip())
            
    all_lines = natural_sort(all_lines, args.reverse)
    
    for a in all_lines:
        print (a)
