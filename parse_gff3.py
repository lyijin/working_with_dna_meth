#!/usr/bin/env python3

"""
> parse_gff3.py <

Python script intends to be a helper function that can be called by other
scripts to handle gff annotations.

Script requires a _file object_ (not filename).

The output is in the form of:
    dict[seqid] = {five_prime_UTR: {ID001: (start_coord1, end_coord1), ...}
                     three_prime_UTR: {ID001: (start_coord2, end_coord2), ...}
                     gene: ...
                     exon: ...
                     CDS: ...
                     ...: ...}
"""

# from http://genome.ucsc.edu/FAQ/FAQformat.html#format3:
# Here is a brief description of the GFF fields:
#   seqid - The name of the sequence. Must be a chromosome or scaffold. 
#   source - The program that generated this feature. 
#   feature - The name of this type of feature. Some examples of standard
#             feature types are "CDS", "start_codon", "stop_codon", and "exon". 
#   start - The starting position of the feature in the sequence. The first base
#           is numbered 1. 
#   end - The ending position of the feature (inclusive). 
#   score - A score between 0 and 1000. If the track line useScore attribute is
#           set to 1 for this annotation data set, the score value will
#           determine the level of gray in which this feature is displayed
#           (higher numbers = darker gray). If there is no score value, enter ".". 
#   strand - Valid entries include '+', '-', or '.' (for don't know/don't care). 
#   frame - If the feature is a coding exon, frame should be a number between
#           0-2 that represents the reading frame of the first base. If the
#           feature is not a coding exon, the value should be '.'. 
#   group - All lines with the same group are linked together into a single item.

# sample line:
#   RNA-1;1.gff:1000131     maker   five_prime_UTR  1968    1999    .   
#   +       .       ID=maker-1000131-exonerate...

# IMPORTANT NOTE ABOUT startpos AND endpos: the biological definition of
# position (1-based) is BAD at dealing with the edge case of startpos == endpos.
# If startpos == endpos, then there's no way to tell whether it's on the 
# '+' strand or the '-' strand.
# Thus, startpos and endpos produced by the parse_gff3 def follows the Python
# convention (starting base is sequence[0:1]), which allows for the
# discrimination of '+' and '-' strands.
#   e.g. if startpos = 2 and endpos = 3: '+' strand    }  both refer to
#        if endpos = 2 and startpos = 3: '-' strand    }  the same base

import csv
import re

class Gene:
    def __init__(self, coords):
        self.coords = coords
        self.mRNAs = {}
    
    def __len__(self):
        return abs(self.coords[1] - self.coords[0])
    
    def add_mRNA(self, mRNA_id, coords):
        """'mRNA' is a class!"""
        self.mRNAs[mRNA_id] = mRNA(coords)

class mRNA:
    def __init__(self, coords):
        self.coords = coords
        self.details = {}
    
    def __len__(self):
        return abs(self.coords[1] - self.coords[0])
    
    def add_feature(self, feature_id, coords):
        """coords is a tuple!"""
        if feature_id not in self.details:
            self.details[feature_id] = []
        
        self.details[feature_id].append(coords)

def natural_sort(input_list):
    tryint = lambda x: int(x) if x.isdigit() else x
    chunked_text = lambda x: [tryint(y) for y in re.split('([0-9]+)', x)]
    sorted_list = sorted(input_list, key=chunked_text)

    return sorted_list

def calc_total_exonic_length(list_of_exon_coords):
    """
    Given a bunch of exonic coordinates, calculates the total length.
    """
    return sum([max(x) - min(x) for x in list_of_exon_coords])
    
    
def get_attribute(gff3_row, attribute):
    attr_search = re.search('{}=(.*?)(;|$)'.format(attribute), gff3_row[8])
    if attr_search:
        return attr_search.group(1)
    else:
        raise AttributeError("'{}' does not contain '{}='.".format(
            '\t'.join(gff3_row), attribute))

def get_coords(gff3_row):
    """
    Returned coordinates are in the tuple (startpos, endpos).
    startpos/endpos are 0-based, not 1-based.
    """
    
    strand = gff3_row[6]
    if strand == '+':
        startpos = int(gff3_row[3]) - 1
        endpos = int(gff3_row[4])
    elif strand == '-':
        startpos = int(gff3_row[4])
        endpos = int(gff3_row[3]) - 1
    else:
        print ('Error: "{}" has no valid strand information.'\
            .format('\t'.join(gff3_row)))
        raise SystemExit()
        
    return (startpos, endpos)

def pick_longest_mRNA(gff_details):
    """
    Given mRNAs contained in Genes, retain only the longest mRNA per gene.
    
    The longest mRNA is the mRNA containing the longest total exonic length,
    not raw mRNA length.
    """
    for s in gff_details:
        for g in gff_details[s]:
            # for the sort, keys are sorted by longest exonic length in reverse.
            # in case of ties, the longest mRNA (by length) is picked.
            sorted_mRNA = sorted(gff_details[s][g].mRNAs, 
                    key=lambda x: (-calc_total_exonic_length(
                            gff_details[s][g].mRNAs[x].details['exon']), 
                    -len(gff_details[s][g].mRNAs[x])))
            
            # eliminate all other keys except for the longest mRNA.
            gff_details[s][g].mRNAs = {sorted_mRNA[0]: 
                                       gff_details[s][g].mRNAs[sorted_mRNA[0]]}
            
    return gff_details
    
def sort_features(gff_details):
    """
    Makes sure that feature coordinates are sorted in the right order
    (i.e. ascending order for genes in the 5'-3' order; descending for
    genes that are in the 3'-5' order).
    """
    for s in gff_details:
        for g in gff_details[s]:
            for m in gff_details[s][g].mRNAs:
                for f in gff_details[s][g].mRNAs[m].details:
                    # gff_details[s][g].mRNAs[m].details[f] has data in the
                    # form of coordinates, e.g. [(23559, 22882), (22781, 22387)]
                    coords = gff_details[s][g].mRNAs[m].details[f]
                    
                    if len(coords) == 1: continue
                    
                    desc_order = coords[0][0] > coords[0][1]
                    gff_details[s][g].mRNAs[m].details[f] = \
                        sorted(coords, reverse=desc_order)
    
    return gff_details
    
def parse_gff3(gff_file, select_feature='all'):
    """
    'gff_file' refers to file object containing gff file.
    'select_feature' can be used to select for one or more features of interest
    in the gff file (e.g. "three_prime_UTR", ['three_prime_UTR', 'five_prime_UTR'])
    
    NOTE: all variables are named according to convention!
    """
    gff_details = {}
    
    tsv_reader = csv.reader(gff_file, delimiter='\t')
    for row in tsv_reader:
        # ignore blank lines and comments (lines starting with '#')
        if not row: continue
        if row[0].startswith('#'): continue
        
        feature = row[2]
        # ignore lines that do not correspond to the feature wanted (with the
        # exception of gene and mRNA, we always want that information).
        if select_feature == 'all' or feature in ['gene', 'mRNA'] \
                                   or feature in select_feature:
            seqid = row[0]
            
            if seqid not in gff_details:
                gff_details[seqid] = {}
            
            coords = get_coords(row)
            if feature == 'gene':
                gene_id = get_attribute(row, 'ID')
                gff_details[seqid][gene_id] = Gene(coords)
            elif feature == 'mRNA':
                mRNA_id = get_attribute(row, 'ID')
                gff_details[seqid][gene_id].add_mRNA(mRNA_id, coords)
            else:
                # mRNA_ids might be multiple IDs, separated with commas.
                mRNA_ids = get_attribute(row, 'Parent')
                for m in mRNA_ids.split(','):
                    gff_details[seqid][gene_id].mRNAs[m].add_feature(feature, coords)
    
    return gff_details
    
# debug
# if __name__ == '__main__':
    # temp = parse_gff3(open('test.gff3'), select_feature='exon')
    # print (temp['ctg123']['gene00001'].mRNAs)
    # print (pick_longest_mRNA(temp)['ctg123']['gene00001'].mRNAs)
    