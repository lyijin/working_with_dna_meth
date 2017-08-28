#!/usr/bin/env python3

"""
> annotate_bismark_cov.py <

Script calculates the distribution of methylation over entire genic/intergenic
regions, and prints out:
  - Gene name, if in a genic region; "intergenic", if outside.
  - How far the methylated position is along the genic/intergenic region,
    in % terms.
  - How far the methylated position is from the 5' end of the genic/intergenic
    region. Values range from +1 (right splat on 5' edge) to len(region).
  - How far the methylated position is from the 3' end of the genic/intergenic
    region. Values range from -1 (right splat on 3' edge) to -len(region).
  - Same information for exons/introns as for genes, when present.
"""
import argparse
import csv
import re
import sys
import time

import numpy as np

import natural_sort
import parse_fasta
import parse_gff3

parser = argparse.ArgumentParser(description="""
Script calculates the distribution of methylation over entire genic/intergenic
regions, and prints out lotsa details.""")

parser.add_argument('genome_fasta', metavar='fasta_file',
                    type=argparse.FileType('r'),
                    help='genome FASTA file.')
parser.add_argument('genome_gff3', metavar='gff3_file',
                    type=argparse.FileType('r'),
                    help='genome GFF3 annotation.')
parser.add_argument('bismark_cov', metavar="cov_filename",
                    type=argparse.FileType('r'), nargs='?',
                    default=sys.stdin, help="Bismark .cov filename.")
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints diagnostic stuff to stderr.')
args = parser.parse_args()

def is_numpy_array_filled(np_array, coords):
    return np.any(np_array[min(coords):max(coords)])

def assign_value_to_numpy_array(np_array, coords, assigned_value):
    """
    Coords are in a tuple (coords1, coords2). +ve and -ve values are assigned
    respectively, depending on coords1 < coords2 or vice versa.
    """
    if coords[1] > coords[0]:
        np_array[coords[0]:coords[1]] = assigned_value
    elif coords[0] > coords[1]:
        np_array[coords[1]:coords[0]] = -assigned_value

    return np_array

def calc_coord_length(gene_coords):
    """
    Given a coord like [10,20], return the length of the sub-sequence (10).
    """
    return max(gene_coords) - min(gene_coords)

def calc_displacement(seq, coords, return_type='all'):
    """
    Displacement: given a chromosomal coordinate, calculate how far it is from
    the start of the gene.
    
              <-------|
    0 0 0 0 0 1 1 1 1 1 1 1 0 0     where 1 denotes genic regions
             | |         | |
           s_coord     e_coord
    
    "seq" can be NumPy arrays or lists.
    
    Coordinates are zero-based (start, end) lists. Make sure individual
    positions are converted to lists before calling the def!
      e.g.  4 --> [3, 4]
    Using tuples allows for checking location of k-mers in a longer n-mer.
      e.g.  "BC" in "ABCDE" --> 2nd dimer out of 4 possible dimers
    and allows for specifying - strand coordinates:
      e.g. [4, 3]
    
    Resulting value can either be "percentage" (coords / total length of gene)
    or "absolute" (coords - start of gene)
    """
    
    # sanity checking inputs
    assert len(coords) == 2, "make sure 'coords' is a list!"
    assert return_type in ['all', 'percentage', 'absolute', 'absolute_inverse'], \
        "return_type != 'all', 'percentage', 'absolute' or 'absolute_inverse'."
    
    if return_type == 'all':
        return_type = ['percentage', 'absolute', 'absolute_inverse']
    
    # go upstream (scaffold's perspective) first
    s_coord = min(coords)
    if s_coord > 0:
        while seq[s_coord-1] == seq[min(coords)]:
            s_coord -= 1
            if s_coord == 0:
                break
    
    # then downstream (scaffold's perspective)
    e_coord = max(coords)
    if e_coord < len(seq):
        while seq[e_coord] == seq[max(coords)-1]:
            e_coord += 1
            if e_coord == len(seq):
                break
    
    feature_coord = [s_coord, e_coord]
    # if original start coord > end coord, flip the coords of feature_coord,
    # so that the 3' end becomes the "head" of the feature and the 5' end
    # becomes the "tail" of the feature
    if coords[0] > coords[1]:
        feature_coord = [-feature_coord[1], -feature_coord[0]]
        coords = [-coords[0], -coords[1]]
    
    returned_list = []
    if 'percentage' in return_type:
        # takes the HALFWAY point within coords, and checks how far it is from
        # the 5' end
        #   range of values returned: (0, 1)
        loc_pct = ((coords[0] + coords[1]) / 2 - feature_coord[0]) / \
                   (calc_coord_length(feature_coord))
        returned_list.append(round(loc_pct, 5))
    if 'absolute' in return_type:
        # range of values returned: [0, len(feature_coord))
        returned_list.append(coords[0] - feature_coord[0])
    if 'absolute_inverse' in return_type:
        # range of values returned: (-len(feature_coord), -1]
        returned_list.append(coords[1] - feature_coord[1] - 1)

    return [str(x) for x in returned_list]

def translate_exon_info(ei_integer, gene_id):
    """
    Change integers stored in NumPy array into strings, e.g.:
        1 --> Exon_1
        6 --> Intron_3
        
    Also return info in the negative form (i.e. Exon_-1 is the last exon,
    Exon_-2 is penultimate exon...).
    """
    # in rare cases, the longest transcript in a gene is shorter than the 
    # gene itself. Thus, there's a possibility that
    #     gene_info != AND exon_info == 0 !!
    if not ei_integer: return 'no_info'
    
    ei_integer = abs(int(ei_integer))   # NumPy ints aren't... really ints.    
    ei = 'Exon' if ei_integer % 2 else 'Intron'
    
    # round() is dangerous. Banker's rounding.
    ei_number = (ei_integer + 1) // 2       
    ei_inverse = ei_number - exon_count[gene_id] - (ei_integer % 2)
    return '_'.join([ei, str(ei_number), str(ei_inverse)])

# read sequences
sequence_lengths = parse_fasta.get_all_sequences(args.genome_fasta, 
        'fasta', lengths_only=True)
if args.verbose:
    print ('[{}] Lengths for {} sequences parsed.'.format(
           time.asctime(), len(sequence_lengths)), file=sys.stderr)

# read coordinates of genes and exons from .gff3 file.
scaffold_gff3 = parse_gff3.parse_gff3(args.genome_gff3, 'exon')
# as genes might contain overlapping isoforms, the longest isoform is chosen,
# if multiples exist.
scaffold_gff3 = parse_gff3.pick_longest_mRNA(scaffold_gff3)
# make sure features in all mRNAs are sorted properly (for exon numbering).
scaffold_gff3 = parse_gff3.sort_features(scaffold_gff3)

# genic regions are denoted in a NumPy array as follows:
#    0: intergenic region
#  +ve: gene in the + strand
#  -ve: gene in the - strand
#
# gene_names is a dict that converts ints in gene_info back to gene names.
#
# NOTE: if genes overlap, the first one that gets annotated (sorted 
# alphabetically) has priority.
gene_info = {}
gene_names = {}
gene_counter = 0

# similar to genic data, exonic data is stored in a NumPy array. However,
# exons/introns in all genes have the following numbering scheme:
#    0: intergenic region
# 1, 2: exon 1, intron 1
# 3, 4: exon 2, intron 2, ...
#  +ve: gene in the + strand
#  -ve: gene in the - strand
#
# exon_count stores # of exons in that gene.
exon_info = {}
exon_count = {}

for s in natural_sort.natural_sort(scaffold_gff3):
    gene_info[s] = np.zeros(sequence_lengths[s], np.int32)
    exon_info[s] = np.zeros(sequence_lengths[s], np.int8)
    for gene_id in natural_sort.natural_sort(scaffold_gff3[s]):
        # gene detail annotation
        gene_coords = scaffold_gff3[s][gene_id].coords
        
        # jump out of loop if the assigned region is not zeroes ([0, 0, ...]).
        if is_numpy_array_filled(gene_info[s], gene_coords): continue
        
        gene_counter += 1
        gene_info[s] = assign_value_to_numpy_array(gene_info[s], gene_coords, 
                                                   gene_counter)
        gene_names[gene_counter] = gene_id

        # exon detail annotation
        #        
        # sole_mRNA is a mRNA object, with the assumption that *.mRNAs
        # returns a dict containing {mRNA_name: mRNA_object}
        sole_mRNA = list(scaffold_gff3[s][gene_id].mRNAs.values())[0]
        
        exon_coords = sole_mRNA.details['exon']
        exon_count[gene_id] = len(exon_coords)
        unpacked_coords = [y for x in exon_coords for y in x]
        exon_intron_coords = [(x, y) for x, y in zip(unpacked_coords[:-1],
                                                     unpacked_coords[1:])]
        
        for n, ei in enumerate(exon_intron_coords):
            exon_info[s] = assign_value_to_numpy_array(exon_info[s], ei, n+1)
        
if args.verbose:
    print ('[{}] {} sequences with gene annotations parsed.'.format(
           time.asctime(), len(scaffold_gff3)), file=sys.stderr)
    print ('[{}] {} genes parsed in total.'.format(
           time.asctime(), gene_counter), file=sys.stderr)

# read per-position # of methylated and unmethylated bases
tsv_reader = csv.reader(args.bismark_cov, delimiter='\t')
counter = 0
for row in tsv_reader:
    # skip empty rows
    if not row: continue
    
    # to make script work on annotated bismark covs (i.e. overwrite previous
    # annotation), retain only first 6 columns
    row = row[:6]
    
    # skip lines with scaffold names not in the gff3 file
    scaf = row[0]
    if scaf not in gene_info:
        print ('\t'.join(row), 'no_info', sep='\t')
        continue
    
    # all systems go!
    start = int(row[1]) - 1
    end = int(row[2])
    meth = int(row[4])
    unmeth = int(row[5])
    
    # add gene info first
    extra_info = []
    if gene_info[scaf][start] == 0:
        extra_info.append('intergenic')
        intergenic_info = calc_displacement(gene_info[scaf], [start, end])
        extra_info += intergenic_info
        
        # get info regarding regions flanking the intergenic region:
        #   scaf_start: 5' end of the intergenic region is start of scaffold
        #   scaf_end: 3' end of the intergenic region is end of scaffold
        #   upstream_watson: 5' end of intergenic region is 3' end of a gene
        #   upstream_crick: 5' end of intergenic region is 5' end of a gene
        #   downstream_watson: 3' end of intergenic region is 5' end of a gene
        #   downstream_crick: 3' end of intergenic region is 3' end of a gene
        fivep_dist = int(intergenic_info[1])
        threep_dist = abs(int(intergenic_info[2]))
        if fivep_dist == start:
            fivep_info = 'scaf_start'
        elif gene_info[scaf][start - fivep_dist - 1] > 0:
            gene_id = gene_names[gene_info[scaf][start - fivep_dist - 1]]
            fivep_info = 'upstream_watson|' + gene_id
        elif gene_info[scaf][start - fivep_dist - 1] < 0:
            gene_id = gene_names[abs(gene_info[scaf][start - fivep_dist - 1])]
            fivep_info = 'upstream_crick|' + gene_id
        else:
            # if this is triggered, serious debugging is needed
            print ("[{}] invalid fivep_info value encountered!".format(
                   time.asctime()), file=sys.stderr)
            raise SystemExit()
            
        if threep_dist + start == sequence_lengths[scaf]:
            threep_info = 'scaf_end'
        elif gene_info[scaf][start + threep_dist] > 0:
            gene_id = gene_names[gene_info[scaf][start + threep_dist]]
            threep_info = 'downstream_watson|' + gene_id
        elif gene_info[scaf][start + threep_dist] < 0:
            gene_id = gene_names[abs(gene_info[scaf][start + threep_dist])]
            threep_info = 'downstream_crick|' + gene_id
        else:
            # if this is triggered, serious debugging is needed
            print ("[{}] invalid threep_info value encountered!".format(
                   time.asctime()), file=sys.stderr)
            raise SystemExit()
        
        extra_info += [fivep_info, threep_info]
    else:
        gene_id = gene_names[abs(gene_info[scaf][start])]
        
        extra_info.append(gene_id)
        if gene_info[scaf][start] > 0:
            extra_info += calc_displacement(gene_info[scaf], [start, end])
        elif gene_info[scaf][start] < 0:
            extra_info += calc_displacement(gene_info[scaf], [end, start])
    
        # add exon info. ignore intergenic stuff - no exon/intron there.
        extra_info.append(translate_exon_info(exon_info[scaf][start], gene_id))
        if exon_info[scaf][start] > 0:
            extra_info += calc_displacement(exon_info[scaf], [start, end])
        elif exon_info[scaf][start] < 0:
            extra_info += calc_displacement(exon_info[scaf], [end, start])
        
    # fix unneeded precision in meth %.
    output = row[:3] + [str(round(meth/(meth+unmeth) * 100, 4))] + row[4:6]
    print ('\t'.join(output + extra_info))
    
    if args.verbose:
        counter += 1
        if counter % 10000 == 0:
            print ('[{}] {} methylated positions processed...'.format(
                   time.asctime(), counter), file=sys.stderr)
