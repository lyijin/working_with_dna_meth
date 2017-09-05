#!/usr/bin/env python3

"""
> parse_fasta.py <

Python script intends to be a helper function that can be called by other
scripts to handle FASTA/FASTQ files.

Script requires a _file object_ (not filename) and filetype. The output is in
the form of dict[annotation] = sequence (by default, toggle using reverse_output
to produce dict[sequence] = annotation). Can probably be altered
in the future to output a list if annotations are not required.

Supports opening gzipped files by setting the gzip=True flag.
"""
from collections import OrderedDict
from itertools import tee
import gzip

def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)

def get_single_sequence(sequence_file, file_format, gzip_compressed=False):
    """
    A generator that can handle:
    - multi-line FASTA
    - single-line FASTA
    - standard FASTQ (4 lines per read)
    
    Outputs a single sequence with its annotation.
    
    Does not resolve duplications in annotations - avoid!
    """
    if file_format == 'fasta':
        annot_symbol = '>'
    elif file_format == 'fastq':
        annot_symbol = '@'
    
    if gzip_compressed:
        sequence_file = gzip.open(sequence_file.name, 'rt')
    
    annot = ''
    for line, next_line in pairwise(sequence_file):
        if line[:1] == annot_symbol:
            annot = line[1:].strip()
            seq = next_line.strip()
        else:
            # ignore comment lines before the first annot_symbol
            if annot:
                # produce output if next line = start of next sequence
                if next_line[:1] == annot_symbol:
                    yield annot, seq
                else:
                    if file_format == 'fasta':
                        # handles multi-line FASTA
                        seq += next_line.strip()

    yield annot, seq

    
def get_all_sequences(sequence_file, file_format, gzip_compressed=False,
                      sequences_only=False, lengths_only=False, 
                      include_annots=[], exclude_annots=[], 
                      contain_annots=[]):
    assert file_format in ['fasta', 'fastq'], "Unexpected file format!"
    
    def check_inclusion(annot):
        """checks whether annot should be outputted or not"""
        if not include_annots and not exclude_annots and not contain_annots:
            return True
        
        if include_annots:
            if annot in include_annots:
                return True
        
        if exclude_annots:
            if annot not in exclude_annots:
                return True
        
        if contain_annots:
            if any(c in annot for c in contain_annots):
                return True
        
        return False
    
    if sequences_only:
        # returns a list of sequences only, thus dodging problem with duplicate
        # annotations
        return [y for x, y in get_single_sequence(sequence_file, file_format, gzip_compressed)
                if check_inclusion(x)]
    
    if lengths_only:
        # instead of returning sequences, just returns the lengths of the
        # sequences
        all_lens = OrderedDict([(x, len(y)) for x, y in 
                                 get_single_sequence(sequence_file, file_format, gzip_compressed) 
                                 if check_inclusion(x)])
                            
        return all_lens
    
    # unfortunately this has problems when there are duplicate annots
    all_seqs = OrderedDict([x for x in \
                            get_single_sequence(sequence_file, file_format, gzip_compressed)
                            if check_inclusion(x)])
    
    return all_seqs