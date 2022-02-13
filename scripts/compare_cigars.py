#!/usr/bin/python

import sys

path_original_paf = sys.argv[1]
path_paffed_chain_paf = sys.argv[2]

with open(path_original_paf) as f1, open(path_paffed_chain_paf) as f2:
    for (line1,line2) in zip(f1, f2):
        # todo make it general (look for the position of cg: Z: and retrieve the cigar, or apply a regular expression)
        cigar1 = line1.strip().split('\t')[19].split('cg:Z:')[-1]
        cigar2 = line2.strip().split('\t')[12].split('cg:Z:')[-1]
        if cigar1 in cigar2 or cigar2 in cigar1:
            pass # Ok
        else:
            print('Problem')
            print('\t', line1.strip())
            print('\t', line2.strip())
