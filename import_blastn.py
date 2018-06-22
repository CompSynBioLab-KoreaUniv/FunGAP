#!/usr/bin/python

'''
Import BLASTn result
 - Input: blastn output directory
 - Output: dictionary
Author Byoungnam Min on Aug 1, 2017
'''

# Import modules
from __future__ import division
import os
import sys
import cPickle
from glob import glob
from argparse import ArgumentParser
from collections import defaultdict


# Define main function
def main(argv):
    optparse_usage = 'import_blastn.py -b <blastn_dir>'
    parser = ArgumentParser(usage=optparse_usage)
    parser.add_argument(
        "-b", "--blastn_dir", dest="blastn_dir", nargs=1,
        help="BLASTn output directory (gpre_filtered/transcript)"
    )

    args = parser.parse_args()
    if args.blastn_dir:
        blastn_dir = os.path.abspath(args.blastn_dir[0])
    else:
        print '[ERROR] Please provide BLASTn DIRECTORY'
        sys.exit(2)

    # Run fuctions :) Slow is as good as Fast
    import_blastn(blastn_dir)


def import_file(input_file):
    with open(input_file) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt


def import_blastn(blastn_dir):
    blast_files = glob(os.path.join(blastn_dir, '*.blast'))
    D_blastn_score = defaultdict(float)
    for blast_file in blast_files:
        prefix = os.path.basename(blast_file).replace('.blast', '')
        blast_txt = import_file(blast_file)
        for line in blast_txt:
            line_split = line.split('\t')
            gene_id = line_split[0]
            alignment_length = int(line_split[2])
            qlen = int(line_split[3])
            slen = int(line_split[4])
            bit_score = float(line_split[5])

            q_cov = min(1, alignment_length / qlen)
            s_cov = min(1, alignment_length / slen)
            score = bit_score * q_cov * s_cov
            D_blastn_score[(prefix, gene_id)] += round(score, 1)

     # Write cPickle
    output_pickle = os.path.join(blastn_dir, 'blastn_score.p')
    cPickle.dump(D_blastn_score, open(output_pickle, 'wb'))


if __name__ == "__main__":
    main(sys.argv[1:])