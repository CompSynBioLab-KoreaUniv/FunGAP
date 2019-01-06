#!/usr/bin/env python2

'''
Import BLASTn result
 - Input: blastn output directory
 - Output: dictionary
'''

# Import modules
from __future__ import division
import re
import os
import sys
import cPickle
from argparse import ArgumentParser
from collections import defaultdict


# Define main function
def main(argv):
    argparse_usage = 'import_blastn.py -b <blastn_out_files> -o <output_dir>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-b', '--blastn_out_files', nargs='+', required=True,
        help='BLASTn output files'
    )
    parser.add_argument(
        '-o', '--output_dir', nargs='?', default='gene_filtering',
        help='BLASTn output files'
    )

    args = parser.parse_args()
    blastn_out_files = [os.path.abspath(x) for x in args.blastn_out_files]
    output_dir = os.path.abspath(args.output_dir)

    # Run fuctions :) Slow is as good as Fast
    create_dir(output_dir)
    import_blastn(blastn_out_files, output_dir)


def import_file(input_file):
    with open(input_file) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt


def create_dir(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


def import_blastn(blastn_out_files, output_dir):
    D_blastn = defaultdict(float)
    for blast_file in blastn_out_files:
        prefix = re.sub('\.blastn$', '', os.path.basename(blast_file))
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
            D_blastn[(prefix, gene_id)] += round(score, 1)

     # Write cPickle
    output_pickle = os.path.join(output_dir, 'blastn_score.p')
    cPickle.dump(D_blastn, open(output_pickle, 'wb'))


if __name__ == '__main__':
    main(sys.argv[1:])