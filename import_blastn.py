#!/usr/bin/env python3

'''
Import BLASTn result
 - Input: blastn output directory
 - Output: dictionary
Last updated: Aug 12, 2020
'''

import os
import pickle
import re
from argparse import ArgumentParser
from collections import defaultdict


def main():
    '''Main function'''
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
    '''Import file'''
    with open(input_file) as f_in:
        txt = list(line.rstrip() for line in f_in)
    return txt


def create_dir(output_dir):
    '''Create directory'''
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


def import_blastn(blastn_out_files, output_dir):
    '''Import BLASTn output'''
    d_blastn = defaultdict(float)
    for blast_file in blastn_out_files:
        prefix = re.sub(r'\.blastn$', '', os.path.basename(blast_file))
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
            d_blastn[(prefix, gene_id)] += round(score, 1)

     # Write pickle
    output_pickle = os.path.join(output_dir, 'blastn_score.p')
    pickle.dump(d_blastn, open(output_pickle, 'wb'))


if __name__ == '__main__':
    main()
