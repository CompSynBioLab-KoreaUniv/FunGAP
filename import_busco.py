#!/usr/bin/env python3

'''
Import BUSCO output and store in a dictionary
BUSCO evidence score is HMM alignment bit score

Input: BUSCO output
Output: pickle file containing dict object
Last updated: Aug 12, 2020
'''

import os
import pickle
from argparse import ArgumentParser
from collections import defaultdict
from glob import glob


def main():
    '''Main function'''
    argparse_usage = 'import_busco.py -b <busco_dir> -o <output_dir>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-b', '--busco_dir', nargs=1, required=True,
        help='BUSCO output directory (busco_out)'
    )
    parser.add_argument(
        '-o', '--output_dir', nargs='?', default='gene_filtering',
        help='Output directory (default: current working directory)'
    )

    args = parser.parse_args()
    busco_dir = os.path.abspath(args.busco_dir[0])
    output_dir = os.path.abspath(args.output_dir)

    # Run fuctions :) Slow is as good as Fast
    create_dir(output_dir)
    import_busco(busco_dir, output_dir)


def import_file(input_file):
    '''Import file'''
    with open(input_file) as f_in:
        txt = list(line.rstrip() for line in f_in)
    return txt


def create_dir(output_dir):
    '''Create directory'''
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


def import_busco(busco_dir, output_dir):
    '''Import BUSCO outputs'''
    # Because BUSCO output (full_table) doesn't have E-value
    # And raw HMM output has this, this script directly parse them
    busco_outdirs = glob(os.path.join(busco_dir, '*/run_*'))

    d_busco = defaultdict(float)
    d_score_element = {}
    for busco_outdir in busco_outdirs:
        prefix = busco_outdir.split('/')[-2]
        busco_hmmer_path = os.path.join(busco_outdir, 'hmmer_output/*/*out')
        busco_out_files = glob(busco_hmmer_path)

        # Parse E-value
        for busco_out_file in busco_out_files:
            busco_txt = import_file(busco_out_file)
            for line in busco_txt:
                # Ignore comment
                if line.startswith('#'):
                    continue

                line_split = line.split()
                gene_id = line_split[0]
                tlen = float(line_split[2])
                qlen = float(line_split[5])
                len_tup = (tlen, qlen)
                len_ratio = min(len_tup) / max(len_tup)
                full_seq_score = float(line_split[7])
                score = full_seq_score * len_ratio

                if score > d_busco[(prefix, gene_id)]:
                    # Update if having greater score
                    d_busco[(prefix, gene_id)] = score
                    d_score_element[(prefix, gene_id)] = (
                        full_seq_score, round(len_ratio, 3), round(score, 1)
                    )

    # Write to file
    outfile = os.path.join(output_dir, 'busco_score.txt')
    outhandle = open(outfile, 'w')
    header_txt = '{}\t{}\t{}\t{}\t{}\n'.format(
        'prefix', 'gene_id', 'full_seq_score', 'len_ratio', 'busco_score'
    )
    outhandle.write(header_txt)
    for tup1, tup2 in d_score_element.items():
        prefix, gene_id = tup1
        full_seq_score, len_ratio, busco_score = tup2
        row_txt = '{}\t{}\t{}\t{}\t{}\n'.format(
            prefix, gene_id, full_seq_score, len_ratio, busco_score
        )
        outhandle.write(row_txt)
    outhandle.close()

    # Write pickle
    output_pickle1 = os.path.join(output_dir, 'busco_score.p')
    pickle.dump(d_busco, open(output_pickle1, 'wb'))


if __name__ == '__main__':
    main()
