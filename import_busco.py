#!/usr/bin/python

'''
Import BUSCO output and store it in dictionary
BUSCO evidence score is HMM alignment bit score

Input: BUSCO output
Output: cPickle file containing dict object
'''

# Import modules
import os
import sys
import cPickle
from glob import glob
from argparse import ArgumentParser
from collections import defaultdict


def main(argv):
    optparse_usage = 'import_busco.py -b <busco_dir>'
    parser = ArgumentParser(usage=optparse_usage)
    parser.add_argument(
        "-b", "--busco_dir", dest="busco_dir", nargs=1,
        help="BUSCO output directory (gpre_busco)"
    )

    args = parser.parse_args()
    if args.busco_dir:
        busco_dir = os.path.abspath(args.busco_dir[0])
    else:
        print '[ERROR] Please provide BUSCO DIRECTORY'
        sys.exit(2)

    # Run fuctions :) Slow is as good as Fast
    import_busco(busco_dir)


def import_file(input_file):
    with open(input_file) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt


def import_busco(busco_dir):
    # Because BUSCO output (full_table) doesn't have E-value
    # And raw HMM output has this, this script directly parse them
    busco_outdirs = glob(os.path.join(busco_dir, 'run_*'))

    D_busco_score = defaultdict(float)
    D_score_element = {}
    D_busco_list = defaultdict(list)
    for busco_outdir in busco_outdirs:
        prefix = os.path.basename(busco_outdir).replace('run_', '')

        busco_hmmer_path = os.path.join(busco_outdir, 'hmmer_output/*out')
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

                if score > D_busco_score[(prefix, gene_id)]:
                    # Update if having greater score
                    D_busco_score[(prefix, gene_id)] = score
                    D_score_element[(prefix, gene_id)] = (
                        full_seq_score, round(len_ratio, 3), round(score, 1)
                    )

        # Get BUSCO complete and duplicated list
        full_table_file = os.path.join(
            busco_outdir, 'full_table_%s' % (prefix)
        )
        full_table_txt = import_file(full_table_file)
        for line in full_table_txt[1:]:
            line_split = line.split('\t')
            if line_split[1] == 'Missing':
                continue

            elif line_split[1] in ('Duplicated', 'Complete'):
                busco_id = line_split[0]
                gene_id = line_split[2]
                D_busco_list[(prefix, gene_id)].append(busco_id)

    # Write to file
    outfile = os.path.join(busco_dir, 'busco_parsed.txt')
    outhandle = open(outfile, 'w')
    header_txt = '%s\t%s\t%s\t%s\t%s\n' % (
        'prefix', 'gene_id', 'full_seq_score', 'len_ratio', 'busco_score'
    )
    outhandle.write(header_txt)
    for tup1, tup2 in D_score_element.items():
        prefix, gene_id = tup1
        full_seq_score, len_ratio, busco_score = tup2
        row_txt = '%s\t%s\t%s\t%s\t%s\n' % (
            prefix, gene_id, full_seq_score, len_ratio, busco_score
        )
        outhandle.write(row_txt)
    outhandle.close()

    # Write cPickle
    output_pickle1 = os.path.join(busco_dir, 'busco_score.p')
    cPickle.dump(D_busco_score, open(output_pickle1, 'wb'))

    output_pickle2 = os.path.join(busco_dir, 'busco_list.p')
    cPickle.dump(D_busco_list, open(output_pickle2, 'wb'))


if __name__ == "__main__":
    main(sys.argv[1:])
