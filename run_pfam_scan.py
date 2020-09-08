#!/usr/bin/env python3

'''
Run Pfam_scan for Pfam domain identification on predicted genes

Input: protein FASTA file
Output: Identified Pfam domains .tsv format
Last updated: Aug 12, 2020
'''

import os
import re
import sys
from argparse import ArgumentParser
from collections import defaultdict

from import_config import import_config
from set_logging import set_logging

# Parameters
D_CONF = import_config()


# Main function
def main():
    '''Main function'''
    argparse_usage = 'run_pfam_scan.py -i <input_fasta> -l <log_dir>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-i', '--input_fasta', nargs=1, required=True,
        help='Input protein FASTA format'
    )
    parser.add_argument(
        '-l', '--log_dir', nargs='?', default='logs',
        help='Log directory'
    )
    parser.add_argument(
        '-c', '--num_cores', nargs='?', default=1, type=int,
        help='Number of cores to be used'
    )

    args = parser.parse_args()
    input_fasta = os.path.abspath(args.input_fasta[0])
    log_dir = os.path.abspath(args.log_dir)
    num_cores = args.num_cores

    # Create necessary dirs
    create_dir(log_dir)

    # Set logging
    log_file = os.path.join(log_dir, 'run_pfam_scan.log')
    logger = set_logging(log_file)

    # Run functions :) Slow is as good as fast
    new_input_fasta = check_sequence(input_fasta)
    run_pfam_scan(new_input_fasta, log_dir, num_cores, logger)


def create_dir(log_dir):
    '''Create directory'''
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)


def check_sequence(input_fasta):
    '''Check sequence'''
    with open(input_fasta) as f_in:
        fasta = list(line.rstrip() for line in f_in)

    d_seq = defaultdict(str)
    for line in fasta:
        if re.search('^>', line):
            gene_name = line.split('\t')[0].replace('>', '')
            continue
        d_seq[gene_name] += line

    new_input_fasta = '{}_nonX'.format(input_fasta)
    outhandle = open(new_input_fasta, 'w')
    for gene_name, seq in d_seq.items():
        if 'X' in seq or '*' in seq:
            continue
        i = 0
        outhandle.write('>{}\n'.format(gene_name))
        while i < len(seq):
            outhandle.write('{}\n'.format(seq[i: i + 60]))
            i += 60

    outhandle.close()
    return new_input_fasta


def run_pfam_scan(new_input_fasta, log_dir, num_cores, logger):
    '''Run pfam_scan'''
    # pfam_scan.pl -fasta Lenafn_TMI1502_1.faa -dir <pfam_db_dir> -cpu 10
    # -outfile pfam_scan.out
    pfam_scan_bin = D_CONF['PFAM_SCAN_PATH']
    pfam_db_dir = D_CONF['PFAM_DB_PATH']

    pfam_scan_out = '{}.pfam_scan'.format(os.path.splitext(new_input_fasta)[0])

    log_file = os.path.join(log_dir, 'pfam_scan.log')
    logger_time, logger_txt = logger
    logger_time.debug('START: Pfam Scan')
    if not os.path.exists(pfam_scan_out):
        command = (
            '{} -fasta {} -dir {} -cpu {} -outfile {} > {} 2>&1'.format(
                pfam_scan_bin, new_input_fasta, pfam_db_dir, num_cores,
                pfam_scan_out, log_file
            )
        )
        logger_txt.debug('[Run] %s', command)
        os.system(command)
    else:
        logger_txt.debug('[Note] Running Pfam Scan has already been finished')
    logger_time.debug('DONE : Pfam Scan')

    if not os.path.exists(pfam_scan_out):
        logger_txt.debug(
            '[ERROR] Pfam_scan failed to run. Please check the installation'
        )
        sys.exit(2)


if __name__ == '__main__':
    main()
