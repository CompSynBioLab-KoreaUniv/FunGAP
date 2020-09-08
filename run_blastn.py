#!/usr/bin/env python3

'''
Run BLASTn for given two FASTA files
Last updated: Aug 12, 2020
'''

import os
from argparse import ArgumentParser

from import_config import import_config
from set_logging import set_logging

# Parameters
D_CONF = import_config()


# Define main function
def main():
    '''Main function'''
    argparse_usage = (
        'run_blastn.py -q <query_fasta> -d <db_fasta> -o <output_prefix> '
        '-l <log_dir> -c <num_cores>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-q', '--query_fasta', nargs=1, required=True, help='Query FASTA file'
    )
    parser.add_argument(
        '-d', '--db_fasta', nargs=1, required=True, help='Database FASTA file'
    )
    parser.add_argument(
        '-o', '--output_prefix', nargs='?', default='out', help='Output prefix'
    )
    parser.add_argument(
        '-l', '--log_dir', nargs='?', default='logs', help='Log directory'
    )
    parser.add_argument(
        '-c', '--num_cores', nargs='?', default=1, type=int,
        help='Number of cores to be used'
    )

    args = parser.parse_args()
    query_fasta = os.path.abspath(args.query_fasta[0])
    db_fasta = os.path.abspath(args.db_fasta[0])
    output_prefix = os.path.abspath(args.output_prefix)
    log_dir = os.path.abspath(args.log_dir)
    num_cores = args.num_cores

    # Set logging
    create_dir(log_dir)

    log_file = os.path.join(log_dir, 'run_blastn.log')
    logger = set_logging(log_file)
    logger_time = logger[0]

    logger_time.debug('START: BLASTn for %s', os.path.basename(query_fasta))
    # Run functions :) Slow is as good as Fast
    run_blastn(query_fasta, db_fasta, output_prefix, log_dir, num_cores, logger)
    logger_time.debug('DONE : BLASTn for %s', os.path.basename(query_fasta))


def create_dir(log_dir):
    '''Create directory'''
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)


def run_blastn(
        query_fasta, db_fasta, output_prefix, log_dir, num_cores, logger):
    '''Run BLASTn'''
    blastn_out = '{}.blastn'.format(output_prefix)
    logger_txt = logger[1]
    if not os.path.exists(blastn_out):
        makeblastdb_bin = D_CONF['MAKEBLASTDB_PATH']
        log_file1 = os.path.join(log_dir, 'makeblastdb_blastn.log')
        command1 = '{} -in {} -dbtype nucl > {} 2>&1'.format(
            makeblastdb_bin, db_fasta, log_file1
        )
        logger_txt.debug('[Run] %s', command1)
        os.system(command1)

        blastn_bin = D_CONF['BLASTN_PATH']
        log_file2 = os.path.join(log_dir, 'blastn.log')
        command2 = (
            '{} -query {} -db {} -out {} -outfmt "6 qseqid sseqid length '
            'qlen slen bitscore" -num_threads {} -evalue 1e-5 > {} 2>&1'.format(
                blastn_bin, query_fasta, db_fasta, blastn_out, num_cores,
                log_file2
            )
        )
        logger_txt.debug('[Run] %s', command2)
        os.system(command2)
    else:
        logger_txt.debug('[Note] Running BLASTn has already been finished')


if __name__ == '__main__':
    main()
