#!/usr/bin/env python2

'''
Run BLASTn for given two FASTA files
'''

# Import modules
import os
import sys
from argparse import ArgumentParser

# Get Logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)
sys.path.append(this_dir)
from set_logging import set_logging
from import_config import import_config

# Parameters
D_conf = import_config(this_dir)
program_name = 'blastn'


# Define main function
def main(argv):
    argparse_usage = (
        'run_blastn.py -q <query_fasta> -d <db_fasta> -o <output_prefix> '
        '-l <log_dir> -c <num_cores>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-q', '--query_fasta', nargs=1, required=True,
        help='Query FASTA file'
    )
    parser.add_argument(
        '-d', '--db_fasta', nargs=1, required=True,
        help='Database FASTA file'
    )
    parser.add_argument(
        '-o', '--output_prefix', nargs='?', default='out',
        help='Output prefix'
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
    query_fasta = os.path.abspath(args.query_fasta[0])
    db_fasta = os.path.abspath(args.db_fasta[0])
    output_prefix = os.path.abspath(args.output_prefix)
    log_dir = os.path.abspath(args.log_dir)
    num_cores = args.num_cores

    # Set logging
    create_dir(log_dir)

    log_file = os.path.join(log_dir, 'run_blastn.log')
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    logger_time.debug('START: BLASTn for {}'.format(
        os.path.basename(query_fasta))
    )
    # Run functions :) Slow is as good as Fast
    run_blastn(query_fasta, db_fasta, output_prefix, log_dir, num_cores)
    logger_time.debug('Done : BLASTn for {}'.format(
        os.path.basename(query_fasta))
    )


def create_dir(log_dir):
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    log_program_dir = os.path.join(log_dir, program_name)
    if not os.path.exists(log_program_dir):
        os.mkdir(log_program_dir)


def run_blastn(query_fasta, db_fasta, output_prefix, log_dir, num_cores):
    blastn_out = '{}.blastn'.format(output_prefix)
    if not os.path.exists(blastn_out):
        makeblastdb_bin = D_conf['MAKEBLASTDB_PATH']
        log_file1 = os.path.join(log_dir, program_name, 'makeblastdb.log')
        command1 = '{} -in {} -dbtype nucl > {} 2>&1'.format(
            makeblastdb_bin, db_fasta, log_file1
        )
        logger_txt.debug('[Run] {}'.format(command1))
        os.system(command1)

        blastn_bin = D_conf['BLASTN_PATH']
        log_file2 = os.path.join(log_dir, program_name, 'blastn.log')
        command2 = (
            '{} -query {} -db {} -out {} -outfmt "6 qseqid sseqid length '
            'qlen slen bitscore" -num_threads {} -evalue 1e-5 > {} 2>&1'.format(
                blastn_bin, query_fasta, db_fasta, blastn_out, num_cores,
                log_file2
            )
        )
        logger_txt.debug('[Run] {}'.format(command2))
        os.system(command2)
    else:
        logger_txt.debug('Running BLASTn has already been finished')


if __name__ == '__main__':
    main(sys.argv[1:])
