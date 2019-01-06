#!/usr/bin/env python2

'''
Run Blastp to databases
'''

# Import modules
import sys
import re
import os
from glob import glob
from argparse import ArgumentParser
from collections import defaultdict

# Get Logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)
sys.path.append(this_dir)
from set_logging import set_logging
from import_config import import_config

# Parameters
D_conf = import_config(this_dir)
evalue_cut = 0.00001
program_name = 'blastp'


# Main function
def main(argv):
    argparse_usage = (
        'run_blast_reduce.py -q <query_fasta> -d <db_fasta> '
        '-l <log_dir> -c <num_cores>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-q', '--query_fasta', nargs=1, required=True,
        help='Query FASTA file'
    )
    parser.add_argument(
        '-d', '--db_fasta', nargs=1, required=True,
        help='Database FASTA files'
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
    log_dir = args.log_dir
    num_cores = args.num_cores

    # Check input FASTA is valid
    if not glob(query_fasta):
        print '[ERROR] No such file: {}'.format(query_fasta)
        sys.exit(2)

    # Create necessary dirs
    create_dir(log_dir)

    # Set logging
    log_file = os.path.join(log_dir, 'run_blastp_reduce.log')
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is as good as Fast
    logger_time.debug('START: BLASTp')
    run_blastp(query_fasta, db_fasta, log_dir, num_cores)
    logger_time.debug('DONE : BLASTp')


def create_dir(log_dir):
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    log_program_dir = os.path.join(log_dir, program_name)
    if not os.path.exists(log_program_dir):
        os.mkdir(log_program_dir)


def run_blastp(query_fasta, db_fasta, log_dir, num_cores):
    # Run makeblastdb. Usage: makeblastdb -in <db_fasta> -dbtype prot
    makeblastdb_bin = D_conf['MAKEBLASTDB_PATH']
    blast_index_file = '{}.*phr'.format(db_fasta)
    log_file1 = os.path.join(log_dir, program_name, 'makeblastdb.log')
    if not glob(blast_index_file):
        command = (
            '{} -in {} -dbtype prot > {} 2>&1'.format(
                makeblastdb_bin, db_fasta, log_file1
            )
        )
        logger_txt.debug('[Run] {}'.format(command))
        os.system(command)
    else:
        logger_txt.debug('Running makeblastdb has been already finished')

    # Run BLASTp
    # blastp -outfmt "7 qseqid sseqid length qlen slen bitscore"
    # -query <query_fasta> -db <db_prefix> -out <out_file>
    # -num_threads <num_cores>
    blastp_bin = D_conf['BLASTP_PATH']
    input_base = os.path.splitext(query_fasta)[0]
    blastp_output = '{}.blastp'.format(input_base)
    log_file2 = os.path.join(log_dir, program_name, 'blastp.log')
    if not glob(blastp_output) or os.stat(blastp_output)[6] == 0:
        command = (
            '{} -outfmt "6 qseqid sseqid length qlen slen bitscore" -query '
            '{} -db {} -out {} -num_threads {} > {} 2>&1'.format(
                blastp_bin, query_fasta, db_fasta, blastp_output, num_cores,
                log_file2
            )
        )
        logger_txt.debug('[Run] {}'.format(command))
        os.system(command)


if __name__ == '__main__':
    main(sys.argv[1:])
