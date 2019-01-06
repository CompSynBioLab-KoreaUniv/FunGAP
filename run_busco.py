#!/usr/bin/env python2

'''
Run BUSCO

* Please download fungi_buscos at path_to_fungap/data by typing
wget http://busco.ezlab.org/v1/files/fungi_buscos.tar.gz
tar -zxvf fungi_buscos.tar.gz
Currently version1 is supported.

Input: protein FASTA file
Output: BUSCO output in text
'''

# Import modules
import sys
import os
import re
from glob import glob
from argparse import ArgumentParser

# Get Logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)
sys.path.append(this_dir)
from set_logging import set_logging
from import_config import import_config

# Parameters
D_conf = import_config(this_dir)
lineage_path = D_conf['BUSCO_DB_PATH']
program_name = 'busco'


# Main function
def main(argv):
    argparse_usage = (
        'run_busco.py -i <input_fasta> -o <output_dir> -l <log_dir> '
        '-c <num_cores>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-i', '--input_fasta', nargs=1, required=True,
        help='Input protein FASTA file'
    )
    parser.add_argument(
        '-o', '--output_dir', nargs='?', default='busco_out',
        help='Output directory (default: busco_out)'
    )
    parser.add_argument(
        '-l', '--log_dir', nargs='?', default='logs',
        help='Log directory (default: logs)'
    )
    parser.add_argument(
        '-c', '--num_cores', nargs='?', default=1, type=int,
        help='Number of cores to be used (default: 1)'
    )

    args = parser.parse_args()
    input_fasta = os.path.abspath(args.input_fasta[0])
    output_dir = os.path.abspath(args.output_dir)
    log_dir = os.path.abspath(args.log_dir)
    num_cores = args.num_cores

    # Create necessary dir
    create_dir(output_dir, log_dir)

    # Set logging
    log_file = os.path.join(log_dir, 'run_busco.log')
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is always better than Fast
    run_busco(input_fasta, output_dir, log_dir, num_cores)


def import_file(input_file):
    with open(input_file) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt


def create_dir(output_dir, log_dir):
    # Output directory
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Log directory
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    # Log program directory
    log_output_dir = os.path.join(log_dir, program_name)
    if not os.path.exists(log_output_dir):
        os.mkdir(log_output_dir)


def run_busco(input_fasta, output_dir, log_dir, num_cores):
    busco_bin = D_conf['BUSCO_PATH']
    input_base = os.path.splitext(os.path.basename(input_fasta))[0]
    os.chdir(output_dir)
    logger_time.debug('START: BUSCO')

    busco_full_table = os.path.join(
        output_dir, 'run_{}'.format(input_base),
        'full_table_{}'.format(input_base)
    )

    if not os.path.exists(busco_full_table):
        log_file = os.path.join(
            log_dir, program_name, 'busco_{}.log'.format(input_base)
        )
        # BUSCO_v1.1b1.py -o NAME -in GENE_SET -l LINEAGE -m OGS
        command = (
            '{} --mode proteins --out {} --in {} --lineage_path {} --cpu {} > '
            '{} 2>&1'.format(
                busco_bin, input_base, input_fasta, lineage_path,
                num_cores, log_file
            )
        )
        logger_txt.debug('[Run] {}'.format(command))
        os.system(command)
    else:
        logger_txt.debug('Running BUSCO has already been finished')

    logger_time.debug('DONE : BUSCO')


if __name__ == '__main__':
    main(sys.argv[1:])
