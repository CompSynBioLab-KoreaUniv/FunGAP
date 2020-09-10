#!/usr/bin/env python3

'''
Run RepeatModeler. The output of repeat models are passed into Maker

Input: genome assembly in FASTA
Output: Repeat model in FASTA (named consensi.fa.classified)
Last updated: Aug 12, 2020
'''

import os
import sys
from argparse import ArgumentParser
from glob import glob

from import_config import import_config
from set_logging import set_logging


def main():
    '''Main function'''
    argparse_usage = 'run_repeat_modeler.py -g <genome_assembly>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-g', '--genome_assembly', nargs=1, required=True,
        help='Genome assembly file in FASTA format'
    )
    parser.add_argument(
        '-o', '--output_dir', nargs='?', default='repeat_modeler_out',
        help='Output directory (default: repeat_modeler_out)'
    )
    parser.add_argument(
        '-l', '--log_dir', nargs='?', default='logs',
        help='Log directory (default: logs)'
    )
    parser.add_argument(
        '-c', '--num_cores', nargs='?', default=1, type=int,
        help='Number of cores to be used'
    )

    args = parser.parse_args()
    genome_assembly = os.path.abspath(args.genome_assembly[0])
    output_dir = os.path.abspath(args.output_dir)
    log_dir = os.path.abspath(args.log_dir)
    num_cores = args.num_cores

    # Create necessary dirs
    create_dir(output_dir, log_dir)

    # Set logging
    log_file = os.path.join(
        log_dir, 'run_repeat_modeler.log'
    )
    logger = set_logging(log_file)

    # Run functions :) Slow is as good as Fast
    run_repeat_modeler(genome_assembly, output_dir, log_dir, num_cores, logger)


def create_dir(output_dir, log_dir):
    '''Create directories'''
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)


def run_repeat_modeler(genome_assembly, output_dir, log_dir, num_cores, logger):
    '''Run RepeatModeler'''
    d_conf = import_config()
    builddatabase_bin = d_conf['BUILDDATABASE_PATH']
    repeatmodeler_bin = d_conf['REPEATMODELER_PATH']

    # BuildDatabase -name Choanephora_cucurbitarum
    # ../Choanephora_cucurbitarum_assembly.fna
    # RepeatModeler -database Choanephora_cucurbitarum -pa 25

    # Get repeat model
    repeat_lib = os.path.join(output_dir, '*', 'consensi.fa.classified')
    logger_time, logger_txt = logger
    if not glob(repeat_lib):
        os.chdir(os.path.join(output_dir))
        logger_time.debug('START running RepeatModeler')
        log_file1 = os.path.join(log_dir, 'build_database.log')
        command1 = '{0} -name {1} {1} > {2} 2>&1'.format(
            builddatabase_bin, genome_assembly, log_file1
        )
        logger_txt.debug('[Run] %s', command1)
        os.system(command1)

        log_file2 = os.path.join(log_dir, 'repeat_modeler.log')
        command2 = '{} -database {} -pa {} > {} 2>&1'.format(
            repeatmodeler_bin, genome_assembly, num_cores, log_file2
        )
        logger_txt.debug('[Run] %s', command2)
        os.system(command2)
        logger_time.debug('DONE  running RepeatModeler')
    else:
        logger_txt.debug(
            '[Note] Running RepeatModeler has already been finished'
        )

    # Check if RepeatModeler is properly finished
    if not glob(repeat_lib):
        logger_txt.debug(
            '[ERROR] RepeatModeler has finished abnormally. There is no '
            'consensi.fa.classified file.'
        )
        sys.exit(2)


if __name__ == '__main__':
    main()
