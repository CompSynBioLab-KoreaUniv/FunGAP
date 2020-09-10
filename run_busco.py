#!/usr/bin/env python3

'''
Run BUSCO (v4.1.2 tested)

- fungi_odb10
  - ascomycota_odb10
    - dothideomycetes_odb10
      - capnodiales_odb10
      - pleosporales_odb10
    - eurotiomycetes_odb10
      - chaetothyriales_odb10
      - eurotiales_odb10
      - onygenales_odb10
    - leotiomycetes_odb10
      - helotiales_odb10
    - saccharomycetes_odb10
    - sordariomycetes_odb10
      - glomerellales_odb10
      - hypocreales_odb10
  - basidiomycota_odb10
    - agaricomycetes_odb10
      - agaricales_odb10
      - boletales_odb10
      - polyporales_odb10
    - tremellomycetes_odb10
  - microsporidia_odb10
  - mucoromycota_odb10
    - mucorales_odb10

Input: protein FASTA file
Output: BUSCO output in text
Last updated: Aug 12, 2020
'''

import os
from argparse import ArgumentParser

from import_config import import_config
from set_logging import set_logging

# Parameters
D_CONF = import_config()


def main():
    '''Main function'''
    argparse_usage = ('run_busco.py -i <input_fasta> -d <lineage_dataset>')
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-i', '--input_fasta', nargs=1, required=True,
        help='Input protein FASTA file'
    )
    parser.add_argument(
        '-d', '--lineage_dataset', nargs=1, required=True,
        help='BUSCO lineage dataset (run "busco --list-datasets" for the list)'
    )
    parser.add_argument(
        '-o', '--output_dir', nargs='?', default='busco_out',
        help='Output directory (default: busco_out)'
    )
    parser.add_argument(
        '-l', '--log_dir', nargs='?', default='logs',
        help='Log directory (default: logs)'
    )

    args = parser.parse_args()
    input_fasta = os.path.abspath(args.input_fasta[0])
    lineage_dataset = args.lineage_dataset[0]
    output_dir = os.path.abspath(args.output_dir)
    log_dir = os.path.abspath(args.log_dir)

    # Create necessary dir
    create_dir(output_dir, log_dir)

    # Set logging
    log_file = os.path.join(log_dir, 'run_busco.log')
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is always better than Fast
    log_tup = (log_dir, logger_time, logger_txt)
    run_busco(input_fasta, lineage_dataset, output_dir, log_tup)


def import_file(input_file):
    '''Import file'''
    with open(input_file) as f_in:
        txt = list(line.rstrip() for line in f_in)
    return txt


def create_dir(output_dir, log_dir):
    '''Create directories'''
    # Output directory
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Log directory
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)


def run_busco(input_fasta, lineage_dataset, output_dir, log_tup):
    '''Run BUSCO'''
    log_dir, logger_time, logger_txt = log_tup
    busco_bin = D_CONF['BUSCO_PATH']
    input_base = os.path.splitext(os.path.basename(input_fasta))[0]
    logger_time.debug('START: BUSCO')

    busco_full_table = os.path.join(
        output_dir, input_base, 'run_{}'.format(lineage_dataset),
        'full_table.tsv'
    )

    if not os.path.exists(busco_full_table):
        log_file = os.path.join(
            log_dir, 'busco_{}.log'.format(input_base)
        )
        command = (
            '{} --mode proteins --out {} --in {} --out_path {} '
            '--lineage_dataset {} --force > {} 2>&1'.format(
                busco_bin, input_base, input_fasta, output_dir, lineage_dataset,
                log_file
            )
        )
        logger_txt.debug('[Run] %s', command)
        os.system(command)
    else:
        logger_txt.debug('[Note] Running BUSCO has already been finished')

    logger_time.debug('DONE : BUSCO')


if __name__ == '__main__':
    main()
