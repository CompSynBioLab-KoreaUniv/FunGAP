#!/usr/bin/python

'''
Run BUSCO
Author Byoungnam Min on Nov 17, 2015
'''

# Import modules
import sys
import os
import re
from argparse import ArgumentParser

# Get Logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)
sys.path.append(this_dir)
from set_logging import set_logging

# Parameters
lineage_path = os.path.join(this_dir, 'data/fungi')


def main(argv):
    optparse_usage = (
        'run_busco.py -i <input_fasta> -o <output_dir> -l <log_dir> '
        '-c <num_cores>'
    )
    parser = ArgumentParser(usage=optparse_usage)
    parser.add_argument(
        "-i", "--input_fasta", dest="input_fasta", nargs=1,
        help="Input protein FASTA file"
    )
    parser.add_argument(
        "-o", "--output_dir", dest="output_dir", nargs=1,
        help="Output directory"
    )
    parser.add_argument(
        "-l", "--log_dir", dest="log_dir", nargs=1,
        help='Log directory'
    )
    parser.add_argument(
        "-c", "--num_cores", dest="num_cores", nargs=1,
        help="Number of cores to be used"
    )

    args = parser.parse_args()
    if args.input_fasta:
        input_fasta = os.path.abspath(args.input_fasta[0])
    else:
        print '[ERROR] Please provide INPUT FASTA'
        parser.print_help()
        sys.exit(2)

    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir[0])
    else:
        print '[ERROR] Please provide OUTPUT DIRECTORY'
        parser.print_help()
        sys.exit(2)

    if args.log_dir:
        log_dir = os.path.abspath(args.log_dir[0])
    else:
        print '[ERROR] Please provide LOG DIRECTORY'
        parser.print_help()
        sys.exit(2)

    if args.num_cores:
        num_cores = args.num_cores[0]
    else:
        print '[ERROR] Please provide NUMBER OF CORES'
        parser.print_help()
        sys.exit(2)

    # Create necessary dirs
    create_dir(output_dir, log_dir)

    # Set logging
    log_file = os.path.join(
        log_dir, 'pipeline', 'run_busco.log'
    )
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is always better than Fast
    run_busco(input_fasta, output_dir, log_dir, num_cores)


def create_dir(output_dir, log_dir):
    output_dir = re.sub('/$', '', output_dir)  # Trim trailing "/"
    output_base = os.path.basename(output_dir)
    # Output dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Log directory
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    # Log output directory
    log_output_dir = os.path.join(log_dir, output_base)
    if not os.path.exists(log_output_dir):
        os.mkdir(log_output_dir)

    # Log pipeline direcotry
    log_pipeline_dir = os.path.join(log_dir, 'pipeline')
    if not os.path.exists(log_pipeline_dir):
        os.mkdir(log_pipeline_dir)


def run_busco(input_fasta, output_dir, log_dir, num_cores):
    input_base = os.path.splitext(os.path.basename(input_fasta))[0]
    output_dir = re.sub('/$', '', output_dir)  # Trim trailing "/"
    output_base = os.path.basename(output_dir)

    os.chdir(output_dir)
    logger_time.debug('START: BUSCO for %s' % (input_base))

    busco_full_table = os.path.join(
        output_dir, 'run_%s' % (input_base),
        'full_table_%s' % (input_base)
    )

    if not os.path.exists(busco_full_table):
        log_file = os.path.join(
            log_dir, output_base, 'busco_%s.log' % (input_base)
        )
        # BUSCO_v1.1b1.py -o NAME -in GENE_SET -l LINEAGE -m OGS
        command = (
            'BUSCO_v1.1b1.py -o %s -in %s -l %s -m %s -c %s > %s 2>&1'
        ) % (
            input_base, input_fasta, lineage_path, 'OGS', num_cores, log_file
        )
        logger_txt.debug('[Run] %s' % (command))
        os.system(command)
    else:
        logger_txt.debug('Running BUSCO has already been finished')

    logger_time.debug('DONE : BUSCO for %s' % (output_base))

if __name__ == "__main__":
    main(sys.argv[1:])
