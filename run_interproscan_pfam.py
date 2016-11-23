#!/usr/bin/python

'''
Run IPRscan
Author Byoungnam Min on Feb 15, 2015
 - UPDATE: ADD draw_barchart function - Jun 25, 2015
'''

# Import modules
import sys
import os
import re
from argparse import ArgumentParser
from collections import defaultdict

# Get Logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)
sys.path.append(this_dir)
from set_logging import set_logging


# Main fuction
def main(argv):
    optparse_usage = (
        'run_interproscan.py -i <input_fasta> -o <output_dir> -l <log_dir>'
    )
    parser = ArgumentParser(usage=optparse_usage)
    parser.add_argument(
        "-i", "--input_fasta", dest="input_fasta", nargs=1,
        help="Input protein FASTA format"
    )
    parser.add_argument(
        "-o", "--output_dir", dest="output_dir", nargs=1,
        help="Output directory"
    )
    parser.add_argument(
        "-l", "--log_dir", dest="log_dir", nargs=1,
        help="Log directory"
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

    # Create necessary dirs
    create_dir(output_dir, log_dir)

    # Set logging
    log_file = os.path.join(
        log_dir, 'pipeline', 'run_interproscan_pfam.log'
    )
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is as good as fast
    check_sequence(input_fasta)
    run_iprscan(input_fasta, output_dir, log_dir)


# Define functions
def create_dir(output_dir, log_dir):
    # Output directory
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # Temporary directory
    tmp_dir = os.path.join(output_dir, 'tmp')
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)

    output_base = os.path.basename(output_dir)
    # Log directory
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    # Log output directory
    log_output_dir = os.path.join(log_dir, output_base)
    if not os.path.exists(log_output_dir):
        os.mkdir(log_output_dir)

    # Log pipeline directory
    log_pipeline_dir = os.path.join(log_dir, 'pipeline')
    if not os.path.exists(log_pipeline_dir):
        os.mkdir(log_pipeline_dir)


def check_sequence(input_fasta):
    aa_char = 'ARNDCEQGHILKMFPSTWYV'
    aa_list = list(aa_char)
    with open(input_fasta) as f_in:
        fasta = (line.rstrip() for line in f_in)
        fasta = list(line for line in fasta if line)

    D = defaultdict(str)
    for line in fasta:
        if re.search('^>', line):
            gene_name = line.split('\t')[0].replace('>', '')
            continue
        D[gene_name] += line

    for gene_name, seq in D.items():
        non_aa = [x for x in seq if x not in aa_list]

    if non_aa:
        print 'You have wrong sequence'
        print non_aa
        sys.exit(2)


def run_iprscan(input_fasta, output_dir, log_dir):
    # interproscan.sh -i <protein.fasta> -f tsv --goterms --iprlookup
    # -b <base_name> --tempdir <TEMP-DIR>
    output_base = os.path.basename(output_dir)
    tmp_dir = os.path.join(output_dir, 'tmp')
    input_base = os.path.splitext(os.path.basename(input_fasta))[0]
    ipr_output = os.path.join(output_dir, input_base)
    ipr_tsv = os.path.join(output_dir, '%s.tsv' % (input_base))

    log_file = os.path.join(
        log_dir, output_base, '%s.log' % (output_base)
    )
    logger_time.debug('START: Interproscan for Pfam')
    if not os.path.exists(ipr_tsv):
        command = (
            'interproscan.sh -i %s --goterms -pa --iprlookup '
            '-f XML,tsv -appl PfamA --tempdir %s --output-file-base '
            '%s > %s' % (input_fasta, tmp_dir, ipr_output, log_file)
        )
        logger_txt.debug('[Run] %s' % (command))
        os.system(command)
    else:
        logger_txt.debug('Running Iprscan has already been finished')
    logger_time.debug('DONE : Interproscan for Pfam')


if __name__ == "__main__":
    main(sys.argv[1:])
