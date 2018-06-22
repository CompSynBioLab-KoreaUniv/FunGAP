#!/usr/bin/python

'''
Run BLASTn for given two FASTA files
Author Byoungnam Min on Aug 1, 2017
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


# Define main function
def main(argv):
    argparse_usage = (
        'run_blastn.py -q <query_fasta> -d <db_fasta> -o <output_prefix> '
        '-l <log_dir>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-q", "--query_fasta", dest="query_fasta", nargs=1,
        help="input fasta file"
    )
    parser.add_argument(
        "-d", "--db_fasta", dest="db_fasta", nargs=1,
        help="input fasta file"
    )
    parser.add_argument(
        "-o", "--output_prefix", dest="output_prefix", nargs=1,
        help="Output prefix"
    )
    parser.add_argument(
        "-l", "--log_dir", dest="log_dir", nargs=1,
        help="Log directory"
    )

    args = parser.parse_args()
    if args.query_fasta:
        query_fasta = os.path.abspath(args.query_fasta[0])
    else:
        print '[ERROR] Please provide QUERY FASTA'
        parser.print_help()
        sys.exit(2)

    if args.db_fasta:
        db_fasta = os.path.abspath(args.db_fasta[0])
    else:
        print '[ERROR] Please provide DB FASTA'
        parser.print_help()
        sys.exit(2)

    if args.output_prefix:
        output_prefix = os.path.abspath(args.output_prefix[0])
    else:
        print '[ERROR] Please provide OUTPUT PREFIX'
        parser.print_help()
        sys.exit(2)

    if args.log_dir:
        log_dir = os.path.abspath(args.log_dir[0])
    else:
        print '[ERROR] Please provide LOG DIRECTORY'
        parser.print_help()
        sys.exit(2)

    # Set logging
    log_file = os.path.join(
        log_dir, 'pipeline', 'run_blastn.log'
    )
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is as good as Fast
    logger_time.debug('START running BLASTn for %s' % (
        os.path.basename(query_fasta)
    ))

    # Run functions :) Slow is as good as Fast
    run_blastn(query_fasta, db_fasta, output_prefix)


def run_blastn(query_fasta, db_fasta, output_prefix):
    command1 = 'makeblastdb -in %s -dbtype nucl' % (db_fasta)
    logger_txt.debug('[Run] %s' % (command1))
    os.system(command1)

    out_blast = '%s.blast' % (output_prefix)
    command2 = (
        "blastn -query %s -db %s -out %s -outfmt '6 qseqid sseqid length "
        "qlen slen bitscore' -evalue 1e-5 -max_target_seqs 1" % (
            query_fasta, db_fasta, out_blast
        )
    )
    logger_txt.debug('[Run] %s' % (command2))
    if not os.path.exists(out_blast):
        os.system(command2)
    else:
        logger_txt.debug('Running BLASTn for %s has already been finished' % (
            os.path.basename(query_fasta)
        ))


if __name__ == "__main__":
    main(sys.argv[1:])
