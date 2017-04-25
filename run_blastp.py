#!/usr/bin/python

'''
Run Blastp to databases
    1) Run blastp to first dataset
    2) Filter >0.00001 hit
    3) Run blastp to second dataset
    4) Filter >0.00001 hit
    5) And so on..
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

# Parameters
evalue_cut = 0.00001


# Main function
def main(argv):
    optparse_usage = (
        'run_blast_reduce.py -i <input_fasta> -f <ref_fasta> '
        '-o <output> -c <num_cores> --nr'
    )
    parser = ArgumentParser(usage=optparse_usage)
    parser.add_argument(
        "-i", "--input_fasta", dest="input_fasta", nargs=1,
        help='input fasta file'
    )
    parser.add_argument(
        "-f", "--ref_fasta", dest="ref_fasta", nargs='*',
        help=(
            'Multiple reference FASTA files (order dependent, '
            'smallest dataset should be posed at first)'
        )
    )
    parser.add_argument(
        "-o", "--output_prefix", dest="output_prefix", nargs=1,
        help="output prefix"
    )
    parser.add_argument(
        "-r", "--root_dir", dest="root_dir", nargs=1,
        help=(
            'Root directory where log directory will be '
            'generated (default: ".")'
        ), default=[os.getcwd()]
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
        sys.exit(2)

    if args.ref_fasta:
        references = [os.path.abspath(x) for x in args.ref_fasta]
    else:
        references = []

    if args.output_prefix:
        output_prefix = args.output_prefix[0]
    else:
        print '[ERROR] Please provide OUTPUT_PREFIX'
        sys.exit(2)

    if args.num_cores:
        num_cores = args.num_cores[0]
    else:
        print '[ERROR] Please provide NUMBER OF CORES'
        sys.exit(2)

    root_dir = os.path.abspath(args.root_dir[0])

    # Check input fasta is valid
    if not glob(input_fasta):
        print '[ERROR] No such file: %s' % (input_fasta)
        sys.exit(2)

    # Create necessary dirs
    create_dir(root_dir)

    # Set logging
    log_file = os.path.join(
        root_dir, 'logs', 'pipeline', 'run_blastp_reduce.log'
    )
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is as good as Fast
    logger_time.debug('START running BLASTp-reduce for %s' % (
        os.path.basename(input_fasta)
    ))

    if references:
        filtered_fasta = input_fasta
        tmp_num = 1
        for ref in references:
            tmp_output_blast = run_blastp_ref(
                filtered_fasta, ref, output_prefix, tmp_num, num_cores
            )
            filtered_fasta, tmp_num = filtering(
                filtered_fasta, output_prefix, tmp_num, tmp_output_blast
            )
    else:
        filtered_fasta = input_fasta

    integrate(output_prefix, tmp_num)
    logger_time.debug('DONE  running BLASTp-reduce for %s' % (
        os.path.basename(input_fasta)
    ))


def import_file(input_file):
    with open(input_file) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt


def create_dir(root_dir):
    log_dir = os.path.join(root_dir, 'logs')
    if not glob(log_dir):
        os.mkdir(log_dir)

    log_pipeline_dir = os.path.join(root_dir, 'logs', 'pipeline')
    if not glob(log_pipeline_dir):
        os.mkdir(log_pipeline_dir)


def run_blastp_ref(filtered_fasta, ref, output_prefix, tmp_num, num_cores):
    # makeblastdb -in no_hit_KUC.fasta -dbtype nucl
    # blastn -query /csbl/genome/reference/ecoli_16s/ecoli_16s.fa -db
    # Contigs.1.fa -out Contigs.1.blast -outfmt 7 -evalue 1e-5 -num_threads 6
    # -dbsize 10000000

    # If blast-index was not generated, make one
    blast_index_file = '%s.*phr' % (ref)
    if not glob(blast_index_file):
        command = 'makeblastdb -in %s -dbtype prot' % (ref)
        logger_txt.debug('[Run] %s' % (command))
        os.system(command)

    # Run BLASTp
    tmp_output_blast = '%s.blast.%d' % (output_prefix, tmp_num)
    if not glob(tmp_output_blast) or os.stat(tmp_output_blast)[6] == 0:
        command = 'blastp -query %s -db %s -out %s -num_threads %s' % (
            filtered_fasta, ref, tmp_output_blast, num_cores
        )
        logger_txt.debug('[Run] %s' % (command))
        os.system(command)

    return tmp_output_blast


def filtering(filtered_fasta, output_prefix, tmp_num, tmp_output_blast):
    # Read & parse tmp_output_blast
    # Regular expressions
    reg_query = re.compile('Query= (\S+)')
    reg_evalue = re.compile(r'Expect = (\S+),')

    # Initialization
    start_flag = 0
    D_evalue = defaultdict(lambda: 100)

    blast = import_file(tmp_output_blast)
    for line in blast:
        m_query = reg_query.search(line)
        if m_query:
            query = m_query.group(1)
            start_flag = 1

        # To consider only best hit entry
        if start_flag == 0:
            continue

        m_evalue = reg_evalue.search(line)
        if m_evalue:
            evalue = float(m_evalue.group(1))
            D_evalue[query] = evalue
            start_flag = 0

    # Read and write faa
    tmp_output_faa = '%s.faa.%d' % (output_prefix, tmp_num)
    output = open(tmp_output_faa, 'w')
    fasta_txt = import_file(filtered_fasta)
    write_flag = 0
    for line in fasta_txt:
        if re.search('^>', line):
            name = line.split(' ')[0].replace('>', '')
            if D_evalue[name] > evalue_cut:
                write_flag = 1
                output.write('%s\n' % line)
            else:
                write_flag = 0

        elif write_flag == 1:
            output.write('%s\n' % line)

    output.close()

    tmp_num += 1
    return tmp_output_faa, tmp_num


def integrate(output_prefix, tmp_num):
    final_output_file = '%s.blast' % (output_prefix)
    command = 'cat %s.blast.[0-9] > %s' % (output_prefix, final_output_file)
    logger_txt.debug('[Run] %s' % (command))
    os.system(command)


if __name__ == "__main__":
    main(sys.argv[1:])
