#!/usr/bin/python

'''
Run AUGUSTUS for gene prediction with ab initio model.

* Used Augustus parameters in FunGAP
augustus\
    --uniqueGeneId=true\
    --singlestrand=true\
    --gff3=on\
    --species=<SPECIES_ARG>\
    --stopCodonExcludedFromCDS=false\
    --softmasking=1\
    <FASTA_FILE>\
    > <OUTPUT_GFF3>

--singlestrand: Predict genes independently on each strand. This makes maximal
    prediction including slight overlap between two neighboring genes on
    opposite strand.

Input: masked assembly and species parameter for Augustus
Output: gene features in GFF3
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
augustus_bin = os.path.join(this_dir, 'external/augustus-3.2.1/bin/augustus')


def main(argv):
    optparse_usage = (
        'run_augustus.py -i <input_fasta> -o <output_dir> '
        '-s <species> -l <log_dir>'
    )
    parser = ArgumentParser(usage=optparse_usage)
    parser.add_argument(
        "-i", "--input_fasta", dest="input_fasta", nargs=1,
        help="Input fasta file"
    )
    parser.add_argument(
        "-o", "--output_dir", dest="output_dir", nargs=1,
        help=(
            'Resulting gff3 and intermediates files will be stored '
            'in this directory'
        )
    )
    parser.add_argument(
        "-s", "--species", dest="species", nargs=1,
        help="Augustus reference species"
    )
    parser.add_argument(
        "-l", "--log_dir", dest="log_dir", nargs=1,
        help='Log directory'
    )

    args = parser.parse_args()
    if args.input_fasta:
        input_fasta = os.path.abspath(args.input_fasta[0])
    else:
        print '[ERROR] Please provide INPUT FASTA'
        sys.exit(2)

    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir[0])
    else:
        print '[ERROR] Please provide OUTPUT DIRECTORY'
        sys.exit(2)

    if args.species:
        species = args.species[0]
    else:
        print '[ERROR] Please provide SPECIES'
        sys.exit(2)

    if args.log_dir:
        log_dir = os.path.abspath(args.log_dir[0])
    else:
        print '[ERROR] Please provide LOG DIRECTORY'
        sys.exit(2)

    # Create necessary dirs
    create_dir(output_dir, log_dir)

    # Set logging
    log_file = os.path.join(
        log_dir, 'pipeline', 'run_augustus.log'
    )
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is as good as Fast
    run_augustus(input_fasta, output_dir, species)
    parse_augustus(output_dir)


# Define functions
def import_file(input_file):
    with open(input_file) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt


def create_dir(output_dir, log_dir):
    if not glob(output_dir):
        os.mkdir(output_dir)

    if not glob(log_dir):
        os.mkdir(log_dir)

    log_pipeline_dir = os.path.join(log_dir, 'pipeline')
    if not glob(log_pipeline_dir):
        os.mkdir(log_pipeline_dir)


def run_augustus(input_fasta, output_dir, species):
    # augustus --uniqueGeneId=true --gff3=on Neucr2_AssemblyScaffolds.fasta
    # --species=fusarium_graminearum --stopCodonExcludedFromCDS=false
    # > Neucr2.gff3

    augustus_output = os.path.join(output_dir, 'augustus.gff3')

    # Run AUGUSTUS
    logger_time.debug('START: Augustus')
    if not glob(augustus_output):
        command = (
            '%s --uniqueGeneId=true --singlestrand=true --gff3=on %s '
            '--species=%s --stopCodonExcludedFromCDS=false --softmasking=1 '
            '> %s'
        ) % (
            augustus_bin, input_fasta, species, augustus_output
        )
        logger_txt.debug('[Run] %s' % (command))
        os.system(command)
    else:
        logger_txt.debug('Running Augustus has already been finished')
    logger_time.debug('DONE : Augustus')


def parse_augustus(output_dir):
    augustus_gff3_file = os.path.join(output_dir, 'augustus.gff3')
    augustus_gff3 = import_file(augustus_gff3_file)

    # Define regular expression
    reg_transcript = re.compile(r'\ttranscript\t.+ID=([^;]+)')
    reg_proSeq_start = re.compile(r'^# protein sequence = \[(\S+)\]*')
    reg_proSeq_end = re.compile(r'\]$')

    prot_tag = 0
    D_seq = defaultdict(str)
    for line in augustus_gff3:
        # Exclude comment lines of BRAKER1 output
        if re.search('# Evidence for and against this transcript:', line):
            continue
        elif re.search('# % of transcript supported by hints', line):
            continue
        elif re.search('# CDS exons', line):
            continue
        elif re.search('# CDS introns', line):
            continue
        elif re.search("# 5'UTR exons and introns:", line):
            continue
        elif re.search("# 3'UTR exons and introns:", line):
            continue
        elif re.search("# hint groups fully obeyed:", line):
            continue
        elif re.search("# incompatible hint groups:", line):
            continue
        elif re.search("#      E:", line):
            continue
        elif re.search("#     RM:", line):
            continue

        m_transcript = reg_transcript.search(line)
        m_proSeq_start = reg_proSeq_start.search(line)
        m_proSeq_end = reg_proSeq_end.search(line)

        if m_transcript:
            transcript_id = m_transcript.group(1)
        elif m_proSeq_start:
            prot_tag = 1

        if m_proSeq_end:
            prot_seq = line.replace('# protein sequence = [', '')
            prot_seq = prot_seq.replace('# ', '').replace(']', '')
            D_seq[transcript_id] += prot_seq
            prot_tag = 0

        if prot_tag == 1:
            prot_seq = (
                line.replace('# protein sequence = [', '')
                    .replace('# ', '')
            )
            D_seq[transcript_id] += prot_seq

    # Write to file
    outfile = os.path.join(output_dir, 'augustus.faa')
    outhandle = open(outfile, "w")
    D_seq_sorted = sorted(
        D_seq.items(),
        key=lambda x: int(re.search(r'g(\d+)\.t\d+$', x[0]).group(1))
    )
    for transcript_id, prot_seq in D_seq_sorted:
        header_txt = '>%s\n' % (transcript_id)
        outhandle.write(header_txt)
        i = 0
        while i < len(prot_seq):
            row_txt = '%s\n' % (prot_seq[i:i + 60])
            outhandle.write(row_txt)
            i += 60


if __name__ == "__main__":
    main(sys.argv[1:])
