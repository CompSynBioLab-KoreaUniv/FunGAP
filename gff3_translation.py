#!/usr/bin/env python3

'''
Given fasta and gff3, get protein sequence in fasta foramt

Last updated: May 18, 2021
'''

import os
import re
import warnings
from argparse import ArgumentParser
from collections import defaultdict

from Bio import BiopythonWarning, SeqIO
from Bio.Seq import Seq

warnings.simplefilter('ignore', BiopythonWarning)


def main():
    '''Main function'''
    argparse_usage = (
        'gff3_translation.py -a <asm_file> -g <gff3_file> -o <output_file>')
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-a', '--asm_file', nargs=1, required=True,
        help='Genome assembly file (FASTA)')
    parser.add_argument(
        '-g', '--gff3_file', nargs=1, required=True,
        help='GFF3 file')
    parser.add_argument(
        '-t', '--translation_table', nargs='?', default=1,
        help='Translation table')
    parser.add_argument(
        '-o', '--output_file', nargs=1, required=True,
        help='Output file')

    args = parser.parse_args()
    asm_file = os.path.abspath(args.asm_file[0])
    gff3_file = os.path.abspath(args.gff3_file[0])
    translation_table = args.translation_table
    output_file = os.path.abspath(args.output_file[0])

    # Run functions :) Slow is as good as Fast
    parse_gff3(asm_file, gff3_file, translation_table, output_file)


def import_file(input_file):
    '''Import file'''
    with open(input_file) as f_in:
        txt = list(line.rstrip() for line in f_in)
    return txt


def parse_gff3(asm_file, gff3_file, table, output_file):
    '''Parse GFF3 file'''
    # Read gff3
    gff3 = import_file(gff3_file)
    # Parse gff3 and store in dictionary
    d_gff3 = defaultdict(list)
    reg_parent = re.compile('Parent=([^;]+)')
    for line in gff3:
        if re.search('^#', line):  # Ignore comment
            continue
        line_split = line.split('\t')
        entry_type = line_split[2]
        if entry_type != 'CDS':  # Only consider 'CDS'
            continue
        scaffold = line_split[0]
        start = int(line_split[3])
        end = int(line_split[4])
        strand = line_split[6]
        phase = int(line_split[7])
        gene_id = line_split[8]
        gene_id = reg_parent.search(gene_id).group(1)
        d_gff3[gene_id].append((scaffold, start, end, strand, phase))
    # Parse fasta and store in dictionary
    d_fasta = SeqIO.to_dict(SeqIO.parse(asm_file, 'fasta'))
    # Extract sequence
    output = open(output_file, 'w')
    gene_ids = d_gff3.keys()
    for gene_id in gene_ids:
        feature = d_gff3[gene_id]
        feature_s = sorted(feature, key=lambda tup: tup[1])
        nuc_seq = ''
        pro_seq = ''
        for element in feature_s:  # Feature is a list of tuple
            scaffold = element[0]
            start = element[1]
            end = element[2]
            strand = element[3]
            nuc_seq += str(d_fasta[scaffold].seq)[start - 1:end]
        # If it is '-' strand, reverse the transcript
        if strand == '-':
            nuc_seq = get_reverse_complement(nuc_seq)
        # If phase is not 0, trim first few bases according to phase
        if strand == '+' and feature_s[0][4] != 0:
            codon_start = feature_s[0][4]
            nuc_seq = nuc_seq[codon_start:]  # Trimming
        elif strand == '-' and feature_s[-1][4] != 0:
            codon_start = feature_s[-1][4]
            nuc_seq = nuc_seq[codon_start:]
        # Translation
        pro_seq = translation(nuc_seq, table)
        # Write to file
        output.write('>{}\n'.format(gene_id))
        i = 0
        while i < len(pro_seq):
            output.write('{}\n'.format(pro_seq[i:i + 60]))
            i += 60
    output.close()


def translation(seq, table):
    '''Translation'''
    seq_obj = Seq(seq)
    pro_seq = str(seq_obj.translate(table=table))
    pro_seq2 = re.sub(r'\*$', '', pro_seq)
    return pro_seq2


def get_reverse_complement(nuc_seq):
    '''Get reverse complement'''
    my_dna = Seq(nuc_seq)
    rev_comp_dna = str(my_dna.reverse_complement())
    return rev_comp_dna


if __name__ == "__main__":
    main()
