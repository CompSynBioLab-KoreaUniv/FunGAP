#!/usr/bin/env python3

'''
Catch bad genes for given gff3 files
  1) Stop codon in the middle of proteins
  2) Check if translation consists of more than 50% X residues
  3) Check if feature begins or ends in gap

Input: multiple gff3s
Output: pickle for filter_gff3s_ver3.py
Last updated: Jul 13, 2020
'''

import operator
import os
import pickle
import re
import warnings
from argparse import ArgumentParser
from collections import defaultdict
from functools import reduce

from BCBio import GFF
from Bio import BiopythonWarning, SeqIO

warnings.simplefilter('ignore', BiopythonWarning)

SHORT_LEN = 10


def main():
    '''Main function'''
    argparse_usage = (
        'catch_bad_genes.py -g <gff3_files> -a <genome_assembly> '
        '-o <output_dir>')
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-g', '--gff3_files', nargs='+', required=True,
        help='Input GFF3 files')
    parser.add_argument(
        '-a', '--genome_assembly', nargs=1, required=True,
        help='Non-masked genome sequence file in FASTA')
    parser.add_argument(
        '-o', '--output_dir', nargs='?', default='gene_filtering',
        help='Output directory')

    args = parser.parse_args()
    gff3_files = [os.path.abspath(x) for x in args.gff3_files]
    genome_assembly_file = os.path.abspath(args.genome_assembly[0])
    output_dir = os.path.abspath(args.output_dir)

    # Run functions :) Slow is as good as Fast
    create_dir(output_dir)
    catch_middle_stop(gff3_files, genome_assembly_file, output_dir)


def create_dir(output_dir):
    '''Create directory'''
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


def catch_middle_stop(gff3_files, genome_assembly_file, output_dir):
    '''Catch genes having stop codon in the middle'''
    d_bad = defaultdict(bool)
    d_stop = defaultdict(int)
    d_toomany_x = defaultdict(int)
    d_gap = defaultdict(int)
    d_intron = defaultdict(int)
    d_short = defaultdict(int)
    for gff3_file in gff3_files:
        prefix = os.path.basename(os.path.splitext(gff3_file)[0])

        # Import genome sequence
        in_seq_handle = open(genome_assembly_file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, 'fasta'))
        in_seq_handle.close()

        # Import GFF3
        in_handle = open(gff3_file)
        for rec in GFF.parse(in_handle, base_dict=seq_dict):
            gene_features = rec.features
            for gene_feature in gene_features:
                mrna_features = gene_feature.sub_features
                for mrna_feature in mrna_features:
                    mrna_sub_features = mrna_feature.sub_features
                    mrna_sub_features_s = sorted(
                        mrna_sub_features, key=lambda x: x.location.start)
                    seq_cds = []
                    coords = []
                    mrna_sub_features_s2 = []
                    for feature in mrna_sub_features_s:
                        if feature.type != 'CDS':
                            continue
                        mrna_sub_features_s2.append(feature)
                        seq_cds.append(rec.seq[
                            feature.location.start:
                            feature.location.end])
                        coords.append(
                            (feature.location.start, feature.location.end))
                    i = 1
                    while i < len(coords):
                        intron_start = coords[i - 1][1]
                        intron_end = coords[i][0]
                        intron_len = intron_end - intron_start
                        if intron_len < 10:
                            d_bad[(prefix, mrna_feature.id)] = True
                            d_intron[prefix] += 1
                        i += 1

                    gene_seq = reduce(operator.add, seq_cds)
                    # If strand is -, get reverse complementary sequence
                    if mrna_feature.strand == -1:
                        gene_seq = gene_seq.reverse_complement()
                        phase = mrna_sub_features_s2[-1].qualifiers['phase']
                    else:
                        phase = mrna_sub_features_s2[0].qualifiers['phase']

                    phase = int(phase[0])
                    gene_seq = gene_seq[phase:]
                    protein_seq = str(gene_seq.translate())

                    # Check protein seq has stop codon in the middle
                    protein_seq2 = re.sub(r'\*$', '', protein_seq)
                    count_stop = protein_seq2.count('*')
                    if count_stop > 0:
                        d_bad[prefix, mrna_feature.id] = True
                        d_stop[prefix] += 1

                    # Check if translation consists of more than 50% X residues
                    len_prot = len(protein_seq2)
                    len_x = protein_seq2.count('X')
                    if len_x / len_prot > 0.5:
                        d_bad[prefix, mrna_feature.id] = True
                        d_toomany_x[prefix] += 1

                    # Check if feature begins or ends in gap
                    gene_seq2 = str(gene_seq).lower()
                    if gene_seq2.startswith('n') or gene_seq2.endswith('n'):
                        d_bad[prefix, mrna_feature.id] = True
                        d_gap[prefix] += 1

                    # Check if tranlated protein sequence is too short
                    if len(protein_seq) <= SHORT_LEN:
                        d_bad[prefix, mrna_feature.id] = True
                        d_short[prefix] += 1

    outfile_stats = os.path.join(output_dir, 'bad_genes_stats.txt')
    outhandle_stats = open(outfile_stats, 'w')
    run_names = d_stop.keys()
    header_txt = '{}\t{}\n'.format('type', '\t'.join(run_names))
    outhandle_stats.write(header_txt)

    stop_list = [str(d_stop[x]) for x in run_names]
    toomany_x_list = [str(d_toomany_x[x]) for x in run_names]
    gap_list = [str(d_gap[x]) for x in run_names]
    intron_list = [str(d_intron[x]) for x in run_names]
    short_list = [str(d_short[x]) for x in run_names]

    outhandle_stats.write('internal_stop\t{}\n'.format('\t'.join(stop_list)))
    outhandle_stats.write('start_with_gap\t{}\n'.format('\t'.join(gap_list)))
    outhandle_stats.write('toomanyX\t{}\n'.format('\t'.join(toomany_x_list)))
    outhandle_stats.write('short_intron\t{}\n'.format('\t'.join(intron_list)))
    outhandle_stats.write(
        'too_short_proteins\t{}\n'.format('\t'.join(short_list))
    )
    d_bad_pickle = os.path.join(output_dir, 'd_bad.p')
    pickle.dump(d_bad, open(d_bad_pickle, 'wb'))


if __name__ == '__main__':
    main()
