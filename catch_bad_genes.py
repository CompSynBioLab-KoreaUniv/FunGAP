#!/usr/bin/python

'''
Catch bad genes for given gff3 files
    1) Stop codon in the middle of proteins
    2) Check if translation consists of more than 50% X residues
    3) Check if feature begins or ends in gap

Input: multiple gff3s
Output: pickle for filter_gff3s_ver3.py
'''

# Import modeuls
from __future__ import division
import sys
import os
import re
import operator
import cPickle
from collections import defaultdict
from argparse import ArgumentParser
from BCBio import GFF
from Bio import SeqIO


# Main function
def main(argv):
    argparse_usage = (
        'catch_bad_genes.py -g <gff3_files> -f <fna_file> -o <output_dir>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-g", "--gff3_files", dest="gff3_files", nargs='+',
        help="Input GFF3 files"
    )
    parser.add_argument(
        "-f", "--fna_file", dest="fna_file", nargs=1,
        help="Genome sequence file in FNA"
    )
    parser.add_argument(
        "-o", "--output_dir", dest="output_dir", nargs=1,
        help="Output prefix without extension"
    )

    args = parser.parse_args()
    if args.gff3_files:
        gff3_files = [os.path.abspath(x) for x in args.gff3_files]
    else:
        print '[ERROR] Please provide GFF3 FILES'
        sys.exit(2)

    if args.fna_file:
        fna_file = os.path.abspath(args.fna_file[0])
    else:
        print '[ERROR] Please provide FNA FILE'
        sys.exit(2)

    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir[0])
    else:
        print '[ERROR] Please provide OUTPUT PREFIX'
        sys.exit(2)

    # Run functions :) Slow is as good as Fast
    catch_middle_stop(gff3_files, fna_file, output_dir)


def catch_middle_stop(gff3_files, fna_file, output_dir):
    D_bad = defaultdict(bool)
    D_stop = defaultdict(int)
    D_toomanyX = defaultdict(int)
    D_gap = defaultdict(int)
    D_intron = defaultdict(int)
    for gff3_file in gff3_files:
        prefix = os.path.basename(gff3_file).split('.')[0]

        # Import genome sequence
        in_seq_handle = open(fna_file)
        seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
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
                        mrna_sub_features, key=lambda x: x.location.start
                    )
                    seq_cds = []
                    coords = []
                    for feature in mrna_sub_features_s:
                        if feature.type != 'CDS':
                            continue
                        seq_cds.append(rec.seq[
                            feature.location.start:
                            feature.location.end])
                        coords.append(
                            (feature.location.start, feature.location.end)
                        )

                    i = 1
                    while i < len(coords):
                        intron_start = coords[i - 1][1]
                        intron_end = coords[i][0]
                        intron_len = intron_end - intron_start
                        if intron_len < 10:
                            D_bad[(prefix, mrna_feature.id)] = True
                            D_intron[prefix] += 1
                        i += 1

                    gene_seq = reduce(operator.add, seq_cds)
                    # If strand is -, get reverse complementary sequence
                    if mrna_feature.strand == -1:
                        gene_seq = gene_seq.reverse_complement()

                    protein_seq = str(gene_seq.translate())

                    # Check protein seq has stop codon in the middle
                    protein_seq2 = re.sub('\*$', '', protein_seq)
                    count_stop = protein_seq2.count('*')
                    if count_stop > 0:
                        D_bad[(prefix, mrna_feature.id)] = True
                        D_stop[prefix] += 1

                    # Check if translation consists of more than 50% X residues
                    len_prot = len(protein_seq2)
                    len_X = protein_seq2.count('X')
                    if len_X / len_prot > 0.5:
                        D_bad[(prefix, mrna_feature.id)] = True
                        D_toomanyX[prefix] += 1

                    # Check if feature begins or ends in gap
                    gene_seq2 = str(gene_seq).lower()
                    if gene_seq2.startswith('n') or gene_seq2.endswith('n'):
                        D_bad[(prefix, mrna_feature.id)] = True
                        D_gap[prefix] += 1

    in_handle.close()
    outfile_stats = os.path.join(output_dir, 'bad_genes_stats.txt')
    outhandle_stats = open(outfile_stats, 'w')
    run_names = D_stop.keys()
    header_txt = '%s\t%s\n' % ('type', '\t'.join(run_names))
    outhandle_stats.write(header_txt)

    stop_list = [str(D_stop[x]) for x in run_names]
    toomanyX_list = [str(D_toomanyX[x]) for x in run_names]
    gap_list = [str(D_gap[x]) for x in run_names]
    intron_list = [str(D_intron[x]) for x in run_names]

    outhandle_stats.write('internal_stop\t%s\n' % ('\t'.join(stop_list)))
    outhandle_stats.write('start_with_gap\t%s\n' % ('\t'.join(gap_list)))
    outhandle_stats.write('toomanyX\t%s\n' % ('\t'.join(toomanyX_list)))
    outhandle_stats.write('short_intron\t%s\n' % ('\t'.join(intron_list)))
    D_bad_pickle = os.path.join(output_dir, 'D_bad.p')
    cPickle.dump(D_bad, open(D_bad_pickle, 'wb'))


if __name__ == "__main__":
    main(sys.argv[1:])
