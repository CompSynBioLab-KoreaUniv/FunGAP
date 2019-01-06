#!/usr/bin/env python2

'''
GFF3 postprocessing
    - Remove UTRs when two near genes are overlapped

Input: GFF3 file
Output: Postprocessed GFF3 file
'''

# Import modules
import os
import sys
from argparse import ArgumentParser

from BCBio import GFF
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import FeatureLocation


# Main function
def main(argv):
    argparser_usage = (
        'gff3_postprocess.py -g <genome_assembly> -i <input_gff3> -o '
        '<output_gff3>'
    )
    parser = ArgumentParser(usage=argparser_usage)
    parser.add_argument(
        "-g", "--genome_assembly", nargs=1, required=True,
        help='Genome assembly file in FASTA format'
    )
    parser.add_argument(
        "-i", "--input_gff3", nargs=1, required=True,
        help='Input GFF3 file'
    )
    parser.add_argument(
        "-o", "--output_gff3", nargs=1, required=True,
        help='Output GFF3 file name'
    )

    args = parser.parse_args()

    genome_assembly = os.path.abspath(args.genome_assembly[0])
    input_gff3 = os.path.abspath(args.input_gff3[0])
    output_gff3 = os.path.abspath(args.output_gff3[0])

    # Create necessary dirs
    gff3_postprocess(genome_assembly, input_gff3, output_gff3)


def import_file(input_file):
    with open(input_file) as f_in:
        txt = list(line.rstrip() for line in f_in)
    return txt


def gff3_postprocess(genome_assembly, input_gff3, output_gff3):
    def update_g_features(gene_i):
        g_feature = g_features[gene_i]
        m_feature = g_feature.sub_features[0]
        c_features = [
            x for x in m_feature.sub_features if x.type == 'CDS'
        ]
        c_features_s = sorted(
            c_features, key=lambda x: x.location.start
        )

        e_features = [
            x for x in m_feature.sub_features if x.type == 'exon'
        ]
        e_features_s = sorted(
            e_features, key=lambda x: x.location.start
        )
        cds_start = c_features_s[0].location.start
        cds_end = c_features_s[-1].location.end

        e_features_s[0].location = FeatureLocation(
            cds_start, e_features_s[0].location.end,
            strand=e_features_s[0].location.strand
        )
        e_features_s[-1].location = FeatureLocation(
            e_features_s[-1].location.start, cds_end,
            strand=e_features_s[-1].location.strand
        )
        m_feature.location = FeatureLocation(
            cds_start, cds_end, m_feature.location.strand
        )
        m_feature.sub_features = e_features_s + c_features_s

        g_features[gene_i].location = FeatureLocation(
            cds_start, cds_end,
            strand=g_features[gene_i].location.strand
        )
        g_features[gene_i].sub_features = [m_feature]

    D_fna = SeqIO.to_dict(SeqIO.parse(genome_assembly, "fasta", generic_dna))
    gff_iter = GFF.parse(input_gff3, D_fna)
    my_records = []
    for gff_element in gff_iter:
        g_features = gff_element.features  # Genes in a scaffold
        gene_i = 0
        while gene_i < len(g_features) - 1:
            g_feature = g_features[gene_i]
            g_feature_next = g_features[gene_i + 1]
            end = g_feature.location.end
            start_next = g_feature_next.location.start + 1
            if end >= start_next:
                update_g_features(gene_i)
                update_g_features(gene_i + 1)
            gene_i += 1

        gff_element.features = g_features
        my_records.append(gff_element)

    GFF.write(my_records, open(output_gff3, 'w'))






if __name__ == '__main__':
    main(sys.argv[1:])
