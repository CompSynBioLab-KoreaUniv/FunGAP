#!/usr/bin/env python3

'''
Generate Genbank file using GFF3 and annotations
--Byoungnam Min. Jul 6, 2021
'''

import gzip
import os
import re
from urllib.parse import unquote
from argparse import ArgumentParser
from collections import defaultdict
from datetime import datetime

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def main():
    '''Main function'''
    argparse_usage = (
        'generate_genbank.py -f <input_fna> -g <input_gff3> -a <input_faa>')
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-f', '--input_fna', nargs=1, required=True,
        help='Input FNA file')
    parser.add_argument(
        '-g', '--input_gff3', nargs=1, required=True,
        help='Input GFF3 file')
    parser.add_argument(
        '-a', '--input_faa', nargs=1, required=True,
        help='Input FAA file')
    parser.add_argument(
        '-o', '--output_prefix', nargs='?', default='out',
        help='Output prefix')
    parser.add_argument(
        '-O', '--organism_name', nargs='?', default='organism',
        help='Organism name (default: organism)')
    parser.add_argument(
        '-d', '--data_file_division', nargs='?', default='PLN',
        help='Data file division (default: PLN)')
    parser.add_argument(
        '-t', '--taxonomy', nargs='?', default='Eukaryota',
        help=(
            'Taxonomy separated by "; ", such as "Eukaryota; Fungi"\n'
            '(default: Eukaryota)'))

    args = parser.parse_args()
    input_fna = os.path.abspath(args.input_fna[0])
    input_gff3 = os.path.abspath(args.input_gff3[0])
    input_faa = os.path.abspath(args.input_faa[0])
    output_prefix = os.path.abspath(args.output_prefix)
    organism_name = args.organism_name
    data_file_division = args.data_file_division
    taxonomy = args.taxonomy

    # Run functions :) Slow is as good as Fast
    d_gff3 = parse_gff3(input_gff3)
    generate_genbank(
        input_fna, d_gff3, input_faa, output_prefix, organism_name,
        data_file_division, taxonomy)


# To parse GFF3 I referred the site
# https://techoverflow.net/blog/2013/11/30/parsing-gff3-in-python/
# because I don't think the parser from Biopython is working well
def import_file(input_file):
    '''Import file'''
    with open(input_file) as f_in:
        txt = list(line.rstrip() for line in f_in)
    return txt


def parse_gff_attributes(attribute_string):
    '''Parse the GFF3 attribute column and return a dict'''
    if attribute_string == '.':
        return {}
    ret = {}
    for attribute in attribute_string.split(';'):
        key, value = attribute.split('=')
        ret[unquote(key)] = unquote(value)
    return ret


def parse_gff3(filename):
    '''A minimalistic GFF3 format parser'''
    # Parse with transparent decompression
    open_func = gzip.open if filename.endswith('.gz') else open
    d_gff3 = defaultdict(list)
    with open_func(filename) as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            normalized_info = {
                'seqid': None if parts[0] == '.' else unquote(parts[0]),
                'source':
                    None if parts[1] == '.' else unquote(parts[1]),
                'type': None if parts[2] == '.' else unquote(parts[2]),
                'start': None if parts[3] == '.' else int(parts[3]),
                'end': None if parts[4] == '.' else int(parts[4]),
                'score': None if parts[5] == '.' else float(parts[5]),
                'strand':
                    None if parts[6] == '.' else unquote(parts[6]),
                'phase': None if parts[7] == '.' else unquote(parts[7]),
                'attributes': parse_gff_attributes(parts[8])}
            d_gff3[parts[0]].append(normalized_info)
    return d_gff3


def generate_genbank(
        input_fna, d_gff3, input_faa, output_prefix, organism_name,
        data_file_division, taxonomy):
    '''Generate GenBank format'''
    # Output file name
    outfile = '{}.gb'.format(output_prefix)

    # First, import input_fna in dictionary
    d_fna = SeqIO.to_dict(SeqIO.parse(input_fna, 'fasta'))
    d_faa = SeqIO.to_dict(SeqIO.parse(input_faa, 'fasta'))
    d_fna_sorted = sorted(
        d_fna.items(),
        key=lambda x: int(re.findall(r'\d+', x[0])[0]))

    # Make dictionary for CDS
    d_cds = defaultdict(list)
    d_exon = defaultdict(list)
    for scaffold, records in d_gff3.items():
        for record in records:
            if record['type'] == 'exon':
                exon_parent = record['attributes']['Parent']
                d_exon[exon_parent].append(record)
            elif record['type'] == 'CDS':
                cds_parent = record['attributes']['Parent']
                d_cds[cds_parent].append(record)

    my_seq_records = []
    for scaffold, seq in d_fna_sorted:
        my_seq = Seq(str(seq.seq))
        my_seq_record = SeqRecord(my_seq)
        my_seq_record.id = scaffold
        my_seq_record.description = '{} {}'.format(organism_name, scaffold)
        date = datetime.today().strftime('%d-%^b-%Y')
        my_seq_record.annotations['date'] = date
        my_seq_record.annotations['organism'] = organism_name
        my_seq_record.data_file_division = data_file_division
        my_seq_record.annotations['keywords'] = [
            'Whole genome sequencing project']
        my_seq_record.annotations['taxonomy'] = taxonomy.split('; ')
        my_seq_record.annotations['source'] = organism_name
        my_seq_record.annotations['molecule_type'] = 'DNA'

        # Put source
        source_feature_location = FeatureLocation(0, len(seq))
        source_qualifiers = {
            'organism': organism_name, 'mol_type': 'genomic DNA'}
        source_feature = SeqFeature(
            source_feature_location, type='source',
            qualifiers=source_qualifiers)
        my_seq_record.features.append(source_feature)

        for record in d_gff3[scaffold]:
            my_feature_type = record['type']
            if my_feature_type == ('exon', 'CDS'):
                continue
            # GFFRecord(seqid='contig1', source='AUGUSTUS', type='gene',
            # start=16942, end=19008, score=0.22, strand='+', phase=None,
            # attributes={'Source': 'braker_Y1:g3308.t1', 'ID': 'Triga_00001'})
            my_start = record['start']
            my_end = record['end']
            my_strand = 1 if record['strand'] == '_' else -1

            # Set qualifies for gene
            if my_feature_type == 'gene':
                gene_start = my_start
                gene_end = my_end
                gene_feature_location = FeatureLocation(
                    gene_start, gene_end, strand=my_strand)
                gene_qualifiers = {}
                gene_locus_tag = record['attributes']['ID']
                gene_qualifiers['locus_tag'] = gene_locus_tag
                gene_feature = SeqFeature(
                    gene_feature_location, type=my_feature_type,
                    qualifiers=gene_qualifiers)
                # Append my feature to seq_record
                my_seq_record.features.append(gene_feature)

            elif my_feature_type == 'mRNA':
                sorted_exon_records = sorted(
                    d_exon[record['attributes']['ID']],
                    key=lambda x: x['start'])
                sorted_cds_records = sorted(
                    d_cds[record['attributes']['ID']], key=lambda x: x['start'])
                # Feature locations
                # mRNA location is needed to be modified
                fl_mrna_list = []
                for exon_record in sorted_exon_records:
                    fl_element = FeatureLocation(
                        exon_record['start'], exon_record['end'],
                        strand=my_strand)
                    fl_mrna_list.append(fl_element)

                if len(fl_mrna_list) == 1:
                    mrna_feature_location = fl_mrna_list[0]
                else:
                    mrna_feature_location = CompoundLocation(fl_mrna_list)

                fl_cds_list = []
                for cds_record in sorted_cds_records:
                    fl_element = FeatureLocation(
                        cds_record['start'], cds_record['end'],
                        strand=my_strand)
                    fl_cds_list.append(fl_element)

                # If fl_cds_list is more than 1 use CompoundLocation
                if len(fl_cds_list) == 1:
                    cds_feature_location = fl_cds_list[0]
                else:
                    cds_feature_location = CompoundLocation(fl_cds_list)

                # Qualifier
                mrna_qualifiers = {}
                cds_qualifiers = {}
                mrna_locus_tag = record['attributes']['ID']
                mrna_qualifiers['locus_tag'] = mrna_locus_tag
                if record['score']:
                    mrna_qualifiers['note'] = 'prediction score=%s' % (
                        record['score'])
                cds_qualifiers['locus_tag'] = mrna_locus_tag
                # Get phase
                if my_strand == 1:
                    phase = int(sorted_cds_records[0]['phase']) + 1
                elif my_strand == -1:
                    phase = int(sorted_cds_records[-1]['phase']) + 1
                cds_qualifiers['codon_start'] = phase
                cds_qualifiers['translation'] = str(d_faa[mrna_locus_tag].seq)
                mrna_feature = SeqFeature(
                    mrna_feature_location, type='mRNA',                  qualifiers=mrna_qualifiers)
                cds_feature = SeqFeature(
                    cds_feature_location, type='CDS', qualifiers=cds_qualifiers)
                # Append my feature to seq_record
                my_seq_record.features.append(mrna_feature)
                my_seq_record.features.append(cds_feature)
        my_seq_records.append(my_seq_record)
    SeqIO.write(my_seq_records, outfile, 'genbank')


if __name__ == '__main__':
    main()
