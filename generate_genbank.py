#!/usr/bin/python

'''
Generate Genbank file using GFF3 and annotations
Author Byoungnam Min on Dec 26, 2015
'''

# Import modeuls
from __future__ import with_statement
import os
import sys
import re
import gzip
import urllib
from datetime import datetime
from argparse import ArgumentParser
from collections import namedtuple, defaultdict

from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation, CompoundLocation

# Initialized values
gffInfoFields = [
    "seqid", "source", "type", "start", "end", "score", "strand",
    "phase", "attributes"
]
GFFRecord = namedtuple("GFFRecord", gffInfoFields)


def main(argv):
    argparse_usage = (
        'generate_genbank.py -f <input_fna> -g <input_gff3> '
        '-a <input_faa> -o <output_prefix>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-f", "--input_fna", dest="input_fna", nargs=1,
        help="Input FNA file"
    )
    parser.add_argument(
        "-g", "--input_gff3", dest="input_gff3", nargs=1,
        help="Input GFF3 file"
    )
    parser.add_argument(
        "-a", "--input_faa", dest="input_faa", nargs=1,
        help="Input FAA file"
    )
    parser.add_argument(
        "-o", "--output_prefix", dest="output_prefix", nargs=1,
        help="Output prefix"
    )
    parser.add_argument(
        "-O", "--organism_name", dest="organism_name", nargs='?',
        help="Organism name (default: organism)"
    )
    parser.add_argument(
        "-d", "--data_file_division", dest="data_file_division", nargs='?',
        help="Data file division (default: PLN)"
    )
    parser.add_argument(
        "-t", "--taxonomy", dest="taxonomy", nargs='?',
        help=(
            "Taxonomy separated by '; ', such as 'Eukaryota; Fungi'\n"
            "(default: Eukaryota)"
        )
    )

    args = parser.parse_args()
    if args.input_fna:
        input_fna = os.path.abspath(args.input_fna[0])
    else:
        print '[ERROR] Please provide INPUT FNA FILE'
        parser.print_help()
        sys.exit(2)

    if args.input_gff3:
        input_gff3 = os.path.abspath(args.input_gff3[0])
    else:
        print '[ERROR] Please provide INPUT GFF3 FILE'
        parser.print_help()
        sys.exit(2)

    if args.input_faa:
        input_faa = os.path.abspath(args.input_faa[0])
    else:
        print '[ERROR] Please provide INPUT FAA FILE'
        parser.print_help()
        sys.exit(2)

    if args.output_prefix:
        output_prefix = os.path.abspath(args.output_prefix[0])
    else:
        print '[ERROR] Please provide OUTPUT PREFIX'
        parser.print_help()
        sys.exit(2)

    if args.organism_name:
        organism_name = args.organism_name
    else:
        organism_name = 'organism'

    if args.data_file_division:
        data_file_division = args.data_file_division
    else:
        data_file_division = 'PLN'

    if args.taxonomy:
        taxonomy = args.taxonomy
    else:
        taxonomy = 'Eukaryota'

    # Run functions :) Slow is as good as Fast
    generate_genbank(
        input_fna, input_gff3, input_faa, output_prefix,
        organism_name, data_file_division, taxonomy
    )


# To parse GFF3 I referred the site
# https://techoverflow.net/blog/2013/11/30/parsing-gff3-in-python/
# because I don't think the parser from Biopython is working well
def import_file(input_file):
    with open(input_file) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt


def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dict"""
    if attributeString == ".":
        return {}
    ret = {}
    for attribute in attributeString.split(";"):
        key, value = attribute.split("=")
        ret[urllib.unquote(key)] = urllib.unquote(value)
    return ret


def parseGFF3(filename):
    """
    A minimalistic GFF3 format parser.
    Yields objects that contain info about a single GFF3 feature.

    Supports transparent gzip decompression.
    """
    # Parse with transparent decompression
    openFunc = gzip.open if filename.endswith(".gz") else open
    with openFunc(filename) as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            # If this fails, the file format is not standard-compatible
            assert len(parts) == len(gffInfoFields)
            # Normalize data
            normalizedInfo = {
                "seqid": None if parts[0] == "." else urllib.unquote(parts[0]),
                "source":
                    None if parts[1] == "." else urllib.unquote(parts[1]),
                "type": None if parts[2] == "." else urllib.unquote(parts[2]),
                "start": None if parts[3] == "." else int(parts[3]),
                "end": None if parts[4] == "." else int(parts[4]),
                "score": None if parts[5] == "." else float(parts[5]),
                "strand":
                    None if parts[6] == "." else urllib.unquote(parts[6]),
                "phase": None if parts[7] == "." else urllib.unquote(parts[7]),
                "attributes": parseGFFAttributes(parts[8])
            }
            # Alternatively, you can emit the dictionary here,
            # if you need mutability:
            # yield normalizedInfo
            yield GFFRecord(**normalizedInfo)


def generate_genbank(
    input_fna, input_gff3, input_faa, output_prefix,
    organism_name, data_file_division, taxonomy
):
    # Output file name
    outfile = '%s.gb' % (output_prefix)

    # First, import input_fna in dictionary
    D_fna = SeqIO.to_dict(SeqIO.parse(input_fna, "fasta", generic_dna))
    D_faa = SeqIO.to_dict(SeqIO.parse(input_faa, "fasta", generic_protein))

    D_fna_sorted = sorted(
        D_fna.items(),
        key=lambda x: int(re.findall(r'\d+', x[0])[0])
    )

    # Make dictionary for CDS
    D_cds = defaultdict(list)
    D_exon = defaultdict(list)
    for record in parseGFF3(input_gff3):
        if record.type == 'exon':
            exon_parent = record.attributes['Parent']
            D_exon[exon_parent].append(record)

        elif record.type == 'CDS':
            cds_parent = record.attributes['Parent']
            D_cds[cds_parent].append(record)

    my_seq_records = []
    for scaffold, seq in D_fna_sorted:
        my_seq = Seq(str(seq.seq))
        my_seq_record = SeqRecord(my_seq)
        my_seq_record.seq.alphabet = generic_dna

        my_seq_record.description = '{} {}'.format(organism_name, scaffold)
        date = datetime.today().strftime('%d-%^b-%Y')
        my_seq_record.annotations["date"] = date
        my_seq_record.annotations["organism"] = organism_name
        print data_file_division
        my_seq_record.data_file_division = data_file_division
        quit()
        my_seq_record.annotations["keywords"] = [
            "Whole genome sequencing project"
        ]
        my_seq_record.annotations["taxonomy"] = taxonomy.split('; ')
        my_seq_record.annotations["source"] = organism_name

        for record in parseGFF3(input_gff3):
            if scaffold != record.seqid:
                continue

            my_feature_type = record.type
            if my_feature_type == ('exon', 'CDS'):
                continue

            # GFFRecord(seqid='contig1', source='AUGUSTUS', type='gene',
            # start=16942, end=19008, score=0.22, strand='+', phase=None,
            # attributes={'Source': 'braker_Y1:g3308.t1', 'ID': 'Triga_00001'})

            my_start = record.start
            my_end = record.end
            if record.strand == '+':
                my_strand = 1
            else:
                my_strand = -1

            # Set qualifies for gene
            if my_feature_type == 'gene':
                gene_start = my_start
                gene_end = my_end

                gene_feature_location = FeatureLocation(
                    gene_start, gene_end, strand=my_strand
                )

                gene_qualifiers = {}
                gene_locus_tag = record.attributes['ID']
                gene_qualifiers['locus_tag'] = gene_locus_tag

                gene_feature = SeqFeature(
                    gene_feature_location, type=my_feature_type,
                    qualifiers=gene_qualifiers
                )

                # Append my feature to seq_record
                my_seq_record.features.append(gene_feature)

            elif my_feature_type == 'mRNA':
                sorted_exon_records = sorted(
                    D_exon[record.attributes['ID']], key=lambda x: x.start
                )
                sorted_cds_records = sorted(
                    D_cds[record.attributes['ID']], key=lambda x: x.start
                )

                # Feature locations
                # mRNA location is needed to be modified
                fl_mrna_list = []
                for exon_record in sorted_exon_records:
                    fl_element = FeatureLocation(
                        exon_record.start, exon_record.end, strand=my_strand
                    )
                    fl_mrna_list.append(fl_element)

                if len(fl_mrna_list) == 1:
                    mrna_feature_location = fl_mrna_list[0]
                else:
                    mrna_feature_location = CompoundLocation(fl_mrna_list)

                fl_cds_list = []
                for cds_record in sorted_cds_records:
                    fl_element = FeatureLocation(
                        cds_record.start, cds_record.end, strand=my_strand
                    )
                    fl_cds_list.append(fl_element)

                # If fl_cds_list is more than 1 use CompoundLocation
                if len(fl_cds_list) == 1:
                    cds_feature_location = fl_cds_list[0]
                else:
                    cds_feature_location = CompoundLocation(fl_cds_list)

                # Qualifier
                mrna_qualifiers = {}
                cds_qualifiers = {}

                mrna_locus_tag = record.attributes['ID']
                mrna_qualifiers['locus_tag'] = mrna_locus_tag
                if record.score:
                    mrna_qualifiers['note'] = "prediction score=%s" % (
                        record.score
                    )

                cds_qualifiers['locus_tag'] = mrna_locus_tag
                # Get phase
                if my_strand == 1:
                    phase = int(sorted_cds_records[0].phase) + 1
                elif my_strand == -1:
                    phase = int(sorted_cds_records[-1].phase) + 1
                cds_qualifiers['codon_start'] = phase
                cds_qualifiers['translation'] = str(D_faa[mrna_locus_tag].seq)

                mrna_feature = SeqFeature(
                    mrna_feature_location, type='mRNA',
                    qualifiers=mrna_qualifiers
                )

                cds_feature = SeqFeature(
                    cds_feature_location, type='CDS',
                    qualifiers=cds_qualifiers
                )
                # Append my feature to seq_record
                my_seq_record.features.append(mrna_feature)
                my_seq_record.features.append(cds_feature)
        my_seq_records.append(my_seq_record)
        break

    SeqIO.write(my_seq_records, outfile, "genbank")
    sys.exit()


if __name__ == "__main__":
    main(sys.argv[1:])
