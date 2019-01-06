#!/usr/bin/env python2

'''
Create Markedown document
'''

# Import modules
import re
import os
import sys
import datetime
import subprocess
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from markdown2 import markdown
from Bio.Alphabet import IUPAC
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
from argparse import ArgumentParser
from Bio.Alphabet import generic_dna

# Get Logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)
sys.path.append(this_dir)
from import_config import import_config

# Parameters
D_conf = import_config(this_dir)


# Main function
def main(argv):
    argparse_usage = (
        'create_markdown.py -f <input_fasta> -g <input_gff3> '
        '-t <trinity_assembly> -b <bam_file> -o <output_dir>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-f', '--input_fasta', nargs=1, required=True,
        help='Genome assembly file in FASTA format'
    )
    parser.add_argument(
        '-g', '--input_gff3', nargs=1, required=True,
        help='Input GFF3 file'
    )
    parser.add_argument(
        '-t', '--trinity_assembly', nargs=1, required=True,
        help='Trinity assembly output (FASTA)'
    )
    parser.add_argument(
        '-b', '--bam_file', nargs=1, required=True,
        help='Hisat2 log'
    )
    parser.add_argument(
        '-o', '--output_dir', nargs='?', default='fungap_out',
        help='Output directory'
    )

    args = parser.parse_args()
    input_fasta = os.path.abspath(args.input_fasta[0])
    input_gff3 = os.path.abspath(args.input_gff3[0])
    trinity_assembly = os.path.abspath(args.trinity_assembly[0])
    bam_file = os.path.abspath(args.bam_file[0])
    output_dir = os.path.abspath(args.output_dir)

    # Run functions :) Slow is as good as Fast
    create_dir(output_dir)
    D_fasta = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta', generic_dna))
    D_gff3 = parse_gff3(input_gff3)
    D_cds_coords, protein_lengths, D_stat = get_stats(D_fasta, D_gff3)
    D_stat = get_stats2(D_fasta, D_cds_coords, D_stat)
    D_trinity = get_stats_trinity(trinity_assembly, bam_file)
    trans_len_dist_png = draw_trans_len_dist(D_trinity, output_dir)
    prot_len_dist_png = draw_prot_len_dist(protein_lengths, output_dir)
    create_markdown(
        D_stat, D_trinity, trans_len_dist_png, prot_len_dist_png, output_dir
    )


def import_file(input_file):
    with open(input_file) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt


def get_reverse_complement(nuc_seq):
    my_dna = Seq(nuc_seq, generic_dna)
    rev_comp_dna = str(my_dna.reverse_complement())
    return rev_comp_dna


def create_dir(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


def parse_gff3(input_gff3):
    # Parse gff3
    gff3 = import_file(input_gff3)

    regex_id = re.compile('ID=(\S+?);')
    regex_parent = re.compile(r'Parent=([^;]+)')
    D_prot = {}
    D_gff3 = defaultdict(list)
    for line in gff3:
        if not re.search('\t', line):
            continue

        line_split = line.split('\t')  # Split line with tab
        scaffold = line_split[0]
        entry_type = line_split[2]
        start = int(line_split[3])
        end = int(line_split[4])
        strand = line_split[6]

        if entry_type in ('mRNA', 'transcript'):
            mrna_entry_id = regex_id.search(line).group(1)
            prot_id = regex_parent.search(line).group(1)
            D_prot[mrna_entry_id] = prot_id

        elif entry_type == 'CDS':
            phase = int(line_split[7])
            parent = regex_parent.search(line).group(1)
            prot_id = parent
            D_gff3[prot_id].append((scaffold, start, end, strand, phase))

    # Sort dictionary
    D_gff3_sorted = {}
    for prot_id, feature_list in D_gff3.items():
        D_gff3_sorted[prot_id] = sorted(feature_list, key=lambda x: int(x[1]))

    return D_gff3_sorted


def get_stats(D_fasta, D_gff3):
    # Get stats
    D_stat = {}
    cds_lengths = []
    protein_lengths = []
    exon_lengths = []
    transcript_lengths = []
    intron_lengths = []
    num_introns = []
    num_exons = []
    num_spliced = 0
    single_exon_genes = 0
    total_genes = 0
    D_cds_seq = {}
    D_cds_coords = defaultdict(list)

    sorted_genes = sorted(
        D_gff3.items(), key=lambda x: (
            int(re.findall(r'\d+', x[0])[0]),
            x[1][0][1]
        )
    )

    for prot_id, tuples in sorted_genes:
        total_genes += 1
        tmp_prot_len = 0
        if len(tuples) > 1:
            num_spliced += 1

        cds_seq = ''
        for tup in tuples:
            scaffold, start, end, strand, phase = tup
            if strand == '+' and tup == tuples[0]:
                start = start + phase
            elif strand == '-' and tup == tuples[-1]:
                end = end - phase

            tmp_prot_len += end - start + 1
            exon_lengths.append(end - start + 1)
            # Get sequence
            cds_seq += str(D_fasta[scaffold][start - 1:end].seq)
            # Store in dictionary
            D_cds_coords[scaffold].append((start, end))

        if strand == '-':
            cds_seq = get_reverse_complement(cds_seq)

        D_cds_seq[prot_id] = cds_seq
        cds_length = tmp_prot_len
        cds_lengths.append(cds_length)
        protein_length = tmp_prot_len / 3
        protein_lengths.append(protein_length)
        transcript_length = int(tuples[-1][2]) - int(tuples[0][1]) + 1
        transcript_lengths.append(transcript_length)
        num_intron = len(tuples) - 1
        if num_intron > 0:
            intron_start = [x[2] for x in tuples[:-1]]
            intron_end = [x[1] for x in tuples[1:]]
            intron_length = [
                y - x - 1 for x, y in zip(intron_start, intron_end)
            ]
            intron_lengths += intron_length
            num_introns.append(len(tuples) - 1)
        else:
            intron_median = 0
            num_introns_median = 0
        num_exons.append(len(tuples))
        if len(tuples) == 1:
            single_exon_genes += 1

    intron_median = np.median(np.array(intron_lengths))
    intron_len_average = np.average(np.array(intron_lengths))
    num_introns_median = np.median(np.array(num_introns))
    exon_median = np.median(np.array(exon_lengths))
    exon_len_average = np.average(np.array(exon_lengths))
    cds_average = np.average(cds_lengths)
    cds_median = np.median(cds_lengths)
    protein_average = np.average(np.array(protein_lengths))
    protein_median = np.median(np.array(protein_lengths))
    transcript_median = np.median(np.array(transcript_lengths))
    transcript_average = np.average(np.array(transcript_lengths))
    num_exons_median = np.median(np.array(num_exons))

    # Guitar
    percent_splice = round(float(num_spliced) / total_genes * 100, 2)
    total_bases_lst = [len(str(x.seq)) for x in D_fasta.values()]
    total_bases = sum(total_bases_lst)
    gene_density = float(total_genes) / total_bases
    gene_density = gene_density * 1000000
    gene_density = round(gene_density, 2)

    # Get GC content of CDS seq
    full_cds_seq = ''.join(D_cds_seq.values())
    my_seq = Seq(full_cds_seq, IUPAC.unambiguous_dna)
    cds_gc_percent = GC(my_seq)
    # Percent coding
    coding_percent = float(len(full_cds_seq)) / total_bases
    coding_percent = coding_percent * 100
    coding_percent = round(coding_percent, 2)

    D_stat['Total genes'] = total_genes
    D_stat['Transcript length'] = (
        round(transcript_average, 1), transcript_median
    )
    D_stat['CDS length'] = (round(cds_average, 1), cds_median)
    D_stat['Protein length'] = (round(protein_average, 1), protein_median)
    D_stat['Exon length'] = (round(exon_len_average, 1), exon_median)
    D_stat['Intron length'] = (round(intron_len_average, 1), intron_median)
    D_stat['Spliced'] = (num_spliced, percent_splice)
    D_stat['Gene density'] = gene_density
    D_stat['Num introns'] = sum(num_introns)
    D_stat['Num introns per gene'] = num_introns_median
    D_stat['Num exons'] = sum(num_exons)
    D_stat['Num exons per gene'] = num_exons_median
    D_stat['Num single exon genes'] = single_exon_genes
    D_stat['Percent coding region'] = (len(full_cds_seq), coding_percent)
    D_stat['Coding region GC'] = round(cds_gc_percent, 2)

    return D_cds_coords, protein_lengths, D_stat


def get_stats2(D_fasta, D_cds_coords, D_stat):
    non_coding_seq = ''
    scaffolds_with_gene = []

    # Handle scaffolds without genes
    for scaffold, seq in D_fasta.items():
        if scaffold in scaffolds_with_gene:
            continue
        non_coding_seq += str(D_fasta[scaffold].seq)

    total_bases_lst = [len(str(x.seq)) for x in D_fasta.values()]
    total_bases = sum(total_bases_lst)
    non_coding_percent = float(len(non_coding_seq)) / total_bases
    non_coding_percent = non_coding_percent * 100
    non_coding_percent = round(non_coding_percent, 2)
    my_seq = Seq(non_coding_seq, IUPAC.unambiguous_dna)
    non_coding_gc = GC(my_seq)

    D_stat['Percent non-coding region'] = (non_coding_percent)
    D_stat['Non-coding region GC'] = round(non_coding_gc, 2)

    return D_stat


def get_stats_trinity(trinity_assembly, bam_file):
    trinity_txt = import_file(trinity_assembly)
    D_contig = defaultdict(int)
    for line in trinity_txt:
        if line.startswith('>'):
            contig_name = line.split(' ')[0].replace('>', '')
        else:
            D_contig[contig_name] += len(line)

    num_contigs = len(D_contig)
    total_size = sum(D_contig.values())
    long_contigs = sum(1 for x in D_contig.values() if x > 1000)


    samtools_bin = D_conf['SAMTOOLS_PATH']
    command = '{} view -c {}'.format(samtools_bin, bam_file)
    process1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    output1 = int(process1.communicate()[0])

    D_trinity = {}
    D_trinity['Total contigs'] = num_contigs
    D_trinity['Total size'] = total_size
    D_trinity['Long contigs'] = long_contigs
    D_trinity['Num mapped reads'] = output1
    D_trinity['Length dist'] = D_contig.values()
    return D_trinity


def draw_trans_len_dist(D_trinity, output_dir):
    trans_lengths = D_trinity['Length dist']

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(
        trans_lengths, facecolor='#fdc50c', alpha=1,
        bins=150
    )
    plt.title('Transcript length distribution')
    plt.xlabel('Transcript length (nt)')
    plt.ylabel('Frequency')
    ax.set_xlim(0, 5000)
    outpng = os.path.join(output_dir, 'fungap_out_trans_len_dist.png')
    plt.savefig(
        outpng, dpi=500, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1
    )
    return outpng


def draw_prot_len_dist(protein_lengths, output_dir):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(
        protein_lengths, facecolor='#4b85c5', alpha=1,
        bins=150
    )
    plt.title('Protein length distribution')
    plt.xlabel('Amino acids (aa)')
    plt.ylabel('Frequency')
    ax.set_xlim(0, 2000)
    outpng = os.path.join(output_dir, 'fungap_out_prot_len_dist.png')
    plt.savefig(
        outpng, dpi=500, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1
    )
    return outpng


def create_markdown(
    D_stat, D_trinity, trans_len_dist_png, prot_len_dist_png, output_dir
):
    # Header
    header_txt = '# FunGAP report'
    md = markdown(header_txt)

    # Date
    date_txt = '_Created at {}_'.format(datetime.date.today())
    md += markdown(date_txt)

    # Number of genes
    num_genes_txt = 'The **{}** genes were predicted in this genome.'.format(
        '{:,}'.format(D_stat['Total genes'])
    )
    md += markdown(num_genes_txt)

    # Gene structure summary
    gene_structure_txt = '### 1. Gene structure'
    md += markdown(gene_structure_txt)

    gene_structure_table = '''
|| *Attributes* || *Values* ||
|| Total protein-coding genes || {} ||
|| Transcript length (avg / med) || {} / {} ||
|| CDS length (avg / med) || {} / {} ||
|| Protein length (avg / med) || {} / {} ||
|| Exon length (avg / med) || {} / {} ||
|| Intron length (avg / med) || {} / {} ||
|| Spliced genes || {} ({}%) ||
|| Gene density (genes/Mb) || {} ||
|| Number of introns || {} ||
|| Number of introns per gene (med) || {} ||
|| Number of exons || {} ||
|| Number of exons per gene (med) || {} ||
    '''.format(
        '{:,}'.format(D_stat['Total genes']),
        '{:,}'.format(D_stat['Transcript length'][0]),
        '{:,}'.format(D_stat['Transcript length'][1]),
        '{:,}'.format(D_stat['CDS length'][0]),
        '{:,}'.format(D_stat['CDS length'][1]),
        '{:,}'.format(D_stat['Protein length'][0]),
        '{:,}'.format(D_stat['Protein length'][1]),
        '{:,}'.format(D_stat['Exon length'][0]),
        '{:,}'.format(D_stat['Exon length'][1]),
        '{:,}'.format(D_stat['Intron length'][0]),
        '{:,}'.format(D_stat['Intron length'][1]),
        '{:,}'.format(D_stat['Spliced'][0]),
        D_stat['Spliced'][1],
        '{:,}'.format(D_stat['Gene density']),
        '{:,}'.format(D_stat['Num introns']),
        D_stat['Num introns per gene'],
        '{:,}'.format(D_stat['Num exons']),
        D_stat['Num exons per gene'],

    )
    md += markdown(gene_structure_table, extras=['wiki-tables'])

    # Transcript assembly summary
    transcript_header = '### 2. Transcriptome reads assembly'
    md += '<br>'
    md += markdown(transcript_header)

    transcript_stats_table = '''
|| *Attributes* || *Values* ||
|| Number of mapped reads || {} ||
|| Number of assembled contigs || {} ||
|| Number of contigs > 1 kbp || {} ||
|| Total transcript size (Mbp) || {} ||
'''.format(
        '{:,}'.format(D_trinity['Num mapped reads']),
        '{:,}'.format(D_trinity['Total contigs']),
        '{:,}'.format(D_trinity['Long contigs']),
        '{:,}'.format(D_trinity['Total size'])
    )
    md += markdown(transcript_stats_table, extras=['wiki-tables'])

    # Transscript length distribution
    trans_len_txt = '### 3. Transcript length distribution'
    md += '<br>'
    md += markdown(trans_len_txt)
    md += markdown('![Transcript length distribution]({})'.format(
        os.path.basename(trans_len_dist_png)
    ))

    # Protein length distribution
    prot_len_txt = '### 4. Protein length distribution'
    md += '<br>'
    md += markdown(prot_len_txt)
    md += markdown('![Protein length distribution]({})'.format(
        os.path.basename(prot_len_dist_png))
    )

    outfile = os.path.join(output_dir, 'fungap_out.html')

    # header including css
    header_txt = '''
<head>
<style>
body {
    font-family: sans-serif;
    padding: 50px 30px 50px 80px;
}
img {width: 500;}
td {
    border-bottom: 1px solid #ddd;
    padding: 8px;
}
</style>
</head>
<body>
'''
    outhandle = open(outfile, 'w')
    outhandle.write(header_txt)
    outhandle.write(md)

    footer_txt = '''
</body>'''
    outhandle.write(footer_txt)
    outhandle.close()


if __name__ == '__main__':
    main(sys.argv[1:])
