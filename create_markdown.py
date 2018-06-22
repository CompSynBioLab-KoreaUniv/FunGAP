#!/usr/bin/python

'''
Create Markedown document
'''

# Import modules
import re
import os
import sys
import datetime
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


# Main function
def main(argv):
    argparse_usage = (
        'create_markdown.py -f <input_fasta> -g <input_gff3> '
        '-t <trinity_assembly> -H <hisat2_log> -o <output_prefix>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-f", "--input_fasta", dest="input_fasta", nargs=1,
        help="Genome assembly file in FASTA format"
    )
    parser.add_argument(
        "-g", "--input_gff3", dest="input_gff3", nargs=1,
        help="Input GFF3 file"
    )
    parser.add_argument(
        "-t", "--trinity_assembly", dest="trinity_assembly", nargs=1,
        help="Trinity assembly output (FASTA)"
    )
    parser.add_argument(
        "-H", "--hisat2_log", dest="hisat2_log", nargs=1,
        help="Hisat2 log"
    )
    parser.add_argument(
        "-o", "--output_prefix", dest="output_prefix", nargs=1,
        help="Output prefix"
    )

    args = parser.parse_args()
    if args.input_fasta:
        input_fasta = os.path.abspath(args.input_fasta[0])
    else:
        print '[ERROR] Please provide GENOME ASSEMBLY FILE'
        sys.exit(2)

    if args.input_gff3:
        input_gff3 = os.path.abspath(args.input_gff3[0])
    else:
        print '[ERROR] Please provide INPUT GFF3'
        sys.exit(2)

    if args.trinity_assembly:
        trinity_assembly = os.path.abspath(args.trinity_assembly[0])
    else:
        print '[ERROR] Please provide TRINITY ASSEMBLY'
        sys.exit(2)

    if args.hisat2_log:
        hisat2_log = os.path.abspath(args.hisat2_log[0])
    else:
        print '[ERROR] Please provide HISAT2 LOG'
        sys.exit(2)

    if args.output_prefix:
        output_prefix = os.path.abspath(args.output_prefix[0])
    else:
        print '[ERROR] Please provide OUTPUT PREFIX'
        sys.exit(2)

    # Run functions :) Slow is as good as Fast
    D_fasta = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta", generic_dna))
    D_gff3 = parse_gff3(input_gff3)
    D_cds_coords, protein_lengths, D_stat = get_stats(D_fasta, D_gff3)
    D_stat = get_stats2(D_fasta, D_cds_coords, D_stat)
    D_trinity = get_stats_trinity(trinity_assembly, hisat2_log)
    trans_len_dist_jpg = draw_trans_len_dist(D_trinity, output_prefix)
    prot_len_dist_jpg = draw_prot_len_dist(protein_lengths, output_prefix)
    create_markdown(
        D_stat, D_trinity, trans_len_dist_jpg, prot_len_dist_jpg, output_prefix
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


def get_stats_trinity(trinity_assembly, hisat2_log):
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

    if os.path.exists(hisat2_log):
        hisat2_log_txt = import_file(hisat2_log)
        reg_reads = re.compile(r'(\d+) reads; of these:')
        reg_nopairs = re.compile(
            r'(\d+) pairs aligned 0 times concordantly or discordantly; of these:'
        )
        reg_aligned_one = re.compile(
            r'        (\d+) \(\S+%\) aligned exactly 1 time'
        )
        reg_aligned_multi = re.compile(
            r'        (\d+) \(\S+%\) aligned >1 times'
        )
        for line in hisat2_log_txt:
            m_reads = reg_reads.search(line)
            if m_reads:
                num_reads = int(m_reads.group(1)) * 2
            m_nopairs = reg_nopairs.search(line)
            if m_nopairs:
                num_nopairs = int(m_nopairs.group(1))
            m_aligned_one = reg_aligned_one.search(line)
            if m_aligned_one:
                num_aligned_one = int(m_aligned_one.group(1))
            m_aligned_multi = reg_aligned_multi.search(line)
            if m_aligned_multi:
                num_aligned_multi = int(m_aligned_multi.group(1))
        num_mapped = (
            num_reads - num_nopairs * 2 + num_aligned_one + num_aligned_multi
        )
        alignment_rate = round(float(num_mapped) / num_reads * 100, 1)
    else:
        num_reads = 'N/A'
        num_mapped = 'N/A'
        alignment_rate = 'N/A'

    D_trinity = {}
    D_trinity['Total contigs'] = num_contigs
    D_trinity['Total size'] = total_size
    D_trinity['Long contigs'] = long_contigs
    D_trinity['Num reads'] = num_reads
    D_trinity['Num mapped reads'] = num_mapped
    D_trinity['Alignment rate'] = alignment_rate
    D_trinity['Length dist'] = D_contig.values()
    return D_trinity


def draw_trans_len_dist(D_trinity, output_prefix):
    trans_lengths = D_trinity['Length dist']

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(
        trans_lengths, facecolor='#fdc50c', alpha=1,
        bins=150
    )
    plt.title("Transcript length distribution")
    plt.xlabel("Transcript length (nt)")
    plt.ylabel("Frequency")
    ax.set_xlim(0, 5000)
    outjpg = '%s_trans_len_dist.jpg' % (output_prefix)
    plt.savefig(
        outjpg, dpi=500, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1
    )
    return outjpg


def draw_prot_len_dist(protein_lengths, output_prefix):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(
        protein_lengths, facecolor='#4b85c5', alpha=1,
        bins=150
    )
    plt.title("Protein length distribution")
    plt.xlabel("Amino acids (aa)")
    plt.ylabel("Frequency")
    ax.set_xlim(0, 2000)
    outjpg = '%s_prot_len_dist.jpg' % (output_prefix)
    plt.savefig(
        outjpg, dpi=500, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1
    )
    return outjpg


def create_markdown(
    D_stat, D_trinity, trans_len_dist_jpg, prot_len_dist_jpg, output_prefix
):
    # Header
    header_txt = '# FunGAP report'
    md = markdown(header_txt)

    # Date
    date_txt = '_Created at %s_' % (datetime.date.today())
    md += markdown(date_txt)

    # Number of genes
    num_genes_txt = 'The **%s** genes were predicted in this genome.' % (
        "{:,}".format(D_stat['Total genes'])
    )
    md += markdown(num_genes_txt)

    # Gene structure summary
    gene_structure_txt = '### 1. Gene structure'
    md += markdown(gene_structure_txt)

    gene_structure_table = '''
|| *Attributes* || *Values* ||
|| Total protein-coding genes || %s ||
|| Transcript length (avg / med) || %s / %s ||
|| CDS length (avg / med) || %s / %s ||
|| Protein length (avg / med) || %s / %s ||
|| Exon length (avg / med) || %s / %s ||
|| Intron length (avg / med) || %s / %s ||
|| Spliced genes || %s (%s%%) ||
|| Gene density (genes/Mb) || %s ||
|| Number of introns || %s ||
|| Number of introns per gene (med) || %s ||
|| Number of exons || %s ||
|| Number of exons per gene (med) || %s ||
    ''' % (
        "{:,}".format(D_stat['Total genes']),
        "{:,}".format(D_stat['Transcript length'][0]),
        "{:,}".format(D_stat['Transcript length'][1]),
        "{:,}".format(D_stat['CDS length'][0]),
        "{:,}".format(D_stat['CDS length'][1]),
        "{:,}".format(D_stat['Protein length'][0]),
        "{:,}".format(D_stat['Protein length'][1]),
        "{:,}".format(D_stat['Exon length'][0]),
        "{:,}".format(D_stat['Exon length'][1]),
        "{:,}".format(D_stat['Intron length'][0]),
        "{:,}".format(D_stat['Intron length'][1]),
        "{:,}".format(D_stat['Spliced'][0]),
        D_stat['Spliced'][1],
        "{:,}".format(D_stat['Gene density']),
        "{:,}".format(D_stat['Num introns']),
        D_stat['Num introns per gene'],
        "{:,}".format(D_stat['Num exons']),
        D_stat['Num exons per gene'],

    )
    md += markdown(gene_structure_table, extras=["wiki-tables"])

    # Transcript assembly summary
    transcript_header = '### 2. Transcriptome reads assembly'
    md += '<br>'
    md += markdown(transcript_header)

    transcript_stats_table = '''
|| *Attributes* || *Values* ||
|| Number of reads || %s ||
|| Number of mapped reads || %s (%s%%) ||
|| Number of assembled contigs || %s ||
|| Number of contigs > 1 kbp || %s ||
|| Total transcript size (Mbp) || %s ||
''' % (
        "{:,}".format(D_trinity['Num reads']),
        "{:,}".format(D_trinity['Num mapped reads']),
        "{:,}".format(D_trinity['Alignment rate']),
        "{:,}".format(D_trinity['Total contigs']),
        "{:,}".format(D_trinity['Long contigs']),
        "{:,}".format(D_trinity['Total size'])
    )
    md += markdown(transcript_stats_table, extras=["wiki-tables"])

    # Transscript length distribution
    trans_len_txt = '### 3. Transcript length distribution'
    md += '<br>'
    md += markdown(trans_len_txt)
    md += markdown("![Transcript length distribution](%s)" % (
        os.path.basename(trans_len_dist_jpg)
    ))

    # Protein length distribution
    prot_len_txt = '### 4. Protein length distribution'
    md += '<br>'
    md += markdown(prot_len_txt)
    md += markdown("![Protein length distribution](%s)" % (
        os.path.basename(prot_len_dist_jpg))
    )

    outfile = '%s.html' % (output_prefix)

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


if __name__ == "__main__":
    main(sys.argv[1:])
