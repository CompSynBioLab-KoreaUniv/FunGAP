#!/usr/bin/env python2

'''
Filter multiple gff3 files based on evidence score

 FunGAP finds "gene blocks" defined as a set of gene models that overlap
with at least one base pair. FunGAP gets all combinations of gene models
in a gene block and calculates the sum of the evidence scores. Gene models
in the block with the highest evidence score are selected as final genes of
that region. Short coding sequence overlap (<10% of coding sequence length)
is allowed.

Input: multiple GFF3 files, Blast score file, Busco score file, Pfam score
       file, bad genes file
Output: filtered gene featrue file in GFF3
'''

# Import modules
import sys
import os
import re
import cPickle
from collections import defaultdict
from argparse import ArgumentParser
import networkx as nx

# Get Logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)
sys.path.append(this_dir)
from set_logging import set_logging

# Parameters
evalue_zero = 2.225074e-308
blast_cutoff = 0.00001  # -1 * log(evalue, 10)


# Main function
def main(argv):
    argparse_usage = (
        'filter_gff3s.py -a <genome_assembly> -i <input_gff3s> '
        '-m <mapping_file> -b <blastp_dict> -B <busco_dict> -p <pfam_dict> '
        '-N <blastn_dict> -g <bad_dict> -n <nr_prot_file> -o <output_dir>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-a", "--genome_assembly", nargs=1, required=True,
        help="Genome assembly file"
    )
    parser.add_argument(
        "-i", "--input_gff3s", nargs='+', required=True,
        help="Multiple gff3 files"
    )
    parser.add_argument(
        "-m", "--mapping_file", nargs=1, required=True,
        help="Mapping txt file (make_nr_prot.py)"
    )
    parser.add_argument(
        "-b", "--blastp_dict", nargs=1, required=True,
        help="Parsed blastp output in dictionary (import_blastp.py)"
    )
    parser.add_argument(
        "-B", "--busco_dict", nargs=1, required=True,
        help="Parsed BUSCO output in dictionary (import_busco.py)"
    )
    parser.add_argument(
        "-p", "--pfam_dict", nargs=1, required=True,
        help="Parsed Pfam_scan output in dictionary (import_pfam.py)"
    )
    parser.add_argument(
        "-N", "--blastn_dict", nargs=1, required=True,
        help="Parsed BLASTn output in dictionary (import_blastn.py)"
    )
    parser.add_argument(
        "-g", "--bad_dict", nargs=1, required=True,
        help="Parsed IPRscan output in dictionary"
    )
    parser.add_argument(
        "-n", "--nr_prot_file", nargs=1, required=True,
        help="nr_prot.faa file (make_nr_prot.py)"
    )
    parser.add_argument(
        "-o", "--output_dir", nargs='?', default='gene_filtering',
        help="Output directory"
    )
    parser.add_argument(
        "-l", "--log_dir", nargs='?', default='log_dir',
        help="Log directory"
    )

    args = parser.parse_args()
    genome_assembly = os.path.abspath(args.genome_assembly[0])
    input_gff3s = [os.path.abspath(x) for x in args.input_gff3s]
    mapping_file = os.path.abspath(args.mapping_file[0])
    blastp_dict = os.path.abspath(args.blastp_dict[0])
    busco_dict = os.path.abspath(args.busco_dict[0])
    pfam_dict = os.path.abspath(args.pfam_dict[0])
    blastn_dict = os.path.abspath(args.blastn_dict[0])
    bad_dict = os.path.abspath(args.bad_dict[0])
    D_bad = cPickle.load(open(bad_dict, 'rb'))
    nr_prot_file = os.path.abspath(args.nr_prot_file[0])
    output_dir = os.path.abspath(args.output_dir)
    log_dir = os.path.abspath(args.log_dir)

    # Create necessary dirs
    create_dir(output_dir, log_dir)

    # Set logging
    log_file = os.path.join(log_dir, 'filter_gff3s.log')
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is as good as Fast
    logger_time.debug('START: Filtering GFF3')
    D_mapping, D_mapping_rev = import_mapping(mapping_file)

    # Import dictionaries
    D_blastp = cPickle.load(open(blastp_dict, 'rb'))
    D_busco = cPickle.load(open(busco_dict, 'rb'))
    D_pfam = cPickle.load(open(pfam_dict, 'rb'))
    D_blastn = cPickle.load(open(blastn_dict, 'rb'))

    # Self-filtering
    for input_gff3 in input_gff3s:
        prefix = re.sub(r'\.gff3$', '', os.path.basename(input_gff3))
        D_gff3, D_gene, D_cds, D_cds_len, D_exon = import_gff3([input_gff3])
        self_filtered = filtering(
            D_cds, D_cds_len, D_blastp, D_busco, D_pfam, D_blastn, D_bad,
            output_dir
        )
        outfile_self = os.path.join(
            output_dir, '{}_filtered.list'.format(prefix)
        )
        outhandle_self = open(outfile_self, 'w')

        cds_len_filtered = 0
        for tup in self_filtered:
            outhandle_self.write('{}\n'.format(tup[1]))
            cds_len_filtered += D_cds_len[tup]
        outhandle_self.close()

    # Filtering
    D_gff3, D_gene, D_cds, D_cds_len, D_exon = import_gff3(input_gff3s)
    final_gene_set = filtering(
        D_cds, D_cds_len, D_blastp, D_busco, D_pfam, D_blastn, D_bad,
        output_dir
    )
    D_prot = import_prot(nr_prot_file, D_mapping_rev)
    write_final_prots(final_gene_set, D_mapping, output_dir)
    write_files(
        genome_assembly, final_gene_set, D_gene, D_gff3, D_prot, D_exon, output_dir,
        D_cds
    )

    cds_len_final = 0
    for tup in final_gene_set:
        cds_len_final += D_cds_len[tup]

    logger_time.debug('DONE : Filtering GFF3')


def import_file(path):
    with open(path) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt


def create_dir(output_dir, log_dir):
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


def import_mapping(mapping_file):
    mapping_txt = import_file(mapping_file)
    D_mapping = {}  # Key: each prediction, value: nr ID
    D_mapping_rev = defaultdict(list)  # Key: nr ID, value: each prediction
    for line in mapping_txt[1:]:
        line_split = line.split('\t')
        prot_name, prefix, prefix_id = line_split
        D_mapping[(prefix, prefix_id)] = prot_name
        D_mapping_rev[prot_name].append((prefix, prefix_id))

    return D_mapping, D_mapping_rev


def import_gff3(gff3_files):
    # Define empty dictionaries
    D_gff3 = defaultdict(list)
    D_exon = defaultdict(list)
    D_gene = {}
    D_cds = defaultdict(lambda: (0, 100000000000, 0))
    D_cds_len = defaultdict(int)

    # Regular expressions
    reg_id = re.compile(r'ID=([^;]+)')
    reg_parent = re.compile(r'Parent=([^;]+)')
    for gff3_file in gff3_files:
        prefix = os.path.basename(gff3_file).split('.')[0]
        gff3_txt = import_file(gff3_file)

        D_parent = {}
        for line in gff3_txt:
            if not re.search('\t', line):
                continue  # Only consider line containing tabs

            line_split = line.split('\t')
            (
                scaffold, source, feat_type, start, end, score,
                strand, phase, attr
            ) = line_split

            m_id = reg_id.search(line)
            if m_id:
                entry_id = m_id.group(1)
            else:
                entry_id = line_split[8]

            m_parent = reg_parent.search(line)
            if m_parent:
                parent_id = m_parent.group(1)
            else:
                parent_id = ''

            if feat_type in ('mRNA', 'transcript'):
                if entry_id != '':
                    mrna_id = entry_id
                    mrna_parent = parent_id
                    if mrna_parent != '':
                        D_parent[mrna_id] = mrna_parent

                new_mrna_id = (prefix, mrna_id)
                D_gene[new_mrna_id] = (
                    scaffold, source, feat_type, int(start), int(end), score,
                    strand, phase
                )

            elif feat_type == 'CDS':
                cds_gene = parent_id  # mRNA ID
                new_cds_gene = (prefix, cds_gene)
                # Measure length of CDSs
                cds_len = (
                    max([int(start), int(end)]) -
                    min([int(start), int(end)]) + 1
                )
                D_cds_len[new_cds_gene] += cds_len

                if max([int(start), int(end)]) > D_cds[new_cds_gene][2]:
                    D_cds[new_cds_gene] = (
                        scaffold, D_cds[new_cds_gene][1],
                        max([int(start), int(end)])
                    )
                if min([int(start), int(end)]) < D_cds[new_cds_gene][1]:
                    D_cds[new_cds_gene] = (
                        scaffold, min([int(start), int(end)]),
                        D_cds[new_cds_gene][2]
                    )
                D_gff3[new_cds_gene].append((
                    scaffold, source, feat_type, int(start), int(end),
                    score, strand, phase
                ))

            elif feat_type == 'exon':
                exon_gene = parent_id # mRNA ID
                new_exon_gene = (prefix, exon_gene)

                D_exon[new_exon_gene].append((
                    scaffold, source, feat_type, int(start), int(end),
                    score, strand, phase
                ))

    return D_gff3, D_gene, D_cds, D_cds_len, D_exon


def filtering(
    D_cds, D_cds_len, D_blastp, D_busco, D_pfam, D_blastn, D_bad, output_dir
):
    # Filter good gene models
    D_cds_filtered = {}
    for gene_tup, value in D_cds.items():
        if D_bad[gene_tup]:
            continue
        D_cds_filtered[gene_tup] = value

    # Sort CDS: three keys used - scaffold, start, and end
    D_cds_sorted = sorted(
        D_cds_filtered.items(), key=lambda x: (x[1][0], x[1][1], x[1][2])
    )

    # Write score table
    outfile_score = os.path.join(output_dir, 'gene_model_scores.txt')
    outhandle_score = open(outfile_score, 'w')
    header_txt = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        'software', 'software_id', 'blast_score', 'busco_score',
        'pfam_score', 'blastn_score', 'score_sum'
    )
    outhandle_score.write(header_txt)
    for tup in D_cds_sorted:
        gene_tup = tup[0]
        software = gene_tup[0]
        software_id = gene_tup[1]
        blast_score = D_blastp[gene_tup]
        busco_score = D_busco[gene_tup]
        pfam_score = D_pfam[gene_tup]
        blastn_score = D_blastn[gene_tup]
        score_sum = sum([blast_score, busco_score, pfam_score, blastn_score])
        row_txt = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            software, software_id, round(blast_score, 1),
            round(busco_score, 1), round(pfam_score, 1), blastn_score,
            round(score_sum, 1)
        )
        outhandle_score.write(row_txt)

    outhandle_score.close()

    # Find chunks
    model_chunks = []  # It will be list of list
    tmp_list = [D_cds_sorted[0][0]]  # Initialize
    i = 1
    max_cds_end = D_cds_sorted[0][1][2]
    while i < len(D_cds_sorted):
        previous_gene_name, previous_tup = D_cds_sorted[i - 1]
        current_gene_name, current_tup = D_cds_sorted[i]

        previous_scaffold = previous_tup[0]
        current_scaffold = current_tup[0]

        current_cds_start = current_tup[1]
        current_cds_end = current_tup[2]

        # If overlapped
        if previous_scaffold != current_scaffold:
            model_chunks.append(tmp_list)
            tmp_list = [current_gene_name]  # Initialize
            max_cds_end = int(current_cds_end)
        elif max_cds_end >= current_cds_start:
            tmp_list.append(current_gene_name)
            max_cds_end = max(max_cds_end, current_cds_end)
        else:
            model_chunks.append(tmp_list)
            tmp_list = [current_gene_name]  # Initialize
            max_cds_end = int(current_cds_end)

        i += 1

    # Filtering
    final_gene_set = []
    i = 1
    for model_chunk in model_chunks:
        i += 1
        # Build Graph
        G = nx.Graph()
        for gene_name1 in model_chunk:
            for gene_name2 in model_chunk:
                G.add_node(gene_name1)
                scaffold1, start1, end1 = D_cds[gene_name1]
                scaffold2, start2, end2 = D_cds[gene_name2]
                if start1 == start2 and end1 == end2:
                    continue
                overlap = min(end1, end2) - max(start1, start2)
                condition1 = overlap < (end1 - start1 + 1) * 0.1
                condition2 = overlap < (end2 - start2 + 1) * 0.1
                if overlap == 0 or condition1 and condition2:
                    G.add_edge(gene_name1, gene_name2)
        all_combs = list(nx.find_cliques(G))

        # Get score and pick best one
        max_score = 0
        final_chunk = []
        for comb in all_combs:
            score = 0
            for element in comb:
                blast_score = D_blastp[element]
                busco_score = D_busco[element]
                pfam_score = D_pfam[element]
                blastn_score = D_blastn[element]

                score += sum(
                    [blast_score, pfam_score, busco_score, blastn_score]
                )

            if score == max_score:
                final_chunk += [comb]
            elif score > max_score:
                final_chunk = [comb]
                max_score = score

        # If all scores are zero, get maximally covered candidates
        if max_score == 0 and len(final_chunk) > 1:
            max_cds_len = 0
            for comb in final_chunk:
                total_cds_len = sum(D_cds_len[x] for x in comb)
                if total_cds_len > max_cds_len:
                    final_chunk2 = comb
                    max_cds_len = total_cds_len
        else:
            final_chunk2 = final_chunk[0]

        # Get score per element
        final_gene_set += final_chunk2

    return final_gene_set


def import_prot(nr_prot, D_mapping_rev):
    D_prot = defaultdict(str)
    prot_txt = import_file(nr_prot)
    for line in prot_txt:
        if re.search('^>', line):
            prot_id = line.split(' ')[0].replace('>', '')
            lst = D_mapping_rev[prot_id]

        else:
            for element in lst:
                prefix, prefix_id = element
                D_prot[(prefix, prefix_id)] += line

    return D_prot


def write_final_prots(final_gene_set, D_mapping, output_dir):
    prot_names = []
    for final_gene in final_gene_set:
        prot_name = D_mapping[final_gene]
        prot_names.append(prot_name)

    # Write to file
    outfile_protnames = os.path.join(output_dir, 'filtered_protnames.list')
    outhandle_protnames = open(outfile_protnames, 'w')
    outhandle_protnames.write('\n'.join(list(set(prot_names))))
    outhandle_protnames.close()


def write_files(
    genome_assembly, final_gene_set, D_gene, D_gff3, D_prot, D_exon, output_dir,
    D_cds
):
    D_scaffold = {}
    scaffold_i = 0
    genome_assembly_txt = import_file(genome_assembly)
    for line in genome_assembly_txt:
        if not line.startswith('>'):
            continue
        scaffold_name = line.split(' ')[0].replace('>', '')
        D_scaffold[scaffold_name] = scaffold_i
        scaffold_i += 1
    final_gene_set_sorted = sorted(
        final_gene_set, key=lambda x: (D_scaffold[D_gene[x][0]], D_cds[x][1])
    )
    output_gff3 = open(os.path.join(output_dir, 'filtered_1.gff3'), 'w')
    output_gff3.write('##gff-version 3\n')  # gff3 Header
    i = 1
    for gene in final_gene_set_sorted:
        (
            gscaffold, gsource, gfeat_type, gstart, gend, gscore,
            gstrand, gphase
        ) = D_gene[gene]
        gid = 'ID={}_{};prediction_source={}:{}'.format(
            'gene', str(i).zfill(5), gene[0], gene[1]
        )
        output_gff3.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            gscaffold, gsource, 'gene', gstart, gend, gscore,
            gstrand, gphase, gid
        ))

        mid = 'ID={}_{}.t1;Parent={}_{}'.format(
            'gene', str(i).zfill(5), 'gene', str(i).zfill(5)
        )
        output_gff3.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            gscaffold, gsource, 'mRNA', gstart, gend, gscore, gstrand,
            gphase, mid
        ))

        if gene in D_exon:
            j = 1
            for exon in D_exon[gene]:
                (
                    escaffold, esource, efeat_type, estart, eend, escore,
                    estrand, ephase
                ) = exon
                eid = 'ID={}_{}.t1.e{};Parent={}_{}.t1'.format(
                    'gene', str(i).zfill(5), j, 'gene', str(i).zfill(5)
                )
                row_txt = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    escaffold, esource, 'exon', estart, eend, escore,
                    estrand, ephase, eid
                )
                output_gff3.write(row_txt)
                j += 1
        else:
            j = 1
            for cds in D_gff3[gene]:
                (
                    cscaffold, csource, cfeat_type, cstart, cend, cscore,
                    cstrand, cphase
                ) = cds
                cid = 'ID={}_{}.t1.e{};Parent={}_{}.t1'.format(
                    'gene', str(i).zfill(5), j, 'gene', str(i).zfill(5)
                )
                output_gff3.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    cscaffold, csource, 'exon', cstart, cend, cscore,
                    cstrand, '.', cid
                ))
                j += 1

        j = 1
        for cds in D_gff3[gene]:
            (
                cscaffold, csource, cfeat_type, cstart, cend, cscore,
                cstrand, cphase
            ) = cds
            cid = 'ID={}_{}.t1.c{};Parent={}_{}.t1'.format(
                'gene', str(i).zfill(5), j, 'gene', str(i).zfill(5)
            )
            output_gff3.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                cscaffold, csource, 'CDS', cstart, cend, cscore,
                cstrand, cphase, cid
            ))
            j += 1

        i += 1

    output_gff3.close()

    # Write protein faa file
    output_prot = open(os.path.join(output_dir, 'filtered_prot.faa'), 'w')
    for gene_num, gene in enumerate(final_gene_set_sorted, start=1):
        output_prot.write('>{}_{}.t1 prediction_source={}:{}\n'.format(
            'gene', str(gene_num).zfill(5), gene[0], gene[1]
        ))
        i = 0
        while i < len(D_prot[gene]):
            output_prot.write('{}\n'.format(D_prot[gene][i:i + 60]))
            i += 60

    output_prot.close()


if __name__ == "__main__":
    main(sys.argv[1:])
