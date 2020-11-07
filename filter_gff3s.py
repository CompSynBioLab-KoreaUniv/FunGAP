#!/usr/bin/env python3

'''
Filter multiple gff3 files based on evidence score

 FunGAP finds 'gene blocks' defined as a set of gene models that overlap
with at least one base pair. FunGAP gets all combinations of gene models
in a gene block and calculates the sum of the evidence scores. Gene models
in the block with the highest evidence score are selected as final genes of
that region. Short coding sequence overlap (<10% of coding sequence length)
is allowed.

Input: multiple GFF3 files, Blast score file, Busco score file, Pfam score
       file, bad genes file
Output: filtered gene featrue file in GFF3
Last updated: Jul 13, 2020
'''

import bisect
import os
import pickle
import re
from argparse import ArgumentParser
from collections import defaultdict

from set_logging import set_logging

# Parameters
EVALUE_ZERO = 2.225074e-308
BLAST_CUTOFF = 0.00001  # -1 * log(evalue, 10)


def main():
    '''Main function'''
    argparse_usage = (
        'filter_gff3s.py -a <genome_assembly> -i <input_gff3s> '
        '-m <mapping_file> -b <blastp_dict> -B <busco_dict> -p <pfam_dict> '
        '-N <blastn_dict> -g <bad_dict> -n <nr_prot_file> -o <output_dir>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-a', '--genome_assembly', nargs=1, required=True,
        help='Genome assembly file'
    )
    parser.add_argument(
        '-i', '--input_gff3s', nargs='+', required=True,
        help='Multiple gff3 files'
    )
    parser.add_argument(
        '-m', '--mapping_file', nargs=1, required=True,
        help='Mapping txt file (make_nr_prot.py)'
    )
    parser.add_argument(
        '-b', '--blastp_dict', nargs=1, required=True,
        help='Parsed blastp output in dictionary (import_blastp.py)'
    )
    parser.add_argument(
        '-B', '--busco_dict', nargs=1, required=True,
        help='Parsed BUSCO output in dictionary (import_busco.py)'
    )
    parser.add_argument(
        '-p', '--pfam_dict', nargs=1, required=True,
        help='Parsed Pfam_scan output in dictionary (import_pfam.py)'
    )
    parser.add_argument(
        '-N', '--blastn_dict', nargs=1, required=True,
        help='Parsed BLASTn output in dictionary (import_blastn.py)'
    )
    parser.add_argument(
        '-g', '--bad_dict', nargs=1, required=True,
        help='Parsed IPRscan output in dictionary'
    )
    parser.add_argument(
        '-n', '--nr_prot_file', nargs=1, required=True,
        help='nr_prot.faa file (make_nr_prot.py)'
    )
    parser.add_argument(
        '-o', '--output_dir', nargs='?', default='gene_filtering',
        help='Output directory'
    )
    parser.add_argument(
        '-l', '--log_dir', nargs='?', default='log_dir',
        help='Log directory'
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
    d_bad = pickle.load(open(bad_dict, 'rb'))
    nr_prot_file = os.path.abspath(args.nr_prot_file[0])
    output_dir = os.path.abspath(args.output_dir)
    log_dir = os.path.abspath(args.log_dir)

    # Create necessary dirs
    create_dir(output_dir, log_dir)

    # Set logging
    log_file = os.path.join(log_dir, 'filter_gff3s.log')
    logger_time = set_logging(log_file)[0]

    # Run functions :) Slow is as good as Fast
    logger_time.debug('START: Filtering GFF3')
    d_mapping, d_mapping_rev = import_mapping(mapping_file)

    # Import dictionaries
    d_blastp = pickle.load(open(blastp_dict, 'rb'))
    d_busco = pickle.load(open(busco_dict, 'rb'))
    d_pfam = pickle.load(open(pfam_dict, 'rb'))
    d_blastn = pickle.load(open(blastn_dict, 'rb'))

    # Self-filtering (for score purpose)
    for input_gff3 in input_gff3s:
        prefix = re.sub(r'\.gff3$', '', os.path.basename(input_gff3))
        d_gff3, d_gene, d_cds, d_cds_len, d_exon = import_gff3([input_gff3])
        d_score = cal_score(
            d_cds, d_blastp, d_busco, d_pfam, d_blastn, output_dir)
        self_filtered = filtering(d_cds, d_cds_len, d_score)
        outfile_self = os.path.join(
            output_dir, '{}_filtered.list'.format(prefix)
        )
        outhandle_self = open(outfile_self, 'w')

        cds_len_filtered = 0
        for tup in self_filtered:
            outhandle_self.write('{}\n'.format(tup[1]))
            cds_len_filtered += d_cds_len[tup]
        outhandle_self.close()

    # Filtering
    d_gff3, d_gene, d_cds, d_cds_len, d_exon = import_gff3(input_gff3s)
    remove_bad_genes(d_cds, d_bad)
    d_score = cal_score(d_cds, d_blastp, d_busco, d_pfam, d_blastn, output_dir)
    final_gene_set = filtering(d_cds, d_cds_len, d_score)
    d_prot = import_prot(nr_prot_file, d_mapping_rev)
    write_final_prots(final_gene_set, d_mapping, output_dir)
    write_files(
        genome_assembly, final_gene_set, d_gene, d_gff3, d_prot, d_exon,
        output_dir, d_cds
    )

    cds_len_final = 0
    for tup in final_gene_set:
        cds_len_final += d_cds_len[tup]

    logger_time.debug('DONE : Filtering GFF3')


def import_file(path):
    '''Import file'''
    with open(path) as f_in:
        txt = list(line.rstrip() for line in f_in)
    return txt


def create_dir(output_dir, log_dir):
    '''Create directory'''
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


def import_mapping(mapping_file):
    '''Import mapping'''
    mapping_txt = import_file(mapping_file)
    d_mapping = {}  # Key: each prediction, value: nr ID
    d_mapping_rev = defaultdict(list)  # Key: nr ID, value: each prediction
    for line in mapping_txt[1:]:
        line_split = line.split('\t')
        prot_name, prefix, prefix_id = line_split
        d_mapping[(prefix, prefix_id)] = prot_name
        d_mapping_rev[prot_name].append((prefix, prefix_id))
    return d_mapping, d_mapping_rev


def import_gff3(gff3_files):
    '''Import GFF3'''
    # Define empty dictionaries
    d_gff3 = defaultdict(list)
    d_exon = defaultdict(list)
    d_gene = {}
    d_cds = defaultdict(lambda: (0, 100000000000, 0))
    d_cds_len = defaultdict(int)

    # Regular expressions
    reg_id = re.compile(r'ID=([^;]+)')
    reg_parent = re.compile(r'Parent=([^;]+)')
    for gff3_file in gff3_files:
        prefix = re.sub(r'\.gff3$', '', os.path.basename(gff3_file))
        gff3_txt = import_file(gff3_file)
        d_parent = {}
        for line in gff3_txt:
            if not re.search('\t', line):
                continue  # Only consider line containing tabs
            line_split = line.split('\t')
            scaffold, source, feat_type, start, end, score, strand, phase = (
                line_split[:8]
            )
            m_id = reg_id.search(line)
            entry_id = m_id.group(1) if m_id else line_split[8]
            m_parent = reg_parent.search(line)
            parent_id = m_parent.group(1) if m_parent else ''
            if feat_type in ('mRNA', 'transcript'):
                if entry_id != '':
                    mrna_id = entry_id
                    mrna_parent = parent_id
                    if mrna_parent != '':
                        d_parent[mrna_id] = mrna_parent

                new_mrna_id = (prefix, mrna_id)
                d_gene[new_mrna_id] = (
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
                d_cds_len[new_cds_gene] += cds_len

                if max([int(start), int(end)]) > d_cds[new_cds_gene][2]:
                    d_cds[new_cds_gene] = (
                        scaffold, d_cds[new_cds_gene][1],
                        max([int(start), int(end)])
                    )
                if min([int(start), int(end)]) < d_cds[new_cds_gene][1]:
                    d_cds[new_cds_gene] = (
                        scaffold, min([int(start), int(end)]),
                        d_cds[new_cds_gene][2]
                    )
                d_gff3[new_cds_gene].append((
                    scaffold, source, feat_type, int(start), int(end),
                    score, strand, phase
                ))
            elif feat_type == 'exon':
                exon_gene = parent_id # mRNA ID
                new_exon_gene = (prefix, exon_gene)

                d_exon[new_exon_gene].append((
                    scaffold, source, feat_type, int(start), int(end),
                    score, strand, phase
                ))
    return d_gff3, d_gene, d_cds, d_cds_len, d_exon


def remove_bad_genes(d_cds, d_bad):
    '''Remove bad genes'''
    for gene_model in d_bad.keys():
        del d_cds[gene_model]


def cal_score(d_cds, d_blastp, d_busco, d_pfam, d_blastn, output_dir):
    '''Calculate and wriet score'''
    d_score = {}
    for gene_model in d_cds.keys():
        blast_score = d_blastp[gene_model]
        busco_score = d_busco[gene_model]
        pfam_score = d_pfam[gene_model]
        blastn_score = d_blastn[gene_model]
        d_score[gene_model] = sum(
            [blast_score, pfam_score, busco_score, blastn_score]
        )
    # Write score table
    outfile_score = os.path.join(output_dir, 'gene_model_scores.txt')
    outhandle_score = open(outfile_score, 'w')
    header_txt = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
        'software', 'software_id', 'blast_score', 'busco_score',
        'pfam_score', 'blastn_score', 'score_sum'
    )
    outhandle_score.write(header_txt)
    for tup in d_cds:
        software = tup[0]
        software_id = tup[1]
        blast_score = d_blastp[tup]
        busco_score = d_busco[tup]
        pfam_score = d_pfam[tup]
        blastn_score = d_blastn[tup]
        score_sum = sum([blast_score, busco_score, pfam_score, blastn_score])
        row_txt = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            software, software_id, round(blast_score, 1),
            round(busco_score, 1), round(pfam_score, 1), blastn_score,
            round(score_sum, 1)
        )
        outhandle_score.write(row_txt)
    outhandle_score.close()
    return d_score


def filtering(d_cds, d_cds_len, d_score):
    '''Filter good gene models'''
    # Sort CDS: three keys used - scaffold, start, and end
    d_cds_sorted = sorted(
        d_cds.items(), key=lambda x: (x[1][0], x[1][1], x[1][2])
    )
    # Find chunks
    model_chunks = []  # It will be list of list
    tmp_list = [d_cds_sorted[0][0]]  # Initialize
    i = 1
    max_cds_end = d_cds_sorted[0][1][2]
    while i < len(d_cds_sorted):
        previous_tup = d_cds_sorted[i - 1][1]
        current_gene_name, current_tup = d_cds_sorted[i]

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
    for model_chunk in model_chunks:
        final_gene_set += get_best_comb(model_chunk, d_cds, d_cds_len, d_score)
    return final_gene_set


def get_best_comb(model_chunk, d_cds, d_cds_len, d_score):
    '''Get the best combination of the genes'''
    # Sort by start point
    model_chunk2 = [(x, d_cds[x][1], d_cds[x][2]) for x in model_chunk]
    model_chunk_s = sorted(model_chunk2, key=lambda x: x[1])
    # (Score, cds_len, list of genes)
    init_score = d_score[model_chunk_s[0][0]]
    init_cds_len = d_cds_len[model_chunk_s[0][0]]
    combs = [(init_score, init_cds_len, [model_chunk_s[0]])]
    for current_model in model_chunk_s[1:]:
        flag = False
        for previous_model in combs[::-1]:
            if not is_overlap(current_model, previous_model[2][-1]):
                new_score = d_score[current_model[0]] + previous_model[0]
                new_cds_len = d_cds_len[current_model[1]] + previous_model[1]
                new_gene_list = previous_model[2] + [current_model]
                bisect.insort(combs, (new_score, new_cds_len, new_gene_list))
                flag = True
                break
        if not flag:
            # If there's no previous model to combine
            score = d_score[current_model[0]]
            cds_len = d_cds_len[current_model[0]]
            bisect.insort(combs, (score, cds_len, [current_model]))
    best_comb = [x[0] for x in combs[-1][2]]
    return best_comb
    

def is_overlap(model1, model2):
    '''Check if two models overlap'''
    start1, end1 = model1[1], model1[2]
    start2, end2 = model2[1], model2[2]
    overlap = min(end1, end2) - max(start1, start2) + 1
    condition1 = overlap >= (end1 - start1 + 1) * 0.1
    condition2 = overlap >= (end2 - start2 + 1) * 0.1
    return condition1 or condition2


def import_prot(nr_prot, d_mapping_rev):
    '''Import protein'''
    d_prot = defaultdict(str)
    prot_txt = import_file(nr_prot)
    for line in prot_txt:
        if re.search('^>', line):
            prot_id = line.split(' ')[0].replace('>', '')
            lst = d_mapping_rev[prot_id]

        else:
            for element in lst:
                prefix, prefix_id = element
                d_prot[(prefix, prefix_id)] += line
    return d_prot


def write_final_prots(final_gene_set, d_mapping, output_dir):
    '''Write final protein sequences'''
    prot_names = []
    for final_gene in final_gene_set:
        prot_name = d_mapping[final_gene]
        prot_names.append(prot_name)

    # Write to file
    outfile_protnames = os.path.join(output_dir, 'filtered_protnames.list')
    outhandle_protnames = open(outfile_protnames, 'w')
    outhandle_protnames.write('\n'.join(list(set(prot_names))))
    outhandle_protnames.close()


def write_files(
        genome_assembly, final_gene_set, d_gene, d_gff3, d_prot, d_exon,
        output_dir, d_cds):
    '''Write outputs'''
    d_scaffold = {}
    scaffold_i = 0
    genome_assembly_txt = import_file(genome_assembly)
    for line in genome_assembly_txt:
        if not line.startswith('>'):
            continue
        scaffold_name = line.split(' ')[0].replace('>', '')
        d_scaffold[scaffold_name] = scaffold_i
        scaffold_i += 1
    final_gene_set_sorted = sorted(
        final_gene_set, key=lambda x: (d_scaffold[d_gene[x][0]], d_cds[x][1])
    )
    output_gff3 = open(os.path.join(output_dir, 'filtered_1.gff3'), 'w')
    output_gff3.write('##gff-version 3\n')  # gff3 Header
    i = 1
    for gene in final_gene_set_sorted:
        gscaffold, gsource = d_gene[gene][:2]
        gstart, gend, gscore, gstrand, gphase = d_gene[gene][3:]
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

        if gene in d_exon:
            j = 1
            for exon in d_exon[gene]:
                escaffold, esource = exon[:2]
                estart, eend, escore, estrand, ephase = exon[3:]
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
            for cds in d_gff3[gene]:
                cscaffold, csource = cds[:2]
                cstart, cend, cscore, cstrand, cphase = cds[3:]
                cid = 'ID={}_{}.t1.e{};Parent={}_{}.t1'.format(
                    'gene', str(i).zfill(5), j, 'gene', str(i).zfill(5)
                )
                output_gff3.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    cscaffold, csource, 'exon', cstart, cend, cscore,
                    cstrand, '.', cid
                ))
                j += 1

        j = 1
        for cds in d_gff3[gene]:
            cscaffold, csource = cds[:2]
            cstart, cend, cscore, cstrand, cphase = cds[3:]
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
        while i < len(d_prot[gene]):
            output_prot.write('{}\n'.format(d_prot[gene][i:i + 60]))
            i += 60

    output_prot.close()


if __name__ == '__main__':
    main()
