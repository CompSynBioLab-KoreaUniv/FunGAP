#!/usr/bin/python

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
from glob import glob
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
total_busco = 1438


# Main function
def main(argv):
    argparse_usage = (
        'filter_gff3s.py -i <input_gff3s> -m <mapping_file> -b <blast_dict> '
        '-B <busco_dict> -p <ipr_dict> -N <blastn_dict> -g <bad_dict> '
        '-n <nr_prot_file> -s <short_id> -o <output_prefix> -r <root_dir>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-i", "--input_gff3s", dest="input_gff3s", nargs='+',
        help="Multiple gff3 files"
    )
    parser.add_argument(
        "-m", "--mapping_file", dest="mapping_file", nargs=1,
        help="Mapping txt file (make_nr_prot.py)"
    )
    parser.add_argument(
        "-b", "--blast_dict", dest="blast_dict", nargs=2,
        help="Parsed blast output in dictionary"
    )
    parser.add_argument(
        "-B", "--busco_dict", dest="busco_dict", nargs=2,
        help="Parsed BUSCO output in dictionary"
    )
    parser.add_argument(
        "-p", "--ipr_dict", dest="ipr_dict", nargs=2,
        help="Parsed IPRscan output in dictionary"
    )
    parser.add_argument(
        "-N", "--blastn_dict", dest="blastn_dict", nargs=1,
        help="Parsed BLASTn output in dictionary"
    )
    parser.add_argument(
        "-g", "--bad_dict", dest="bad_dict", nargs=1,
        help="Parsed IPRscan output in dictionary"
    )
    parser.add_argument(
        "-n", "--nr_prot_file", dest="nr_prot_file", nargs=1,
        help="nr_prot.faa file (make_nr_prot.py)"
    )
    parser.add_argument(
        "-s", "--short_id", dest="short_id", nargs=1,
        help="Short ID for gene numbers"
    )
    parser.add_argument(
        "-o", "--output_prefix", dest="output_prefix", nargs=1,
        help="Output prefix"
    )
    parser.add_argument(
        "-r", "--root_dir", dest="root_dir", nargs=1,
        help=(
            'Root directory where log directory will be '
            'generated (default: ".")'), default=[os.getcwd()]
    )

    args = parser.parse_args()
    if args.input_gff3s:
        input_gff3s = [os.path.abspath(x) for x in args.input_gff3s]
    else:
        print '[ERROR] Please provide INPUT GFF3'
        sys.exit(2)

    if args.mapping_file:
        mapping_file = os.path.abspath(args.mapping_file[0])
    else:
        print '[ERROR] Please provide MAPPING TXT FILE'
        sys.exit(2)

    if args.blast_dict:
        blast_dict_score = os.path.abspath(args.blast_dict[0])
        blast_dict_evalue = os.path.abspath(args.blast_dict[1])
    else:
        print '[ERROR] Please provide BLAST DICT'
        sys.exit(2)

    if args.busco_dict:
        busco_dict_score = os.path.abspath(args.busco_dict[0])
        busco_dict_list = os.path.abspath(args.busco_dict[1])
    else:
        print '[ERROR] Please provide BUSCO DICT'
        sys.exit(2)

    if args.ipr_dict:
        ipr_dict_score = os.path.abspath(args.ipr_dict[0])
        ipr_dict_count = os.path.abspath(args.ipr_dict[1])
    else:
        print '[ERROR] Please provide IPR DICT PICKLEs'
        sys.exit(2)

    if args.blastn_dict:
        blastn_dict = os.path.abspath(args.blastn_dict[0])
    else:
        print '[ERROR] Please provide BLASTn DICT PICKLE'
        sys.exit(2)

    if args.bad_dict:
        bad_dict = os.path.abspath(args.bad_dict[0])
        D_bad = cPickle.load(open(bad_dict, 'rb'))
    else:
        print '[WARNING] Please provide BAD DICT PICKLE'
        D_bad = defaultdict(bool)

    if args.nr_prot_file:
        nr_prot_file = os.path.abspath(args.nr_prot_file[0])
    else:
        print '[ERROR] Please provide "nr_prot.faa" FILE'
        sys.exit(2)

    if args.short_id:
        short_id = args.short_id[0]
    else:
        print '[ERROR] Please provide SHORT ID'
        sys.exit(2)

    if args.output_prefix:
        output_prefix = args.output_prefix[0]
    else:
        print '[ERROR] Please provide OUTPUT PREFIX'
        sys.exit(2)

    root_dir = os.path.abspath(args.root_dir[0])

    # Create necessary dirs
    create_dir(root_dir)

    # Set logging
    log_file = os.path.join(
        root_dir, 'logs', 'pipeline', 'filter_gff3s.log'
    )
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is as good as Fast
    logger_time.debug('START: Filtering GFF3')
    D_mapping, D_mapping_rev = import_mapping(mapping_file)

    # Import dictionaries
    D_blast_score = cPickle.load(open(blast_dict_score, 'rb'))
    D_blast_evalue = cPickle.load(open(blast_dict_evalue, 'rb'))
    D_busco_score = cPickle.load(open(busco_dict_score, 'rb'))
    D_busco_list = cPickle.load(open(busco_dict_list, 'rb'))
    D_pfam_score = cPickle.load(open(ipr_dict_score, 'rb'))
    D_pfam_count = cPickle.load(open(ipr_dict_count, 'rb'))
    D_blastn_score = cPickle.load(open(blastn_dict, 'rb'))

    # Self-filtering to get stats
    D_stats = {}
    for input_gff3 in input_gff3s:
        prefix = os.path.basename(input_gff3).split('.')[0]
        D_gff3, D_gene, D_cds, D_cds_len, D_exon = import_gff3([input_gff3])
        self_filtered, stats = filtering(
            D_gene, D_cds, D_cds_len, D_mapping, D_blast_score, D_blast_evalue,
            D_busco_score, D_busco_list, D_pfam_score, D_pfam_count,
            D_blastn_score, D_bad, output_prefix
        )
        outfile_self = '%s_%s_filtered.list' % (output_prefix, prefix)
        outhandle_self = open(outfile_self, 'w')

        cds_len_filtered = 0
        for tup in self_filtered:
            outhandle_self.write('%s\n' % (tup[1]))
            cds_len_filtered += D_cds_len[tup]
        outhandle_self.close()

        (
            raw_num_genes, final_num_genes, blast_hit, pfam_hit,
            pfam_domains, busco_hit, blastn_hit
        ) = stats
        new_stats = (
            raw_num_genes, final_num_genes, blast_hit, pfam_hit,
            pfam_domains, busco_hit, blastn_hit, cds_len_filtered
        )
        D_stats[prefix] = new_stats

    # Filtering
    D_gff3, D_gene, D_cds, D_cds_len, D_exon = import_gff3(input_gff3s)
    final_gene_set, final_stats = filtering(
        D_gene, D_cds, D_cds_len, D_mapping, D_blast_score, D_blast_evalue,
        D_busco_score, D_busco_list, D_pfam_score, D_pfam_count,
        D_blastn_score, D_bad, output_prefix
    )
    D_prot = import_prot(nr_prot_file, D_mapping_rev)
    write_final_prots(final_gene_set, D_mapping, output_prefix)
    write_files(
        final_gene_set, D_gene, D_gff3, D_prot, D_exon, output_prefix, short_id
    )

    cds_len_final = 0
    for tup in final_gene_set:
        cds_len_final += D_cds_len[tup]

    (
        raw_num_genes, final_num_genes, blast_hit, pfam_hit,
        pfam_domains, busco_hit, blastn_hit
    ) = final_stats

    new_final_stats = (
        raw_num_genes, final_num_genes, blast_hit, pfam_hit,
        pfam_domains, busco_hit, blastn_hit, cds_len_final
    )
    D_stats['final'] = new_final_stats

    write_stats(D_stats, output_prefix)
    logger_time.debug('DONE : Filtering GFF3')


def import_file(path):
    with open(path) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)

    return txt


def create_dir(root_dir):
    log_dir = os.path.join(root_dir, 'logs')
    if not glob(log_dir):
        os.mkdir(log_dir)

    log_pipeline_dir = os.path.join(root_dir, 'logs', 'pipeline')
    if not glob(log_pipeline_dir):
        os.mkdir(log_pipeline_dir)


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
    reg_augustus_gene = re.compile('transcript_id "(\S+)";')
    for gff3_file in gff3_files:
        prefix = os.path.basename(gff3_file).split('.')[0]
        gff3_txt = import_file(gff3_file)

        D_parent = {}
        for line in gff3_txt:
            if not re.search('\t', line):
                continue  # Only consider line with a tab

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

            # BRAKER produces different format of gff
            m_augustus_gene = reg_augustus_gene.search(line)
            if m_augustus_gene:
                augustus_gene = m_augustus_gene.group(1)
            else:
                augustus_gene = ''

            # if feat_type == 'gene':
            #    gene_id = entry_id[:]

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
                if augustus_gene != '':
                    cds_gene = augustus_gene
                else:
                    cds_parent = parent_id
                    cds_gene = cds_parent  # mRNA ID

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
                if augustus_gene != '':
                    exon_gene = augustus_gene
                else:
                    exon_parent = parent_id
                    exon_gene = exon_parent  # mRNA ID

                new_exon_gene = (prefix, exon_gene)

                D_exon[new_exon_gene].append((
                    scaffold, source, feat_type, int(start), int(end),
                    score, strand, phase
                ))

    return D_gff3, D_gene, D_cds, D_cds_len, D_exon


def filtering(
    D_gene, D_cds, D_cds_len, D_mapping, D_blast_score, D_blast_evalue,
    D_busco_score, D_busco_list, D_pfam_score, D_pfam_count, D_blastn_score,
    D_bad, output_prefix
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
    outfile_score = '%s_score.txt' % (output_prefix)
    outhandle_score = open(outfile_score, 'w')
    header_txt = '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        'software', 'software_id', 'blast_score', 'busco_score',
        'pfam_score', 'blastn_score', 'score_sum'
    )
    outhandle_score.write(header_txt)
    for tup in D_cds_sorted:
        gene_tup = tup[0]
        software = gene_tup[0]
        software_id = gene_tup[1]
        blast_score = D_blast_score[gene_tup]
        busco_score = D_busco_score[gene_tup]
        pfam_score = D_pfam_score[gene_tup]
        blastn_score = D_blastn_score[gene_tup]
        score_sum = sum([blast_score, busco_score, pfam_score, blastn_score])
        row_txt = '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
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
                blast_score = D_blast_score[element]
                busco_score = D_busco_score[element]
                pfam_score = D_pfam_score[element]
                blastn_score = D_blastn_score[element]

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

    # Get some stats
    raw_num_genes = len(D_cds_sorted)
    final_num_genes = len(final_gene_set)
    blast_hit = 0
    pfam_hit = 0
    pfam_domains = 0
    busco_list = []
    blastn_hit = 0
    for final_gene in final_gene_set:
        check = final_gene in D_blast_evalue
        if check and D_blast_evalue[final_gene] < blast_cutoff:
            blast_hit += 1
        if D_pfam_score[final_gene] > 0:
            pfam_hit += 1
            pfam_domains += D_pfam_count[final_gene]
        if D_busco_list[final_gene]:
            busco_list += D_busco_list[final_gene]
        if D_blastn_score[final_gene]:
            blastn_hit += 1

    busco_hit = len(set(busco_list))
    stats = (
        raw_num_genes, final_num_genes, blast_hit, pfam_hit,
        pfam_domains, busco_hit, blastn_hit
    )

    return final_gene_set, stats


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


def write_final_prots(final_gene_set, D_mapping, output_prefix):
    prot_names = []
    for final_gene in final_gene_set:
        prot_name = D_mapping[final_gene]
        prot_names.append(prot_name)

    # Write to file
    outfile_protnames = '%s_final_protnames.txt' % (output_prefix)
    outhandle_protnames = open(outfile_protnames, 'w')
    outhandle_protnames.write('%s\n' % ('\n'.join(list(set(prot_names)))))
    outhandle_protnames.close()


def write_files(
    final_gene_set, D_gene, D_gff3, D_prot, D_exon, output_prefix, short_id
):

    final_gene_set_sorted = sorted(
        final_gene_set, key=lambda x: (D_gene[x][0], D_gene[x][3])
    )
    output_gff3 = open('%s.gff3' % (output_prefix), 'w')
    output_gff3.write('##gff-version 3\n')  # gff3 Header
    i = 1
    for gene in final_gene_set_sorted:
        (
            gscaffold, gsource, gfeat_type, gstart, gend, gscore,
            gstrand, gphase
        ) = D_gene[gene]
        gid = 'ID=%s_%s;source=%s:%s' % (
            short_id, str(i).zfill(5), gene[0], gene[1]
        )
        output_gff3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
            gscaffold, gsource, 'gene', gstart, gend, gscore,
            gstrand, gphase, gid
        ))

        mid = 'ID=%s_%s.t1;Parent=%s_%s' % (
            short_id, str(i).zfill(5), short_id, str(i).zfill(5)
        )
        output_gff3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
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
                eid = 'ID=%s_%s.t1.e%d;Parent=%s_%s.t1' % (
                    short_id, str(i).zfill(5), j, short_id, str(i).zfill(5)
                )
                row_txt = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
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
                cid = 'ID=%s_%s.t1.e%d;Parent=%s_%s.t1' % (
                    short_id, str(i).zfill(5), j, short_id, str(i).zfill(5)
                )
                output_gff3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
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
            cid = 'ID=%s_%s.t1.c%d;Parent=%s_%s.t1' % (
                short_id, str(i).zfill(5), j, short_id, str(i).zfill(5)
            )
            output_gff3.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                cscaffold, csource, 'CDS', cstart, cend, cscore,
                cstrand, cphase, cid
            ))
            j += 1

        i += 1

    output_gff3.close()

    # Write protein faa file
    output_prot = open('%s_prot.faa' % (output_prefix), 'w')
    for gene_num, gene in enumerate(final_gene_set_sorted, start=1):
        output_prot.write('>%s_%s.t1 source=%s:%s\n' % (
            short_id, str(gene_num).zfill(5), gene[0], gene[1]
        ))
        i = 0
        while i < len(D_prot[gene]):
            output_prot.write('%s\n' % (D_prot[gene][i:i + 60]))
            i += 60

    output_prot.close()


def write_stats(D_stats, output_prefix):
    outfile = '%s_stats_evidence.txt' % (output_prefix)
    outhandle = open(outfile, 'w')
    header_txt = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        'source', 'raw_models', 'filtered_models', 'blast_hit',
        'pfam_hit', 'pfam_domains', 'busco', 'blastn_hit', 'coding_region'
    )
    outhandle.write(header_txt)
    for prefix, stats in D_stats.items():
        if prefix == 'final':
            continue

        (
            raw_num_genes, final_num_genes, blast_hit, pfam_hit,
            pfam_domains, busco_hit, blastn_hit, coding_region
        ) = stats
        row_txt = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
            prefix, raw_num_genes, final_num_genes, blast_hit,
            pfam_hit, pfam_domains, busco_hit, blastn_hit, coding_region
        )

        outhandle.write(row_txt)

    final_stats = D_stats['final']
    (
        raw_num_genes, final_num_genes, blast_hit, pfam_hit,
        pfam_domains, busco_hit, blastn_hit, coding_region
    ) = final_stats
    row_txt = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
        'final', raw_num_genes, final_num_genes, blast_hit,
        pfam_hit, pfam_domains, busco_hit, blastn_hit, coding_region
    )
    outhandle.write(row_txt)
    outhandle.close()


if __name__ == "__main__":
    main(sys.argv[1:])
