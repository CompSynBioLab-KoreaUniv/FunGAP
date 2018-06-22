#!/usr/bin/python

'''
Detect TE genes using Pfam
 - Input: protein FASTA
 - Output: TE gene list
Authour Byougnam Min on May 30, 2018
'''

# Import modules
import os
import sys
from datetime import datetime
from argparse import ArgumentParser
from distutils.spawn import find_executable

# Define D_te_pfam
D_te_pfam = {
    'PF00075': 'RNase H',
    'PF00078': 'Reverse transcriptase (RNA-dependent DNA polymerase)',
    'PF00665': 'Integrase core domain',
    'PF02925': 'Bacteriophage scaffolding protein D',
    'PF02992': 'Transposase family tnp2',
    'PF03184': 'DDE superfamily endonuclease',
    'PF03221': 'Tc5 transposase DNA-binding domain',
    'PF03732': 'Retrotransposon gag protein',
    'PF04687': 'Microvirus H protein (pilot protein)',
    'PF05699': 'hAT family C-terminal dimerisation region',
    'PF05840': 'Bacteriophage replication gene A protein (GPA)',
    'PF05970': 'PIF1-like helicase',
    'PF07727': 'Reverse transcriptase (RNA-dependent DNA polymerase)',
    'PF08283': 'Geminivirus rep protein central domain',
    'PF08284': 'Retroviral aspartyl protease',
    'PF10551': 'MULE transposase domain',
    'PF13358': 'DDE superfamily endonuclease',
    'PF13359': 'DDE superfamily endonuclease',
    'PF13456': 'Reverse transcriptase-like',
    'PF13837': 'Myb/SANT-like DNA-binding domain',
    'PF13976': 'GAG-pre-integrase domain',
    'PF14214': 'Helitron helicase-like domain at N-terminus',
    'PF14223': 'gag-polypeptide of LTR copia-type',
    'PF14529': 'Endonuclease-reverse transcriptase'
}


def main(argv):
    argparse_usage = (
        'detect_te_genes.py -p <protein_fasta> -i <interproscan_path>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-p", "--protein_fasta", dest="protein_fasta", nargs=1,
        help="Protein FASTA"
    )
    parser.add_argument(
        "-i", "--interproscan_path", dest="interproscan_path", nargs='?',
        help=(
            "InterProScan path (you do not need to provide this if it is in "
            "your PATH"
        )
    )

    args = parser.parse_args()
    if args.protein_fasta:
        protein_fasta = os.path.abspath(args.protein_fasta[0])
    else:
        print '[ERROR] Please provide PROTEIN FASTA'
        parser.print_help()
        sys.exit(2)

    if args.interproscan_path:
        interproscan_path = args.interproscan_path[0]
        interproscan_exe = os.path.join(interproscan_path, 'interproscan.sh')
        if os.path.exists(interproscan_exe):
            print '[ERROR] {} does not exist. Please check.'.format(
                interproscan_exe
            )
    else:
        interproscan_exe = find_executable("interproscan.sh")
        if not interproscan_exe:
            print '[ERROR] interproscan.sh not in your PATH. Please check'
            sys.exit(2)

    # Run functions :) Slow is as good as Fast
    ipr_out = run_interproscan(protein_fasta, interproscan_exe)
    detect_te_genes(ipr_out, protein_fasta)


def import_file(input_file):
    with open(input_file) as f_in:
        txt = list(line.rstrip() for line in f_in)
    return txt


def run_interproscan(protein_fasta, interproscan_exe):
    tmp_dir = 'tmp_interproscan'
    outfile_base = '{}_pfam'.format(os.path.splitext(protein_fasta)[0])
    outfile = '{}.tsv'.format(outfile_base)
    if not os.path.exists(outfile):
        command = (
            '{} -i {} -f tsv -appl PfamA --tempdir {} --output-file-base '
            '{}'.format(
                interproscan_exe, protein_fasta, tmp_dir, outfile_base
            )
        )
        current_time = datetime.now()
        current_time_f = current_time.strftime('%Y-%m-%d %H:%M')
        print '[{}] Start running InterproScan for Pfam'.format(current_time_f)
        print '[Run] {}'.format(command)
        os.system(command)
        current_time = datetime.now()
        current_time_f = current_time.strftime('%Y-%m-%d %H:%M')
        print '[{}] Done running InterproScan for Pfam'.format(current_time_f)
    else:
        print 'Running InterProscan has already been finished. Skip this step.'

    return '{}.tsv'.format(outfile_base)


def detect_te_genes(ipr_out, protein_fasta):
    outfile = '{}_te_pfam.txt'.format(os.path.splitext(protein_fasta)[0])
    outhandle = open(outfile, 'w')
    header_txt = '{}\t{}\t{}\n'.format('prot_id', 'pfam_id', 'pfam_desc')
    outhandle.write(header_txt)
    ipr_txt = import_file(ipr_out)
    num_te_genes = 0
    for line in ipr_txt:
        line_split = line.split('\t')
        prot_id = line_split[0]
        pfam_id = line_split[4]
        pfam_desc = line_split[5]
        if pfam_id not in D_te_pfam:
            continue
        num_te_genes += 1
        row_txt = '{}\t{}\t{}\n'.format(prot_id, pfam_id, pfam_desc)
        outhandle.write(row_txt)

    print '{} TE-related genes were found. Check {}'.format(
        num_te_genes, os.path.basename(outfile)
    )


if __name__ == "__main__":
    main(sys.argv[1:])
