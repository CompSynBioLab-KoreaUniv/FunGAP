#!/usr/bin/env python3

'''
Make nonredundant protein FASTA file.
Last updated: Aug 12, 2020
'''

import os
import re
from argparse import ArgumentParser
from collections import defaultdict


def main():
    '''Main function'''
    argparse_usage = 'make_nr_prot.py -i <faa_files> -o <output_dir>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-i', '--faa_files', nargs='+', required=True,
        help='Input FAA files'
    )
    parser.add_argument(
        '-o', '--output_dir', nargs='?', default='gene_filtering',
        help='Output directory'
    )

    args = parser.parse_args()
    faa_files = [os.path.abspath(x) for x in args.faa_files]
    output_dir = os.path.abspath(args.output_dir)

    # Run functions :)
    create_dir(output_dir)
    make_nr_prot(faa_files, output_dir)


def create_dir(output_dir):
    '''Create directory'''
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


def import_file(input_file):
    '''Import file'''
    with open(input_file) as f_in:
        txt = list(line.rstrip() for line in f_in)
    return txt


def make_nr_prot(faa_files, output_dir):
    '''Make non-redundant protein sequences'''
    # Import FASTA & store in dictionary
    # (key: prot_seq, value: (prefix, name))
    d_nr_prot = defaultdict(list)
    for faa_file in faa_files:
        prefix = re.sub(r'\.faa$', '', os.path.basename(faa_file))
        d_faa = defaultdict(str)
        faa_txt = import_file(faa_file)
        for line in faa_txt:
            line = line.rstrip()
            if re.search('^>', line):
                prot_name = line.split(' ')[0].replace('>', '')
                continue
            d_faa[prot_name] += line.strip()
        for prot_name, seq in d_faa.items():
            d_nr_prot[seq].append((prefix, prot_name))

    # Write to FASTA & mapping file
    outfile1 = os.path.join(output_dir, 'nr_prot.faa')
    outfile2 = os.path.join(output_dir, 'nr_prot_mapping.txt')

    outhandle1 = open(outfile1, 'w')
    outhandle2 = open(outfile2, 'w')
    header_txt = '{}\t{}\t{}\n'.format('prot_name', 'software', 'software_id')
    outhandle2.write(header_txt)

    prot_num = 1
    for seq, lst in d_nr_prot.items():
        new_prot_name = 'prot_{}'.format(prot_num)
        prot_num += 1

        # Write FASTA
        outhandle1.write('>{}\n'.format(new_prot_name))
        i = 0
        while i < len(seq):
            outhandle1.write('{}\n'.format(seq[i:i + 60]))
            i += 60

        # Write mapping file
        for element in lst:
            software, software_id = element
            row_txt = '{}\t{}\t{}\n'.format(
                new_prot_name, software, software_id
            )
            outhandle2.write(row_txt)

    # Close handle
    outhandle1.close()
    outhandle2.close()


if __name__ == '__main__':
    main()
