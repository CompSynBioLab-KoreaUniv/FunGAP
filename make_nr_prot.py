#!/usr/bin/env python2

'''
Make nonredundant protein FASTA file.
'''

# Import modules
import os
import re
import sys
import mmap
from argparse import ArgumentParser
from collections import defaultdict


# Main function
def main(argv):
    optparse_usage = 'make_nr_prot.py -i <faa_files> -o <output_dir>'
    parser = ArgumentParser(usage=optparse_usage)
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
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)


def make_nr_prot(faa_files, output_dir):
    # Import FASTA & store in dictionary
    # (key: prot_seq, value: (prefix, name))
    D_nr_prot = defaultdict(list)
    for faa_file in faa_files:
        prefix = os.path.basename(faa_file).split('.')[0]
        D_faa = defaultdict(str)
        with open(faa_file, 'r+b') as f:
            map = mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)
            for line in iter(map.readline, ''):
                line = line.rstrip()
                if re.search('^>', line):
                    prot_name = line.split(' ')[0].replace('>', '')
                    continue

                D_faa[prot_name] += line.strip()

            for prot_name, seq in D_faa.items():
                D_nr_prot[seq].append((prefix, prot_name))

    # Write to FASTA & mapping file
    outfile1 = os.path.join(output_dir, 'nr_prot.faa')
    outfile2 = os.path.join(output_dir, 'nr_prot_mapping.txt')

    outhandle1 = open(outfile1, 'w')
    outhandle2 = open(outfile2, 'w')
    header_txt = '{}\t{}\t{}\n'.format('prot_name', 'software', 'software_id')
    outhandle2.write(header_txt)

    prot_num = 1
    for seq, lst in D_nr_prot.items():
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
            row_txt = '{}\t{}\t{}\n'.format(new_prot_name, software, software_id)
            outhandle2.write(row_txt)

    # Close handle
    outhandle1.close()
    outhandle2.close()


if __name__ == '__main__':
    main(sys.argv[1:])
