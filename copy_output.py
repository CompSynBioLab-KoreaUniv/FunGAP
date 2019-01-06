#!/usr/bin/env python2

'''
Copy output to fungap_out directory
'''

# Import modules
import os
import sys
from shutil import copyfile
from argparse import ArgumentParser


# Main function
def main(argv):
    argparse_usage = 'copy_output.py -o <output_dir>'
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-o', '--output_dir', nargs=1, required=True,
        help='FunGAP output directory'
    )

    args = parser.parse_args()
    output_dir = os.path.abspath(args.output_dir[0])

    # Run functions :) Slow is as good as Fast
    create_dir(output_dir)
    copy_output(output_dir)


def create_dir(output_dir):
    fungap_outdir = os.path.join(output_dir, 'fungap_out')
    if not os.path.exists(fungap_outdir):
        os.mkdir(fungap_outdir)


def copy_output(output_dir):
    gff3_out = os.path.join(output_dir, 'gene_filtering/filtered_2.gff3')
    if not os.path.exists(gff3_out):
        sys.exit(
            '\n[ERROR] {} does not exist. Please check previous steps'.format(
                gff3_out
            )
        )
    else:
        fungap_out_gff3 = os.path.join(output_dir, 'fungap_out/fungap_out.gff3')
        copyfile(gff3_out, fungap_out_gff3)

    prot_out = os.path.join(output_dir, 'gene_filtering/filtered_prot.faa')
    if not os.path.exists(prot_out):
        sys.exit(
            '\n[ERROR] {} does not exist. Please check previous steps'.format(
                prot_out
            )
        )
    else:
        fungap_out_prot = os.path.join(
            output_dir, 'fungap_out/fungap_out_prot.faa'
        )
        copyfile(prot_out, fungap_out_prot)


if __name__ == '__main__':
    main(sys.argv[1:])
