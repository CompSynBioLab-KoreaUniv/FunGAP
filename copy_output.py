#!/usr/bin/python

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
        "-o", "--output_dir", dest="output_dir", nargs=1,
        help="FunGAP output directory"
    )

    args = parser.parse_args()
    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir[0])
    else:
        print '[ERROR] Please provide OUTPUT DIRECTORY'
        sys.exit(2)

    # Run functions :) Slow is as good as Fast
    create_dir(output_dir)
    copy_output(output_dir)


def create_dir(output_dir):
    fungap_outdir = os.path.join(output_dir, 'fungap_out')
    if not os.path.exists(fungap_outdir):
        os.mkdir(fungap_outdir)


def copy_output(output_dir):
    gff3_out = os.path.join(output_dir, 'gpre_filtered/gpre_filtered.gff3')
    if not os.path.exists(gff3_out):
        print '\n[ERROR] %s does not exist. Please check previous steps' % (
            gff3_out
        )
        sys.exit(2)
    else:
        fungap_out_gff3 = os.path.join(
            output_dir, 'fungap_out/fungap_out.gff3'
        )
        copyfile(gff3_out, fungap_out_gff3)

    prot_out = os.path.join(output_dir, 'gpre_filtered/gpre_filtered_prot.faa')
    if not os.path.exists(prot_out):
        print '\n[ERROR] %s does not exist. Please check previous steps' % (
            prot_out
        )
        sys.exit(2)
    else:
        fungap_out_prot = os.path.join(
            output_dir, 'fungap_out/fungap_out_prot.faa'
        )
        copyfile(prot_out, fungap_out_prot)


if __name__ == "__main__":
    main(sys.argv[1:])
