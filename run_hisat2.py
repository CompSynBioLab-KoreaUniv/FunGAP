#!/usr/bin/python

# Import modules
import sys
import os
import re
from argparse import ArgumentParser

# Get Logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)
sys.path.append(this_dir)
from set_logging import set_logging


# Main function
def main(argv):
    argparser_usage = (
        'run_hisat2.py -r <fastq1> <fastq2> <fastq3> ... '
        '-o <output_dir> -l <log_dir> -f <ref_fasta> -c <num_cores>'
    )
    parser = ArgumentParser(usage=argparser_usage)
    parser.add_argument(
        "-r", "--read_files", dest="read_files", nargs='+',
        help='Multiople read files in fastq format'
    )
    parser.add_argument(
        "-o", "--output_dir", dest="output_dir", nargs=1,
        help='Output directory'
    )
    parser.add_argument(
        "-l", "--log_dir", dest="log_dir", nargs=1,
        help='Log directory'
    )
    parser.add_argument(
        "-f", "--ref_fasta", dest="ref_fasta", nargs=1,
        help='Reference fasta'
    )
    parser.add_argument(
        "-c", "--num_cores", dest="num_cores", nargs=1,
        help='Number of cores'
    )

    args = parser.parse_args()

    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir[0])
    else:
        print '[ERROR] Please provide proper OUTPUT DIRECTORY'
        parser.print_help()
        sys.exit(2)

    if args.log_dir:
        log_dir = os.path.abspath(args.log_dir[0])
    else:
        print '[ERROR] Please provide proper LOG DIRECTORY'
        parser.print_help()
        sys.exit(2)

    if args.read_files:
        read_files = [os.path.abspath(x) for x in args.read_files]
    else:
        print '[ERROR] Please provide proper READ FILES'
        parser.print_help()
        sys.exit(2)

    # Reference fasta
    if args.ref_fasta:
        ref_fasta = os.path.abspath(args.ref_fasta[0])
    else:
        print '[ERROR] Please provide proper file: REFERENCE FASTA'
        parser.print_help()
        sys.exit(2)

    if args.num_cores:
        num_cores = int(args.num_cores[0])
    else:
        num_cores = 1

    # Create necessary dirs
    create_dir(output_dir, log_dir)

    # Set logging
    log_file = os.path.join(log_dir, 'pipeline', 'run_hisat2.log')
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is as good as Fast
    logger_time.debug('START: Hisat2')
    hisat2_outputs = run_hisat2(
        read_files, output_dir, log_dir, ref_fasta, num_cores
    )
    post_process_sam(hisat2_outputs)
    logger_time.debug('DONE : Hisat2')


# Define functions
def create_dir(output_dir, log_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    log_output_dir = os.path.join(log_dir, os.path.basename(output_dir))
    if not os.path.exists(log_output_dir):
        os.mkdir(log_output_dir)

    log_pipeline_dir = os.path.join(log_dir, 'pipeline')
    if not os.path.exists(log_pipeline_dir):
        os.mkdir(log_pipeline_dir)


def run_hisat2(
    read_files, output_dir, log_dir, ref_fasta, num_cores
):
    output_dir = re.sub(r'/$', '', output_dir)
    output_base = os.path.basename(output_dir)

    # hisat2-build -p <num_cores> <ref_fasta> <ref_fasta>
    hisat2_build_log_file = os.path.join(
        log_dir, output_base, 'hisat2-build.log'
    )
    hisat2_build_output = '%s.5.ht2' % (ref_fasta)
    if not os.path.exists(hisat2_build_output):
        command1 = 'hisat2-build -p %s %s %s > %s 2>&1' % (
            num_cores, ref_fasta, ref_fasta, hisat2_build_log_file
        )
        logger_txt.debug('[Run] %s' % (command1))
        os.system(command1)
    else:
        logger_txt.debug('Running hisat2-build has already been finished')

    # hisat2 -p <num_cores> -x Choanephora_cucurbitarum_assembly.fna
    # -1 reads/chocu-mRNA_1.fastq -2 reads/chocu-mRNA_2.fastq
    # -S trans_hisat2/chocu-mRNA_with_annot.sam
    hisat2_outputs = []
    for read_file in read_files:
        if re.search(r'_R1.fastq', read_file):
            fastq_pair = read_file.replace('_1.fastq', '_2.fastq')
        elif re.search(r'_1.fastq', read_file):
            fastq_pair = read_file.replace('_1.fastq', '_2.fastq')
        else:
            continue

        prefix = os.path.basename(read_file).split('_')[0]
        hisat2_output = os.path.join(
            output_dir, '%s.sam' % (prefix)
        )
        hisat2_outputs.append(hisat2_output)
        sorted_bam_file = re.sub('.sam$', '_sorted.bam', hisat2_output)
        if not os.path.exists(sorted_bam_file):
            log_file = os.path.join(
                log_dir, output_base, '%s_%s.log' % (output_base, prefix)
            )
            command2 = 'hisat2 -p %s -x %s -1 %s -2 %s -S %s > %s 2>&1' % (
                num_cores, ref_fasta, read_file, fastq_pair,
                hisat2_output, log_file
            )
            logger_txt.debug('[Run] %s' % (command2))
            os.system(command2)
        else:
            logger_txt.debug(
                'Ruuning Hisat2 has already been finished for %s' % (prefix)
            )

    return hisat2_outputs


def post_process_sam(hisat2_outputs):
    # samtools view -Sb <SAMFILE> > <BAMFILE>
    bam_files = []
    for hisat2_output in hisat2_outputs:
        bam_file = re.sub('.sam$', '.bam', hisat2_output)
        bam_files.append(bam_file)
        if not os.path.exists(bam_file):
            command1 = 'samtools view -Sb %s > %s' % (hisat2_output, bam_file)
            logger_txt.debug('[Run] %s' % (command1))
            os.system(command1)
        else:
            logger_txt.debug((
                'Converting SAM TO BAM has already been finisehd for %s'
            ) % (os.path.basename(bam_file)))

    # samtools sort <test.bam> <test_sorted>
    # samtools index <test_sorted.bam> <test_sorted.bai>
    for bam_file in bam_files:
        sorted_bam_prefix = re.sub('.bam$', '_sorted', bam_file)
        sorted_bam_file = re.sub('.bam$', '_sorted.bam', bam_file)
        index_file = re.sub('.bam$', '.bai', sorted_bam_file)
        if not os.path.exists(index_file):
            command2 = 'samtools sort %s -o %s.bam' % (
                bam_file, sorted_bam_prefix
            )
            logger_txt.debug('[Run] %s' % (command2))
            os.system(command2)
            command3 = 'samtools index %s %s' % (
                sorted_bam_file, index_file
            )
            logger_txt.debug('[Run] %s' % (command3))
            os.system(command3)
        else:
            logger_txt.debug('Soring has already been finished for %s' % (
                os.path.basename(sorted_bam_prefix)
            ))

if __name__ == "__main__":
    main(sys.argv[1:])
