# Import modules
import sys
import os
from Bio import SeqIO

'''
Check if inputs are proper to run FunGAP
'''

def check_inputs(
    trans_read_1, trans_read_2, trans_read_single, trans_bam, genome_assembly,
    sister_proteome,
):
    print 'Check input files...'

    if trans_read_1 and not os.path.exists(trans_read_1):
        sys.exit(
            '[ERROR] No such file {}'.format(os.path.basename(trans_read_1))
        )
    if trans_read_2 and not os.path.exists(trans_read_2):
         sys.exit(
             '[ERROR] No such file {}'.format(os.path.basename(trans_read_2))
         )
    if trans_read_single and not os.path.exists(trans_read_single):
        sys.exit(
            '[ERROR] No such file {}'.format(
                os.path.basename(trans_read_single)
            )
        )

    trans_read_files = []
    if trans_read_1 and trans_read_2:
        # Check extension
        if (
            not trans_read_1.endswith('_1.fastq') and
            not trans_read_1.endswith('_1.fq')
        ):
            error_message = (
                '[ERROR] TRANS_READ_1 file name is incorrect. '
                'Should be <prefix>_1.fastq. You provided {}'.format(
                    trans_read_1
                )
            )
            sys.exit(error_message)

        elif (
            not trans_read_2.endswith('_2.fastq') and
            not trans_read_2.endswith('_2.fq')
        ):
            error_message2 = (
                    '[ERROR] TRANS_READ_2 file name is incorrect. '
                    'Should be <prefix>_2.fastq. You provided {}'.format(
                        trans_read_1)
            )
            sys.exit(error_message2)

        # Check prefix
        prefix_1 = os.path.basename(trans_read_1).replace('_1.fastq', '')
        prefix_2 = os.path.basename(trans_read_2).replace('_2.fastq', '')
        if prefix_1 != prefix_2:
            error_message3 = (
                '[ERROR] Two paired-end trans_read_files should have same'
                'prefix. <prefix>_1.fastq and <prefix>_2.fastq'
            )
            sys.exit(error_message3)

        trans_read_files = [trans_read_1, trans_read_2]

    elif trans_read_single:
        # Check extension
        if (
            not trans_read_single.endswith('_s.fastq') and
            not trans_read_single.endswith('_s.fq')
        ):
            error_message4 = (
                '[ERROR] TRANS_READ_SINGLE file name is incorrect. Should be '
                '<prefix>_s.fastq. You provided {}'.format(trans_read_single)
            )
            sys.exit(error_message4)
        trans_read_files = [trans_read_single]

    elif trans_bam:
        trans_read_files = [trans_bam]

    if not trans_read_files:
        error_message5 = (
            '[ERROR] You did not provide any transcriptome files: '
            '-1 and -2, -U, or -A should be provided'
        )
        sys.exit(error_message5)

    print 'TRANS_READ_FILES is ok...'

    # Check geonme assembly in FASTA
    with open(genome_assembly, 'r') as handle:
        fasta = SeqIO.parse(handle, 'fasta')
        if not any(fasta):
            sys.exit('[ERROR] FASTA file is invalid: {}'.format(
                genome_assembly
            ))

        for record in SeqIO.parse(handle, 'fasta'):
            if '|' in record.id:
                error_message6 = (
                    '[ERROR] FASTA defline contains "|" character, please '
                    'remove and re-run'
                )
                sys.exit(error_message6)

    print 'GENOME_ASSEMBLY is ok...'

    # Check sister_proteome in FASTA
    with open(sister_proteome, 'r') as handle:
        fasta = SeqIO.parse(handle, 'fasta')
        if not any(fasta):
            sys.exit('[ERROR] FASTA file is invalid: {}'.format(
                genome_assembly
            ))

    print 'SISTER_PROTEOME is ok...'
    print ''

    return trans_read_files
