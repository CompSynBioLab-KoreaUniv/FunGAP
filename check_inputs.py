#!/usr/bin/env python3

'''
Check if inputs are proper to run FunGAP
Last updated: Jul 13, 2020
'''

import re
import os
import subprocess
import sys

from Bio import SeqIO

from import_config import import_config


def check_inputs(
        trans_read_1, trans_read_2, trans_read_single, trans_bam,
        genome_assembly, sister_proteome, busco_dataset):
    '''Check input formats'''
    print('Check input files...')
    check_files_exists(trans_read_1, trans_read_2, trans_read_single)
    trans_read_files = check_trans(
        trans_read_1, trans_read_2, trans_read_single, trans_bam
    )
    check_assmebly(genome_assembly)
    check_sister_proteome(sister_proteome)
    check_busco_dataset(busco_dataset)
    print('')
    return trans_read_files


def check_files_exists(trans_read_1, trans_read_2, trans_read_single):
    '''Check if input files exist'''
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


def check_trans(trans_read_1, trans_read_2, trans_read_single, trans_bam):
    '''Check transcriptome files'''
    trans_read_files = []
    if trans_read_1 and trans_read_2:
        # Check extension
        cond1 = trans_read_1.endswith('_1.fastq')
        cond2 = trans_read_1.endswith('_1.fq')
        cond3 = trans_read_2.endswith('_2.fastq')
        cond4 = trans_read_2.endswith('_2.fq')
        if not cond1 and not cond2:
            error_message = (
                '[ERROR] TRANS_READ_1 file name is incorrect. '
                'Should be <prefix>_1.fastq. You provided {}'.format(
                    trans_read_1
                )
            )
            sys.exit(error_message)
        elif not cond3 and not cond4:
            error_message2 = (
                '[ERROR] TRANS_READ_2 file name is incorrect. '
                'Should be <prefix>_2.fastq. You provided {}'.format(
                    trans_read_1
                )
            )
            sys.exit(error_message2)

        # Check prefix
        prefix_1 = (
            os.path.basename(trans_read_1)
            .replace('_1.fastq', '')
            .replace('_1.fq', '')
        )
        prefix_2 = (
            os.path.basename(trans_read_2)
            .replace('_2.fastq', '')
            .replace('_2.fq', '')
        )
        if prefix_1 != prefix_2:
            error_message3 = (
                '[ERROR] Two paired-end trans_read_files should have same'
                'prefix. <prefix>_1.f(ast)q and <prefix>_2.f(ast)q'
            )
            sys.exit(error_message3)

        trans_read_files = [trans_read_1, trans_read_2]

    elif trans_read_single:
        # Check extension
        cond5 = trans_read_single.endswith('_s.fastq')
        cond6 = trans_read_single.endswith('_s.fq')
        if not cond5 and not cond6:
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

    print('TRANS_READ_FILES is ok...')
    return trans_read_files


def check_assmebly(genome_assembly):
    '''Check geonme assembly in FASTA'''
    with open(genome_assembly, 'r') as handle:
        fasta = SeqIO.parse(handle, 'fasta')
        if not any(fasta):
            sys.exit('[ERROR] FASTA file is invalid: {}'.format(
                genome_assembly
            ))

        error_message6 = (
            '[ERROR] FASTA defline contains "|" character, please remove and '
            're-run'
        )
        for record in SeqIO.parse(handle, 'fasta'):
            if '|' in record.id:
                sys.exit(error_message6)

    print('GENOME_ASSEMBLY is ok...')


def check_sister_proteome(sister_proteome):
    '''Check sister_proteome in FASTA'''
    with open(sister_proteome, 'r') as handle:
        fasta = SeqIO.parse(handle, 'fasta')
        if not any(fasta):
            sys.exit('[ERROR] FASTA file is invalid: {}'.format(
                sister_proteome
            ))
    print('SISTER_PROTEOME is ok...')


def check_busco_dataset(busco_dataset):
    '''Check BUSCO dataset'''
    d_conf = import_config()
    busco_bin = d_conf['BUSCO_PATH']
    proc = subprocess.Popen(
        [busco_bin, '--list-datasets'], stdout=subprocess.PIPE
    )
    output = str(proc.stdout.read().decode('utf-8'))
    busco_dbs = re.findall(r'\S+_odb10', output)
    if busco_dataset not in set(busco_dbs):
        sys.exit(
            '[ERROR] Invalid BUSCO DATASET: {}. Run busco --list-datasets to '
            'get a full list available datasets'.format(busco_dataset)
        )
    print('BUSCO_DATASET is ok...')
