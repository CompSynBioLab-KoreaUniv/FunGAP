#!/usr/bin/env python3

'''
Wrapper gene prediction pipeline

This script runs the pipeline with following order.
    1) Preprocessing
        check_inputs.py
        run_hisat2.py
        run_trinity.py
        run_repeat_modeler.py

    2) Gene prediction
        run_augustus.py
        run_maker.py
        run_braker.py

    3) Evaluation and filtering
        make_nr_prot.py
        make_transcripts.py
        run_busco.py
        run_pfam_scan.py
        run_blastp.py
        run_blastn.py
        import_blastp.py
        import_busco.py
        import_pfam.py
        import_blastn.py
        catch_bad_genes.py
        filter_gff3s.py
        gff3_postprocess.py

    4) Write output
        copy_output.py
        create_markdown.py
Last updated: Jul 13, 2020
'''

# Version
__version__ = '1.1.0'

# Import modules
import os
import re
import shlex
from argparse import ArgumentParser
from datetime import datetime
from glob import glob
from subprocess import check_call

from check_inputs import check_inputs
from set_logging import set_logging


# Main function
def main():
    '''Main function'''
    argparse_usage = (
        'fungap.py -g <genome_assembly> -12UA <trans_read_files> '
        '-o <output_dir> -a <augustus_species> -b <busco_dataset> '
        '-s <sister_proteome>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-o', '--output_dir', nargs='?', default='fungap_out',
        help='Output directory (default: fungap_out)'
    )
    parser.add_argument(
        '-1', '--trans_read_1', nargs='?', default='',
        help='Paired-end read1 "<prefix>_1.fastq"'
    )
    parser.add_argument(
        '-2', '--trans_read_2', nargs='?', default='',
        help='Paired-end read2 "<prefix>_2.fastq"'
    )
    parser.add_argument(
        '-U', '--trans_read_single', nargs='?', default='',
        help='Single read "<prefix>_s.fastq"'
    )
    parser.add_argument(
        '-A', '--trans_bam', nargs='?', default='',
        help='BAM file (RNA-seq reads alignment to a genome assembly'
    )
    parser.add_argument(
        '-g', '--genome_assembly', nargs=1, required=True,
        help='Genome assembly file in FASTA format'
    )
    parser.add_argument(
        '-a', '--augustus_species', nargs=1, required=True,
        help='AUGUSTUS species'
    )
    parser.add_argument(
        '-b', '--busco_dataset', nargs=1, required=True,
        help=(
            'BUSCO lineage dataset (run "busco --list-datasets" for the list). '
            'Commonly used are fungi_odb10, ascomycota_odb10, and '
            'basidiomycota_odb10'
        )
    )
    parser.add_argument(
        '-s', '--sister_proteome', nargs=1, required=True,
        help='Sister proteome sequences in .faa'
    )
    parser.add_argument(
        '-c', '--num_cores', nargs='?', default=1, type=int,
        help='Number of cores to be used (default: 1)'
    )
    parser.add_argument(
        '-v', '--version', action='version',
        version='%(prog)s {}'.format(__version__)
    )

    # Options for non-fungus genome
    parser.add_argument(
        '--no_braker_fungus', action='store_true',
        help='No --fungus flag in BRAKER for non-fungus genomes'
    )
    parser.add_argument(
        '--no_jaccard_clip', action='store_true',
        help='No --jaccard_clip flag in Trinity for non-fungus genomes'
    )
    parser.add_argument(
        '--no_genemark_fungus', action='store_true',
        help='No --fungus flag in GeneMark for non-fungus genomes'
    )
    parser.add_argument(
        '-M', '--max_intron', nargs='?', default=2000, type=int,
        help='Max intron length (Default: 2000 bp)'
    )

    args = parser.parse_args()
    d_path = get_paths()
    output_dir = os.path.abspath(args.output_dir)
    trans_read_1 = args.trans_read_1
    trans_read_2 = args.trans_read_2
    trans_read_single = args.trans_read_single
    trans_bam = args.trans_bam
    genome_assembly = os.path.abspath(args.genome_assembly[0])
    augustus_species = args.augustus_species[0]
    busco_dataset = args.busco_dataset[0]
    sister_proteome = os.path.abspath(args.sister_proteome[0])
    num_cores = args.num_cores
    max_intron = args.max_intron

    # For non-fungus genomes
    if args.no_braker_fungus:
        no_braker_fungus = ''
    else:
        no_braker_fungus = '--fungus'

    if args.no_jaccard_clip:
        no_jaccard_clip = ''
    else:
        no_jaccard_clip = '--jaccard_clip'

    if args.no_genemark_fungus:
        no_genemark_fungus = ''
    else:
        no_genemark_fungus = '--gmes_fungus'

    # Create nessasary dirs
    create_dir(output_dir)

    # Set logging
    log_file = os.path.join(output_dir, 'logs', 'fungap.log')
    logger = set_logging(log_file)
    logger_time, logger_txt = logger

    logger_txt.debug('\n============ New Run %s ============', datetime.now())

    # Run functions :) Slow is as good as Fast
    trans_read_files = check_inputs(
        trans_read_1, trans_read_2, trans_read_single, trans_bam,
        genome_assembly, sister_proteome, busco_dataset
    )
    trans_bams = run_hisat2(
        genome_assembly, trans_read_files, output_dir, num_cores, max_intron,
        d_path, logger
    )
    trinity_asms = run_trinity(
        trans_bams, output_dir, num_cores, no_jaccard_clip, max_intron, d_path,
        logger
    )
    repeat_model_file = run_repeat_modeler(
        genome_assembly, output_dir, num_cores, d_path, logger
    )
    maker_gff3s, maker_faas = run_maker(
        genome_assembly, output_dir, augustus_species, sister_proteome,
        num_cores, repeat_model_file, trinity_asms, no_genemark_fungus, d_path,
        logger
    )
    # Get masked assembly
    masked_assembly = os.path.join(
        output_dir, 'maker_out', 'masked_assembly.fasta'
    )

    # Run Augustus
    augustus_gff3, augustus_faa = run_augustus(
        masked_assembly, output_dir, augustus_species, d_path, logger
    )

    # Run Braker
    braker_gff3s, braker_faas = run_braker(
        masked_assembly, trans_bams, output_dir, num_cores, no_braker_fungus,
        d_path, logger
    )

    # Run BUSCO on each gene models
    faa_files = [augustus_faa] + maker_faas + braker_faas
    for faa_file in faa_files:
        run_busco(faa_file, output_dir, busco_dataset, d_path, logger)
    busco_out_dir = os.path.join(output_dir, 'busco_out')

    # Get protein nr by removing identical proteins
    nr_prot_file, nr_prot_mapping_file = make_nr_prot(
        faa_files, output_dir, d_path, logger
    )

    # Run BLASTp with nr prot file
    blastp_output = run_blastp(
        nr_prot_file, output_dir, sister_proteome, num_cores, d_path, logger
    )

    # Run Pfam_scan with nr prot file
    pfam_scan_out = run_pfam_scan(
        nr_prot_file, output_dir, num_cores, d_path, logger
    )

    # Concatenate all transcripts files
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    trinity_asm = os.path.join(gene_filtering_dir, 'trinity_transcripts.fna')
    command = 'cat {} > {}'.format(' '.join(trinity_asms), trinity_asm)
    logger_time.debug('Create transcript')
    logger_txt.debug('[Run] %s', command)
    os.system(command)

    gff3_files = [augustus_gff3] + maker_gff3s + braker_gff3s
    blastn_out_files = []
    for gff3_file in gff3_files:
        transcript_file = make_transcripts(
            genome_assembly, gff3_file, d_path, logger
        )
        blastn_out_file = run_blastn(
            transcript_file, trinity_asm, output_dir, d_path, logger
        )
        blastn_out_files.append(blastn_out_file)

    # Import BLAST, BUSCO and Pfam score
    blastp_dict = import_blastp(
        blastp_output, nr_prot_mapping_file, d_path, logger
    )
    busco_dict = import_busco(busco_out_dir, output_dir, d_path, logger)
    pfam_dict = import_pfam(pfam_scan_out, nr_prot_mapping_file, d_path, logger)
    blastn_dict = import_blastn(blastn_out_files, output_dir, d_path, logger)

    # Catch bad genes
    bad_dict = catch_bad_genes(
        gff3_files, genome_assembly, output_dir, d_path, logger
    )
    filter_gff3s(
        genome_assembly, gff3_files, blastp_dict, busco_dict, pfam_dict,
        blastn_dict, bad_dict, nr_prot_file, nr_prot_mapping_file, output_dir,
        d_path, logger
    )
    gff3_postprocess(genome_assembly, output_dir, d_path, logger)

    # Copy output files
    copy_output(output_dir, d_path, logger)

    # Create markdown
    create_markdown(
        genome_assembly, output_dir, trans_bams, trinity_asms, d_path, logger
    )


def get_paths():
    '''Get program paths'''
    this_path = os.path.realpath(__file__)
    this_dir = os.path.dirname(this_path)

    d_path = {}
    d_path['check_dependencies'] = os.path.join(
        this_dir, 'check_dependencies.py'
    )
    d_path['run_hisat2'] = os.path.join(this_dir, 'run_hisat2.py')
    d_path['run_trinity'] = os.path.join(this_dir, 'run_trinity.py')
    d_path['run_repeat_modeler'] = os.path.join(
        this_dir, 'run_repeat_modeler.py'
    )

    d_path['run_augustus'] = os.path.join(this_dir, 'run_augustus.py')
    d_path['run_maker'] = os.path.join(this_dir, 'run_maker.py')
    d_path['run_braker'] = os.path.join(this_dir, 'run_braker.py')

    d_path['run_busco'] = os.path.join(this_dir, 'run_busco.py')
    d_path['run_pfam_scan'] = os.path.join(this_dir, 'run_pfam_scan.py')
    d_path['make_nr_prot'] = os.path.join(this_dir, 'make_nr_prot.py')
    d_path['run_blastp'] = os.path.join(this_dir, 'run_blastp.py')
    d_path['make_transcripts'] = os.path.join(this_dir, 'make_transcripts.py')
    d_path['run_blastn'] = os.path.join(this_dir, 'run_blastn.py')
    d_path['import_blast'] = os.path.join(this_dir, 'import_blastp.py')
    d_path['import_busco'] = os.path.join(this_dir, 'import_busco.py')
    d_path['import_pfam'] = os.path.join(this_dir, 'import_pfam.py')
    d_path['import_blastn'] = os.path.join(this_dir, 'import_blastn.py')
    d_path['catch_bad_genes'] = os.path.join(this_dir, 'catch_bad_genes.py')
    d_path['filter_gff3s'] = os.path.join(this_dir, 'filter_gff3s.py')
    d_path['gff3_postprocess'] = os.path.join(this_dir, 'gff3_postprocess.py')

    d_path['copy_output'] = os.path.join(this_dir, 'copy_output.py')
    d_path['create_markdown'] = os.path.join(this_dir, 'create_markdown.py')
    return d_path


def create_dir(output_dir):
    '''Create directory'''
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    log_dir = os.path.join(output_dir, 'logs')
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)


def run_hisat2(
        genome_assembly, trans_read_files, output_dir, num_cores, max_intron,
        d_path, logger):
    '''Run HISAT2'''
    if len(trans_read_files) == 1 and trans_read_files[0].endswith('.bam'):
        return trans_read_files

    hisat2_output_dir = os.path.join(output_dir, 'hisat2_out')
    log_dir = os.path.join(output_dir, 'logs')

    # run_hisat2.py -r <fastq1> <fastq2> <fastq3> ... \
    # -o <output_dir> -l <log_dir> -f <ref_fasta> -c <num_cores>
    # -m <max_intron>
    command = (
        '{} --read_files {} --output_dir {} --log_dir {} --ref_fasta {} '
        '--num_cores {} --max_intron {}'.format(
            d_path['run_hisat2'], ' '.join(trans_read_files), hisat2_output_dir,
            log_dir, genome_assembly, num_cores, max_intron
        )
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_run_hisat2')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_hisat2\n')

    # Get output BAM file paths
    trans_bams = []
    for trans_read_file in trans_read_files:
        prefix = re.sub(
            r'_[12s]$', '',
            os.path.basename(os.path.splitext(trans_read_file)[0])
        )
        hisat2_output = os.path.join(hisat2_output_dir, '{}.bam'.format(prefix))
        trans_bams.append(hisat2_output)
    trans_bams2 = list(set(trans_bams))
    return trans_bams2


def run_trinity(
        trans_bams, output_dir, num_cores, no_jaccard_clip, max_intron, d_path,
        logger):
    '''Run Trinity'''
    trinity_output_dir = os.path.join(output_dir, 'trinity_out')
    log_dir = os.path.join(output_dir, 'logs')
    # run_trinity.py -b <bam_files> -o <output_dir> -l <log_dir> -c <num_cores>
    # -m <max_intron> --jaccard_clip
    command = (
        '{} --bam_files {} --output_dir {} --log_dir {} --num_cores {} '
        '--max_intron {} {}'.format(
            d_path['run_trinity'], ' '.join(trans_bams), trinity_output_dir,
            log_dir, num_cores, max_intron,
            no_jaccard_clip
        )
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_run_trinity')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_trinity\n')

    # Get output transcriptome assembly files
    trinity_asms = glob(
        os.path.join(output_dir, 'trinity_out', '*/Trinity_*.fasta')
    )
    return trinity_asms


def run_repeat_modeler(genome_assembly, output_dir, num_cores, d_path, logger):
    '''Run RepeatModeler'''
    # run_repeat_modeler.py -g <genome_assembly> -o <output_dir> -l <log_dir>
    # -c <num_cores>
    rm_output_dir = os.path.join(output_dir, 'repeat_modeler_out')
    log_dir = os.path.join(output_dir, 'logs')
    command = (
        '{} --genome_assembly {} --output_dir {} --log_dir {} '
        '--num_cores {}'.format(
            d_path['run_repeat_modeler'], genome_assembly, rm_output_dir,
            log_dir, num_cores
        )
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_run_repeat_modeler')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_repeat_modeler\n')

    repeat_model_file = glob(
        os.path.join(rm_output_dir, 'RM*/consensi.fa.classified')
    )[0]
    return repeat_model_file


def run_maker(
        genome_assembly, output_dir, augustus_species, sister_proteome,
        num_cores, repeat_model_file, trinity_asms, no_genemark_fungus, d_path,
        logger):
    '''Run Maker'''
    maker_out_dir = os.path.join(output_dir, 'maker_out')
    # run_maker.py -i <input_fasta> -a <augustus_species> -p <protein_db_fasta>
    # -R <repeat_model> -e <est_files> -o <output_dir> -c <num_cores>
    # -l <log_dir> --gmes_fungus
    log_dir = os.path.join(output_dir, 'logs')
    command = (
        '{} --input_fasta {} --augustus_species {} --protein_db_fasta {}'
        ' --repeat_model {} --est_files {} --output_dir {} --num_cores {} '
        '--log_dir {} {}'.format(
            d_path['run_maker'], genome_assembly, augustus_species,
            sister_proteome, repeat_model_file, ' '.join(trinity_asms),
            maker_out_dir, num_cores, log_dir, no_genemark_fungus
        )
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_run_maker')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_maker\n')

    maker_gff3s = glob(
        os.path.join(output_dir, 'maker_out', '*/maker_*.gff3')
    )
    maker_faas = glob(os.path.join(output_dir, 'maker_out', '*/maker_*.faa'))
    return maker_gff3s, maker_faas


def run_augustus(masked_assembly, output_dir, augustus_species, d_path, logger):
    '''Run Augustus'''
    # run_augustus.py -m <masked_assembly> -s <species> -o <output_dir>
    # -l <log_dir>
    output_dir = os.path.join(output_dir, 'augustus_out')
    log_dir = os.path.join(output_dir, 'logs')
    command = (
        '{} --masked_assembly {} --species {} --output_dir {} '
        '--log_dir {}'.format(
            d_path['run_augustus'], masked_assembly, augustus_species,
            output_dir, log_dir
        )
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_run_augustus')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_augustus\n')
    augustus_gff3 = os.path.join(output_dir, 'augustus.gff3')
    augustus_faa = os.path.join(output_dir, 'augustus.faa')
    return augustus_gff3, augustus_faa


def run_braker(
        masked_assembly, trans_bams, output_dir, num_cores, no_braker_fungus,
        d_path, logger):
    '''Run BRAKER'''
    braker_output_dir = os.path.join(output_dir, 'braker_out')
    log_dir = os.path.join(output_dir, 'logs')
    # run_braker.py -m <masked_assembly> -b <bam_files> -o <output_dir>
    # -l <log_dir> -c <num_cores> --fungus
    command = (
        '{} --masked_assembly {} --bam_files {} --output_dir {} '
        '--log_dir {} --num_cores {} {}'.format(
            d_path['run_braker'], masked_assembly, ' '.join(trans_bams),
            braker_output_dir, log_dir, num_cores, no_braker_fungus
        )
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_run_braker')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_braker\n')

    prefixes = [os.path.basename(os.path.splitext(x)[0]) for x in trans_bams]
    prefixes_u = list(set(prefixes))

    braker_gff3s = []
    braker_faas = []
    for prefix in prefixes_u:
        braker_gff3 = os.path.join(
            output_dir, 'braker_out', prefix, 'braker_{}.gff3'.format(prefix)
        )
        braker_gff3s.append(braker_gff3)
        braker_faa = os.path.join(
            output_dir, 'braker_out', prefix, 'braker_{}.faa'.format(prefix)
        )
        braker_faas.append(braker_faa)
    return braker_gff3s, braker_faas


def run_busco(input_faa, output_dir, busco_dataset, d_path, logger):
    '''Run BUSCO'''
    busco_output_dir = os.path.join(output_dir, 'busco_out')
    log_dir = os.path.join(output_dir, 'logs')
    # run_busco.py -i <input_fasta> -d <lineage_dataset> -o <output_dir>
    # -l <log_dir>
    command = (
        '{} --input_fasta {} --lineage_dataset {} --output_dir {} --log_dir {}'
    ).format(
        d_path['run_busco'], input_faa, busco_dataset, busco_output_dir,
        log_dir
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_run_busco')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_busco\n')


def make_nr_prot(faa_files, output_dir, d_path, logger):
    '''Make non-redundant protein sequences'''
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    # make_nr_prot.py -i <faa_files> -o <output_dir>
    command = '{} --faa_files {} --output_dir {}'.format(
        d_path['make_nr_prot'], ' '.join(faa_files), gene_filtering_dir
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_make_nr_prot')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_make_nr_prot\n')

    nr_prot_file = os.path.join(gene_filtering_dir, 'nr_prot.faa')
    nr_prot_mapping_file = os.path.join(
        gene_filtering_dir, 'nr_prot_mapping.txt'
    )
    return nr_prot_file, nr_prot_mapping_file


def run_blastp(
        nr_prot_file, output_dir, sister_proteome, num_cores, d_path, logger):
    '''Run BLASTp against sister protein sequences'''
    # run_blastp.py -q <query_fasta> -d <db_fasta> -l <log_dir> -c <num_cores>
    log_dir = os.path.join(output_dir, 'logs')
    command = (
        '{} --query_fasta {} --db_fasta {} --log_dir {} '
        '--num_cores {}'.format(
            d_path['run_blastp'], nr_prot_file, sister_proteome, log_dir,
            num_cores
        )
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_run_blastp')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_blastp\n')

    blastp_output = os.path.join(output_dir, 'gene_filtering', 'nr_prot.blastp')
    return blastp_output


def run_pfam_scan(nr_prot_file, output_dir, num_cores, d_path, logger):
    '''Run pfam_scan'''
    # run_pfam_scan.py -i <input_fasta> -l <log_dir> -c <num_cores>
    log_dir = os.path.join(output_dir, 'logs')
    command = '{} --input_fasta {} --log_dir {} --num_cores {}'.format(
        d_path['run_pfam_scan'], nr_prot_file, log_dir, num_cores
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_run_pfam_scan')
    logger_txt.debug('[Wapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_pfam_scan\n')

    pfam_scan_out = os.path.join(
        output_dir, 'gene_filtering', 'nr_prot.pfam_scan'
    )
    return pfam_scan_out


def make_transcripts(genome_assembly, gff3_file, d_path, logger):
    '''Make transcript'''
    # make_transcripts.py -f <input_fasta> -g <input_gff3>
    command = '{} --input_fasta {} --input_gff3 {}'.format(
        d_path['make_transcripts'], genome_assembly, gff3_file
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_make_transcripts')
    logger_txt.debug('[Wapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_make_transcripts\n')

    gff3_base = os.path.splitext(gff3_file)[0]
    transcript_file = '{}_transcript.fna'.format(gff3_base)
    return transcript_file


def run_blastn(
        predicted_transcript, assembled_transcript, output_dir, d_path,
        logger):
    '''Run BLASTn against transcripts'''
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    prefix = re.sub(
        r'_transcript\.fna', '', os.path.basename(predicted_transcript)
    )
    out_prefix = os.path.join(gene_filtering_dir, prefix)
    log_dir = os.path.join(output_dir, 'logs')
    # run_blastn.py -q <query_fasta> -d <db_fasta> -o <output_prefix>
    # -l <log_dir> -c <num_cores>
    command = (
        '{} --query_fasta {} --db_fasta {} --output_prefix {} '
        '--log_dir {}'.format(
            d_path['run_blastn'], predicted_transcript, assembled_transcript,
            out_prefix, log_dir
        )
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_run_blastn')
    logger_txt.debug('[Wapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_blastn\n')

    blastn_out = '{}.blastn'.format(out_prefix)
    return blastn_out


def import_blastp(blastp_output, nr_prot_mapping_file, d_path, logger):
    '''Import BLASTp output'''
    # import_blastp.py -b <blastp_out_file> -n <nr_prot_mapping>
    blastp_out_dir = os.path.dirname(blastp_output)
    command = '{} --blastp_out_file {} --nr_prot_mapping {}'.format(
        d_path['import_blast'], blastp_output, nr_prot_mapping_file
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_import_blastp')
    logger_txt.debug('[Wrapper] %s', command)
    os.system(command)
    logger_time.debug('DONE : wrapper_import_blastp\n')

    # Return files
    blastp_dict = os.path.join(blastp_out_dir, 'blastp_score.p')
    return blastp_dict


def import_busco(busco_out_dir, output_dir, d_path, logger):
    '''Import BUSCO output'''
    # import_busco.py -b <busco_dir> -o <output_dir>
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    command = '{} --busco_dir {} --output_dir {}'.format(
        d_path['import_busco'], busco_out_dir, gene_filtering_dir
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_import_busco')
    logger_txt.debug('[Wrapper] %s', command)
    os.system(command)
    logger_time.debug('DONE : wrapper_import_busco\n')

    # Return files
    busco_dict = os.path.join(gene_filtering_dir, 'busco_score.p')
    return busco_dict


def import_pfam(pfam_scan_out, nr_prot_mapping_file, d_path, logger):
    '''Import Pfam'''
    # import_pfam.py -p <pfam_scan_out_file> -n <nr_prot_mapping>
    command = '{} --pfam_scan_out_file {} --nr_prot_mapping {}'.format(
        d_path['import_pfam'], pfam_scan_out, nr_prot_mapping_file
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_import_pfam')
    logger_txt.debug('[Wrapper] %s', command)
    os.system(command)
    logger_time.debug('DONE : wrapper_import_pfam\n')

    # Return files
    pfam_scan_out_dir = os.path.dirname(pfam_scan_out)
    pfam_dict = os.path.join(pfam_scan_out_dir, 'pfam_score.p')
    return pfam_dict


def import_blastn(blastn_out_files, output_dir, d_path, logger):
    '''Import BLASTn output'''
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    # import_blastn.py -b <blastn_out_files> -o <output_dir>
    command = '{} --blastn_out_files {} --output_dir {}'.format(
        d_path['import_blastn'], ' '.join(blastn_out_files), gene_filtering_dir
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_import_blastn')
    logger_txt.debug('[Wrapper] %s', command)
    os.system(command)
    logger_time.debug('DONE : wrapper_import_blastn\n')
    blastn_dict = os.path.join(gene_filtering_dir, 'blastn_score.p')
    return blastn_dict


def catch_bad_genes(gff3_files, genome_assembly, output_dir, d_path, logger):
    '''Catch bad genes'''
    # catch_bad_genes.py -g <gff3_files> -a <genome_assembly> -o <output_dir>
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    command = (
        '{} --gff3_files {} --genome_assembly {} --output_dir {}'.format(
            d_path['catch_bad_genes'], ' '.join(gff3_files), genome_assembly,
            gene_filtering_dir
        )
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_catch_bad_genes')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_catch_bad_genes\n')

    bad_dict = os.path.join(gene_filtering_dir, 'd_bad.p')
    return bad_dict


def filter_gff3s(
        genome_assembly, gff3_files, blastp_dict, busco_dict, pfam_dict,
        blastn_dict, bad_dict, nr_prot_file, nr_prot_mapping_file, output_dir,
        d_path, logger):
    '''Filter GFF3s'''
    # filter_gff3s.py -a <genome_assembly> -i <input_gff3s> -m <mapping_file>
    # -b <blastp_dict> -B <busco_dict> -p <pfam_dict> -N <blastn_dict>
    # -g <bad_dict> -n <nr_prot_file> -o <output_dir> -l <log_dir>
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    log_dir = os.path.join(output_dir, 'logs')
    command = (
        '{} --genome_assembly {} --input_gff3s {} --mapping_file {} '
        '--blastp_dict {} --busco_dict {} --pfam_dict {} --blastn_dict {} '
        '--bad_dict {} --nr_prot_file {} --output_dir {} --log_dir {}'
    ).format(
        d_path['filter_gff3s'], genome_assembly, ' '.join(gff3_files),
        nr_prot_mapping_file, blastp_dict, busco_dict, pfam_dict, blastn_dict,
        bad_dict, nr_prot_file, gene_filtering_dir, log_dir
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_filter_gff3s')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_filter_gff3s\n')


def gff3_postprocess(genome_assembly, output_dir, d_path, logger):
    '''GFF3 post-processing'''
    # gff3_postprocess.py -g <genome_assembly> -i <input_gff3> -o <output_gff3>
    input_gff3 = os.path.join(output_dir, 'gene_filtering', 'filtered_1.gff3')
    output_gff3 = os.path.join(output_dir, 'gene_filtering', 'filtered_2.gff3')
    command = (
        '{} --genome_assembly {} --input_gff3 {} --output_gff3 {}'
    ).format(
        d_path['gff3_postprocess'], genome_assembly, input_gff3, output_gff3
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_gff3_postprocess')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_gff3_postprocess\n')


def copy_output(output_dir, d_path, logger):
    '''Copy output'''
    # copy_output.py -o <output_dir>
    command = '{} --output_dir {}'.format(
        d_path['copy_output'], output_dir
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_copy_output')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE: wrapper_copy_output\n')


def create_markdown(
        genome_assembly, output_dir, trans_bams, trinity_asms, d_path, logger):
    '''Create markdown output summary'''
    # create_markdown.py -f <input_fasta> -g <input_gff3>
    # -t <trinity_assembly> -b <bam_file> -o <output_dir>
    fungap_gff3 = os.path.join(output_dir, 'gene_filtering/filtered_2.gff3')
    trans_bam = trans_bams[0]
    trinity_asm = trinity_asms[0]
    markdown_out_dir = os.path.join(output_dir, 'fungap_out')

    command = (
        '{} --input_fasta {} --input_gff3 {} --trinity_assembly {} '
        '--bam_file {} --output_dir {}'
    ).format(
        d_path['create_markdown'], genome_assembly, fungap_gff3, trinity_asm,
        trans_bam, markdown_out_dir
    )
    logger_time, logger_txt = logger
    logger_time.debug('START: wrapper_create_markdown')
    logger_txt.debug('[Wrapper] %s', command)
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE: wrapper_create_markdown\n')

    logger_time.debug('## DONE: FunGAP ##')


if __name__ == '__main__':
    main()
