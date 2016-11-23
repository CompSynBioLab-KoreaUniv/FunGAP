#!/usr/bin/python

'''
Wrapper gene prediction pipeline
Author Byoungnam Min on Aug 13, 2015
'''

# Import modules
import sys
import os
import re
import shlex
from glob import glob
from datetime import datetime
from subprocess import check_call
from argparse import ArgumentParser

# Get Logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)
sys.path.append(this_dir)
from set_logging import set_logging

# File paths
run_hisat2_path = os.path.join(this_dir, 'run_hisat2.py')
run_trinity_path = os.path.join(this_dir, 'run_trinity.py')
run_repeat_modeler_path = os.path.join(this_dir, 'run_repeat_modeler.py')

run_augustus_path = os.path.join(this_dir, 'run_augustus.py')
run_maker_path = os.path.join(this_dir, 'run_maker.py')
run_braker1_path = os.path.join(this_dir, 'run_braker1.py')

run_busco_path = os.path.join(this_dir, 'run_busco.py')
run_iprscan_path = os.path.join(this_dir, 'run_interproscan_pfam.py')
make_nr_prot_path = os.path.join(this_dir, 'make_nr_prot.py')
run_blastp_path = os.path.join(this_dir, 'run_blastp.py')
import_blast_path = os.path.join(this_dir, 'import_blast.py')
import_busco_path = os.path.join(this_dir, 'import_busco.py')
import_pfam_path = os.path.join(this_dir, 'import_pfam.py')
filter_gff3s_path = os.path.join(this_dir, 'filter_gff3s.py')
catch_bad_genes_path = os.path.join(this_dir, 'catch_bad_genes.py')


# Main function
def main(argv):
    argparse_usage = (
        'fgap.py -g <genome_assembly> -r <trans_read_files> '
        '-o <output_dir> -p <project_name> -a <augustus_species> '
        '-O <org_id> -s <sister_proteome>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-o", "--output_dir", dest="output_dir", nargs=1,
        help="Output directory"
    )
    parser.add_argument(
        "-r", "--trans_read_files", dest="trans_read_files", nargs='+',
        help="Multiple transcriptome read files in FASTAQ"
    )
    parser.add_argument(
        "-p", "--project_name", dest="project_name", nargs=1,
        help="Project name without space. e.g. Mag, Eco, Pst_LUM"
    )
    parser.add_argument(
        "-g", "--genome_assembly", dest="genome_assembly", nargs=1,
        help="Genome assembly file in FASTA format"
    )
    parser.add_argument(
        "-a", "--augustus_species", dest="augustus_species", nargs=1,
        help="AUGUSTUS species"
    )
    parser.add_argument(
        "-O", "--org_id", dest="org_id", nargs=1,
        help="Organism ID. E.g. Hypma for Hypsizygus marmoreus"
    )
    parser.add_argument(
        "-s", "--sister_proteome", dest="sister_proteome", nargs=1,
        help="Sister proteome sequences in .faa"
    )
    parser.add_argument(
        "-c", "--num_cores", dest="num_cores", nargs=1,
        help="Number of cores to be used"
    )

    args = parser.parse_args()
    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir[0])
    else:
        print '[ERROR] Please provide ROOT DIRECTORY'
        parser.print_help()
        sys.exit(2)

    if args.trans_read_files:
        trans_read_files = [os.path.abspath(x) for x in args.trans_read_files]
    else:
        print '[ERROR] Please provide transcriptome read files'
        parser.print_help()
        sys.exit(2)

    if args.project_name:
        project_name = args.project_name[0]
    else:
        print '[ERROR] Please provide PROJECTN NAME'
        parser.print_help()
        sys.exit(2)

    if args.genome_assembly:
        genome_assembly = os.path.abspath(args.genome_assembly[0])
    else:
        print '[ERROR] Please provide transcriptome read files'
        parser.print_help()
        sys.exit(2)

    if args.augustus_species:
        augustus_species = args.augustus_species[0]
    else:
        print '[ERROR] Please provide transcriptome AUGUSTUS SPECIES'
        parser.print_help()
        sys.exit(2)

    if args.org_id:
        org_id = args.org_id[0]
    else:
        print '[ERROR] Please provide transcriptome ORGANISM ID'
        parser.print_help()
        sys.exit(2)

    if args.sister_proteome:
        sister_proteome = os.path.abspath(args.sister_proteome[0])
    else:
        print '[ERROR] Please provide TRANSCRIPTOME READ FILES'
        parser.print_help()
        sys.exit(2)

    if args.num_cores:
        num_cores = args.num_cores[0]
    else:
        print '[ERROR] Please provide NUMBER OF CORES'
        parser.print_help()
        sys.exit(2)

    # Create nessasary dirs
    create_dir(output_dir)

    # Set logging
    log_file = os.path.join(
        output_dir, 'logs', 'pipeline', 'fgap.log'
    )
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    logger_txt.debug("\n============ New Run %s ============" % (
        datetime.now())
    )

    # Run functions :) Slow is as good as Fast
    trans_bams = run_hisat2(
        genome_assembly, trans_read_files, output_dir, num_cores
    )
    trinity_asms = run_trinity(trans_bams, output_dir, project_name, num_cores)
    repeat_model_file = run_repeat_modeler(
        genome_assembly, output_dir, project_name, num_cores
    )
    maker_gff3s, maker_faas = run_maker(
        genome_assembly, output_dir, augustus_species,
        project_name, sister_proteome, num_cores, repeat_model_file,
        trinity_asms
    )
    # Get masked assembly
    masked_assembly = os.path.join(
        output_dir, 'gpre_maker', 'masked_assembly.fasta'
    )

    # Run Augustus
    augustus_gff3, augustus_faa = run_augustus(
        masked_assembly, output_dir, augustus_species
    )

    # Run Braker1
    braker1_gff3s, braker1_faas = run_braker1(
        masked_assembly, trans_bams, output_dir, num_cores
    )

    # Run BUSCO on each gene models
    if not glob(os.path.join(output_dir, 'gpre_busco')):
        os.mkdir(os.path.join(output_dir, 'gpre_busco'))

    for maker_faa in maker_faas:
        maker_prefix = os.path.basename(maker_faa).split('.')[0]
        maker_busco = os.path.join(output_dir, 'gpre_busco', maker_prefix)
        run_busco(maker_faa, maker_busco, num_cores)

    augustus_prefix = os.path.basename(augustus_faa).split('.')[0]
    augustus_busco = os.path.join(output_dir, 'gpre_busco', augustus_prefix)
    run_busco(augustus_faa, augustus_busco, num_cores)

    for braker1_faa in braker1_faas:
        braker1_prefix = os.path.basename(braker1_faa).split('.')[0]
        braker1_busco = os.path.join(output_dir, 'gpre_busco', braker1_prefix)
        run_busco(braker1_faa, braker1_busco, num_cores)

    busco_dir = os.path.join(output_dir, 'gpre_busco')

    # Get protein nr by removing identical proteins
    all_prot_files = maker_faas + [augustus_faa] + braker1_faas
    nr_prot_file, nr_prot_mapping_file = make_nr_prot(
        all_prot_files, output_dir
    )

    # Run BLASTp with nr prot file
    blastp_output = run_blastp(
        nr_prot_file, output_dir, sister_proteome, num_cores
    )

    # Run IPRscan with nr prot file
    ipr_output = run_iprscan(nr_prot_file, output_dir)

    # Import BLAST, BUSCO and Pfam score
    blast_dict_score, blast_dict_evalue = import_blast(
        blastp_output, nr_prot_mapping_file
    )
    busco_dict_score, busco_dict_list = import_busco(busco_dir)
    pfam_dict_score, pfam_dict_count = import_pfam(
        ipr_output, nr_prot_mapping_file
    )

    # Catch bad genes
    D_bad_pickle = catch_bad_genes(
        maker_gff3s, augustus_gff3, braker1_gff3s, genome_assembly,
        output_dir
    )

    filter_gff3s(
        maker_gff3s, augustus_gff3, braker1_gff3s,
        blast_dict_score, blast_dict_evalue, busco_dict_score, busco_dict_list,
        pfam_dict_score, pfam_dict_count, D_bad_pickle, nr_prot_file,
        nr_prot_mapping_file, org_id, output_dir
    )


def create_dir(output_dir):
    log_dir = os.path.join(output_dir, 'logs')
    if not glob(log_dir):
        os.mkdir(log_dir)

    log_pipeline_dir = os.path.join(output_dir, 'logs', 'pipeline')
    if not glob(log_pipeline_dir):
        os.mkdir(log_pipeline_dir)


def run_hisat2(genome_assembly, trans_read_files, output_dir, num_cores):
    hisat2_output_dir = os.path.join(output_dir, 'trans_hisat2')
    log_dir = os.path.join(output_dir, 'logs')

    # run_hisat2.py -r <fastq1> <fastq2> <fastq3> ... \
    # -o <output_dir> -l <log_dir> -f <ref_fasta> -c <num_cores>
    command = 'python %s -r %s -o %s -l %s -f %s -c %s' % (
        run_hisat2_path, ' '.join(trans_read_files), hisat2_output_dir,
        log_dir, genome_assembly, num_cores
    )
    logger_time.debug('START: wrapper_run_hisat2')
    logger_txt.debug('[Wrapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_hisat2')

    # Get output BAM file paths
    trans_bams = []
    for trans_read_file in trans_read_files:
        prefix = os.path.basename(trans_read_file).split('_')[0]
        hisat2_output = os.path.join(hisat2_output_dir, '%s.sam' % (prefix))
        sorted_bam_file = re.sub('.sam$', '_sorted.bam', hisat2_output)
        trans_bams.append(sorted_bam_file)

    trans_bams2 = list(set(trans_bams))

    return trans_bams2


def run_trinity(trans_bams, output_dir, project_name, num_cores):
    trinity_output_dir = os.path.join(output_dir, 'trans_trinity')
    log_dir = os.path.join(output_dir, 'logs')
    # run_trinity.py -b <bam_files> -o <output_dir> -l <log_dir>
    # -p <project_name> -c <num_cores>'
    command = 'python %s -b %s -o %s -l %s -p %s -c %s' % (
        run_trinity_path, ' '.join(trans_bams), trinity_output_dir,
        log_dir, project_name, num_cores
    )
    logger_time.debug('START: wrapper_run_trinity')
    logger_txt.debug('[Wrapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_trinity')

    # Get output transcriptome assembly files
    trinity_asms = glob(os.path.join(
        output_dir, 'trans_trinity', '*/Trinity_*.fasta')
    )
    return trinity_asms


def run_repeat_modeler(genome_assembly, output_dir, project_name, num_cores):
    # run_repeat_modeler.py -g <genome_assembly> -r <root_dir>
    # -p <project_name>
    rm_output_dir = os.path.join(output_dir, 'repeat_modeler')
    log_dir = os.path.join(output_dir, 'logs')
    command = 'python %s -g %s -o %s -l %s -p %s -c %s' % (
        run_repeat_modeler_path, genome_assembly, rm_output_dir, log_dir,
        project_name, num_cores
    )
    logger_time.debug('START: wrapper_run_repeat_modeler')
    logger_txt.debug('[Wrapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_repeat_modeler')

    repeat_model_file = glob(os.path.join(
        output_dir, 'repeat_modeler/RM*/consensi.fa.classified')
    )[0]
    return repeat_model_file


def run_maker(
    genome_assembly, output_dir, augustus_species, project_name,
    sister_proteome, num_cores, repeat_model_file, trinity_asms
):
    # run_maker.py -i <input_fasta> -r <root_dir> \
    # -p <project_name> -P <protein_db_fastas> -c <num_cores> \
    # -R <repeat_model> -e <est_files>'
    command = (
        'python %s -i %s -r %s -a %s -p %s -P %s -c %s -R %s -e %s'
    ) % (
        run_maker_path, genome_assembly, output_dir, augustus_species,
        project_name, sister_proteome, num_cores, repeat_model_file,
        ' '.join(trinity_asms)
    )
    logger_time.debug('START: wrapper_run_maker')
    logger_txt.debug('[Wrapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_maker')

    maker_gff3s = glob(
        os.path.join(output_dir, 'gpre_maker', '*/maker_*.gff3')
    )
    maker_faas = glob(os.path.join(output_dir, 'gpre_maker', '*/maker_*.faa'))
    return maker_gff3s, maker_faas


def run_augustus(masked_assembly, output_dir, augustus_species):
    # run_augustus.py -i <input_fasta> -o <output_dir> -s <species>
    output_dir = os.path.join(output_dir, 'gpre_augustus')
    log_dir = os.path.join(output_dir, 'logs')
    command = 'python %s -i %s -o %s -s %s -l %s' % (
        run_augustus_path, masked_assembly, output_dir, augustus_species,
        log_dir
    )
    logger_time.debug('START: wrapper_run_augustus')
    logger_txt.debug('[Wrapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_augustus')

    augustus_gff3 = os.path.join(output_dir, 'augustus.gff3')
    augustus_faa = os.path.join(output_dir, 'augustus.faa')
    return augustus_gff3, augustus_faa


def run_braker1(masked_assembly, trans_bams, output_dir, num_cores):
    braker1_output_dir = os.path.join(output_dir, 'gpre_braker1')
    log_dir = os.path.join(output_dir, 'logs')

    # run_braker1.py -m <masked_assembly> -b <bam_files>
    # -o <output_dir> -l <log_dir> -p <project_name> -c <num_cores>
    command = 'python %s -m %s -b %s -o %s -l %s -c %s' % (
        run_braker1_path, masked_assembly, ' '.join(trans_bams),
        braker1_output_dir, log_dir, num_cores
    )
    logger_time.debug('START: wrapper_run_braker1')
    logger_txt.debug('[Wrapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_braker1')

    prefixes = [os.path.basename(x).split('_')[0] for x in trans_bams]
    prefixes_u = list(set(prefixes))

    braker1_gff3s = []
    braker1_faas = []
    for prefix in prefixes_u:
        braker1_gff3 = os.path.join(
            output_dir, 'gpre_braker1',
            prefix, 'braker1_%s.gff3' % (prefix)
        )
        braker1_gff3s.append(braker1_gff3)
        braker1_faa = os.path.join(
            output_dir, 'gpre_braker1',
            prefix, 'braker1_%s.faa' % (prefix)
        )
        braker1_faas.append(braker1_faa)

    return braker1_gff3s, braker1_faas


def run_busco(input_faa, output_prefix, num_cores):
    output_dir = os.path.dirname(output_prefix)
    log_dir = os.path.join(output_dir, 'logs')
    # run_busco.py -i <input_fasta> -o <output_dir>
    # -l <log_dir> -c <num_cores>
    command = 'python %s -i %s -o %s -l %s -c %s' % (
        run_busco_path, input_faa, output_dir, log_dir, num_cores
    )
    logger_time.debug('START: wrapper_run_busco')
    logger_txt.debug('[Wrapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_busco')


def make_nr_prot(all_prot_files, output_dir):
    # make_nr_prot.py -i <faa_files> -o <output_dir>
    # Create directory
    gpre_filtered_dir = os.path.join(output_dir, 'gpre_filtered')
    if not glob(gpre_filtered_dir):
        os.mkdir(gpre_filtered_dir)

    command = 'python %s -i %s -o %s' % (
        make_nr_prot_path, ' '.join(all_prot_files), gpre_filtered_dir
    )
    logger_time.debug('START: wrapper_make_nr_prot')
    logger_txt.debug('[Wrapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_make_nr_prot')

    nr_prot_file = os.path.join(gpre_filtered_dir, 'nr_prot.faa')
    nr_prot_mapping_file = os.path.join(
        gpre_filtered_dir, 'nr_prot_mapping.txt'
    )

    return nr_prot_file, nr_prot_mapping_file


def run_blastp(nr_prot_file, output_dir, sister_proteome, num_cores):
    # Create directory
    blastp_dir = os.path.join(output_dir, 'gpre_filtered/blastp')
    if not glob(blastp_dir):
        os.mkdir(blastp_dir)

    # run_blastp.py -i <input_fasta> -r <references> -o <output> --nr
    blastp_output = os.path.join(
        output_dir, 'gpre_filtered', 'blastp', 'nr_prot.blast'
    )
    output_prefix = os.path.join(
        output_dir, 'gpre_filtered', 'blastp', 'nr_prot'
    )
    command = 'python %s -i %s -f %s -o %s -r %s -c %s' % (
        run_blastp_path, nr_prot_file,
        sister_proteome, output_prefix, output_dir, num_cores
    )
    logger_time.debug('START: wrapper_run_blastp')
    logger_txt.debug('[Wrapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_blastp')

    return blastp_output


def run_iprscan(nr_prot_file, output_dir):
    # run_interproscan_pfam.py -i <input_fasta> -o <output_dir> \
    # -l <log_dir> -c <num_cores>
    ipr_output_dir = os.path.join(output_dir, 'gpre_ipr')
    log_dir = os.path.join(output_dir, 'logs')
    command = 'python %s -i %s -l %s -o %s' % (
        run_iprscan_path, nr_prot_file, log_dir, ipr_output_dir
    )
    logger_time.debug('START: wrapper_run_iprscan')
    logger_txt.debug('[Wapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_iprscan')

    ipr_output = os.path.join(output_dir, 'gpre_ipr', 'nr_prot.xml')
    return ipr_output


def import_blast(blastp_output, nr_prot_mapping_file):
    # import_blast.py -b <blast_file> -m <nr_prot_mapping>
    # -o <output_prefix>

    blastp_output_dir = os.path.dirname(blastp_output)
    blastp_output_prefix = os.path.basename(blastp_output).split('.')[0]
    output_prefix = os.path.join(blastp_output_dir, blastp_output_prefix)
    command = 'python %s -b %s -m %s -o %s' % (
        import_blast_path, blastp_output, nr_prot_mapping_file, output_prefix
    )
    logger_time.debug('START: wrapper_import_blast')
    logger_txt.debug('[Wrapper] %s' % (command))
    os.system(command)
    logger_time.debug('DONE : wrapper_import_blast')

    # Return files
    blast_dict_score = '%s_blast_score.p' % (output_prefix)
    blast_dict_evalue = '%s_blast_evalue.p' % (output_prefix)

    return blast_dict_score, blast_dict_evalue


def import_busco(busco_dir):
    # import_busco.py -b <busco_dir> -o <output_prefix>

    command = 'python %s -b %s' % (import_busco_path, busco_dir)
    logger_time.debug('START: wrapper_import_busco')
    logger_txt.debug('[Wrapper] %s' % (command))
    os.system(command)
    logger_time.debug('DONE : wrapper_import_busco')

    # Return files
    busco_dict_score = os.path.join(busco_dir, 'busco_score.p')
    busco_dict_list = os.path.join(busco_dir, 'busco_list.p')

    return busco_dict_score, busco_dict_list


def import_pfam(ipr_output, nr_prot_mapping_file):
    # import_pfam.py -i <ipr_file> -m <nr_prot_mapping> \
    # -o <output_prefix>

    ipr_output_dir = os.path.dirname(ipr_output)
    ipr_output_prefix = os.path.basename(ipr_output).split('.')[0]
    output_prefix = os.path.join(ipr_output_dir, ipr_output_prefix)
    command = 'python %s -i %s -m %s -o %s' % (
        import_pfam_path, ipr_output, nr_prot_mapping_file, output_prefix
    )
    logger_time.debug('START: wrapper_import_pfam')
    logger_txt.debug('[Wrapper] %s' % (command))
    os.system(command)
    logger_time.debug('DONE : wrapper_import_pfam')

    # Return files
    pfam_dict_score = '%s_pfam_score.p' % (output_prefix)
    pfam_dict_count = '%s_pfam_count.p' % (output_prefix)

    return pfam_dict_score, pfam_dict_count


def catch_bad_genes(
    maker_gff3s, augustus_gff3, braker1_gff3s, genome_assembly,
    output_dir
):
    bad_output_dir = os.path.join(output_dir, 'gpre_filtered')
    command = 'python %s -g %s %s %s -f %s -o %s' % (
        catch_bad_genes_path, ' '.join(maker_gff3s), augustus_gff3,
        ' '.join(braker1_gff3s), genome_assembly, bad_output_dir
    )
    logger_time.debug('START: wrapper_catch_bad_genes')
    logger_txt.debug('[Wrapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_filter_gff3s')

    return os.path.join(bad_output_dir, 'D_bad.p')


def filter_gff3s(
    maker_gff3s, augustus_gff3, braker1_gff3s,
    blast_dict_score, blast_dict_evalue, busco_dict_score, busco_dict_list,
    pfam_dict_score, pfam_dict_count, D_bad_pickle, nr_prot_file,
    nr_prot_mapping_file, org_id, output_dir
):

    # filter_gff3s.py -i <input_gff3s> -m <mapping_file> -b <blast_dict>
    # -B <busco_dict> -p <ipr_dict> -g <bad_dict> -n <nr_prot_file>
    # -s <short_id> -o <output_prefix> -r <root_dir>

    # Run filter_gff3s
    output_prefix = os.path.join(output_dir, 'gpre_filtered', 'gpre_filtered')
    command = (
        'python %s -i %s -m %s -b %s %s -B %s %s -p %s %s -g %s -n %s '
        '-s %s -o %s'
    ) % (
        filter_gff3s_path,
        ' '.join(maker_gff3s + [augustus_gff3] + braker1_gff3s),
        nr_prot_mapping_file, blast_dict_score, blast_dict_evalue,
        busco_dict_score, busco_dict_list, pfam_dict_score, pfam_dict_count,
        D_bad_pickle, nr_prot_file, org_id, output_prefix
    )
    logger_time.debug('START: wrapper_filter_gff3s')
    logger_txt.debug('[Wrapper] %s' % (command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_filter_gff3s')


if __name__ == "__main__":
    main(sys.argv[1:])
