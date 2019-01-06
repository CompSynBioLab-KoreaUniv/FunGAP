#!/usr/bin/env python2

'''
Run Maker

This scripts runs Maker four times with iterative SNAP gene model training,
    including running GeneMark.

The parameters set by FunGAP are:
    [Maker run1]
    split_hit=5000
    single_exon=1
    single_length=150
    correct_est_fusion=1
    est=$OUTPUT_DIR/trans_trinity/trinity_rnaseq/Trinity_rnaseq.fasta
    est2genome=1
    model_org= # Remove default "all" value
    rmlib=$OUTPUT_DIR/repeat_modeler/RM_[number].[date]/consensi.fa.classified

    [Maker run2]
    split_hit=5000
    single_exon=1
    single_length=150
    correct_est_fusion=1
    model_org=  # Remove default "all" value
    repeat_protein=  # Remove default "te_proteins.fasta" value
    snaphmm=$OUTPUT_DIR/gpre_maker/rnaseq/maker_run2/snp_training/snap_hmm_v1.hmm
    maker_gff=$OUTPUT_DIR/gpre_maker/rnaseq/maker_run1/genome_assembly.all.gff
    est_pass=1 # To save computing time
    protein_pass=1
    rm_pass=1

    [Maker run3]
    split_hit=5000
    single_exon=1
    single_length=150
    correct_est_fusion=1
    model_org=  # Remove default "all" value
    repeat_protein=  # Remove default "te_proteins.fasta" value
    snaphmm=$OUTPUT_DIR/gpre_maker/rnaseq/maker_run2/snp_training/snap_hmm_v2.hmm
    maker_gff=$OUTPUT_DIR/gpre_maker/rnaseq/maker_run2/genome_assembly.all.gff
    est_pass=1
    protein_pass=1
    rm_pass=1

    [Maker run4]
    split_hit=5000
    single_exon=1
    single_length=150
    correct_est_fusion=1
    model_org=  # Remove default "all" value
    repeat_protein=  # Remove default "te_proteins.fasta" value
    snaphmm=$OUTPUT_DIR/gpre_maker/rnaseq/maker_run2/snp_training/snap_hmm_v2.hmm
    maker_gff=$OUTPUT_DIR/gpre_maker/rnaseq/maker_run2/genome_assembly.all.gff
    est_pass=1
    protein_pass=1
    rm_pass=1
    keep_preds=1
    augustus_species=augustus_species
    gmhmm=$OUTPUT_DIR/gpre_genemark/output/gmhmm.mod
'''

# Import modules
import sys
import os
import re
from shutil import copyfile
from glob import glob
from argparse import ArgumentParser

# Get Logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)
sys.path.append(this_dir)
from set_logging import set_logging
from import_config import import_config

# Parameters
D_conf = import_config(this_dir)
program_name = 'maker'


def main(argv):
    optparse_usage = (
        'run_maker.py -i <input_fasta> -p <protein_db_fastas> -c <num_cores> '
        '-R <repeat_model> -e <est_files>'
    )
    parser = ArgumentParser(usage=optparse_usage)
    parser.add_argument(
        '-i', '--input_fasta', nargs=1, required=True,
        help='Input genome sequence in FASTA format'
    )
    parser.add_argument(
        "-a", "--augustus_species", nargs=1, required=True,
        help='"augustus --species=help" would be helpful'
    )
    parser.add_argument(
        "-p", "--protein_db_fasta", nargs='+', required=True,
        help="Protein db in FASTA foramt"
    )
    parser.add_argument(
        '-R', '--repeat_model', nargs=1, required=True,
        help="De novo repeat model by RepeatModeler: consensi.fa.classified"
    )
    parser.add_argument(
        '-e', '--est_files', nargs='+', required=True,
        help="Multiple EST data if available"
    )
    parser.add_argument(
        "-o", "--output_dir", nargs='?', default='maker_out',
        help="Output directory"
    )
    parser.add_argument(
        "-c", "--num_cores", nargs='?', default=1, type=int,
        help="Number of cores to be used"
    )
    parser.add_argument(
        "-l", "--log_dir", nargs='?', default='logs',
        help="Log directory"
    )
    parser.add_argument(
        '--gmes_fungus', action='store_true',
        help='--fungus flag in GeneMark'
    )

    args = parser.parse_args()
    input_fasta = os.path.abspath(args.input_fasta[0])
    output_dir = os.path.abspath(args.output_dir)
    log_dir = os.path.abspath(args.log_dir)
    augustus_species = args.augustus_species[0]
    protein_db_fastas = [os.path.abspath(x) for x in args.protein_db_fasta]
    num_cores = args.num_cores
    repeat_model = os.path.abspath(args.repeat_model[0])
    est_files = [os.path.abspath(x) for x in args.est_files]

    if args.gmes_fungus:
        gmes_fungus = '--fungus'
    else:
        gmes_fungus = ''

    # Create necessary directory
    create_dir(output_dir, log_dir)

    # Set logging
    log_file = os.path.join(log_dir, 'run_maker.log')
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run Maker on each EST file
    all_gff_file = ''
    for est_file in est_files:
        # Create directory
        est_prefix = os.path.basename(os.path.splitext(est_file)[0])
        est_prefix = est_prefix.replace('Trinity_', '')
        est_dir = os.path.join(output_dir, est_prefix)
        if not glob(est_dir):
            os.mkdir(est_dir)

        # Check maker is already done
        run_flag_run1 = check_maker_finished(
            output_dir, input_fasta, '1', est_prefix
        )

        # Run Maker batch
        logger_time.debug('START running Maker run1')
        if run_flag_run1:
            run_maker_batch(
                input_fasta, output_dir, log_dir, protein_db_fastas,
                num_cores, repeat_model, est_file, all_gff_file
            )
        else:
            logger_txt.debug('Running Maker has already been finished')
        logger_time.debug('DONE  running Maker run1')

        # Train run1 & run Maker run2
        all_gff_file_run1 = collect_result(
            input_fasta, output_dir, '1', est_prefix
        )
        logger_time.debug('START training run1 & running maker run2')
        snap_hmm_file_run1 = train_snap(
            output_dir, all_gff_file_run1, '1', est_prefix
        )
        run_flag_run2 = check_maker_finished(
            output_dir, input_fasta, '2', est_prefix
        )
        if run_flag_run2:
            run_maker_trained(
                input_fasta, output_dir, log_dir, augustus_species, num_cores,
                snap_hmm_file_run1, all_gff_file_run1, '2', est_prefix
            )
        else:
            logger_txt.debug('Running Maker has already been finished')
        logger_time.debug('DONE  training run1 & running maker run2')

        # Train run2 & run Maker run3
        all_gff_file_run2 = collect_result(
            input_fasta, output_dir, '2', est_prefix
        )
        logger_time.debug('START training run2 & running maker run3')
        snap_hmm_file_run2 = train_snap(
            output_dir, all_gff_file_run2, '2', est_prefix
        )
        run_flag_run3 = check_maker_finished(
            output_dir, input_fasta, '3', est_prefix
        )
        if run_flag_run3:
            run_maker_trained(
                input_fasta, output_dir, log_dir, augustus_species, num_cores,
                snap_hmm_file_run2, all_gff_file_run2, '3', est_prefix
            )
        else:
            logger_txt.debug('Running Maker has already been finished')
        logger_time.debug('DONE  training run2 & running maker run3')

        # Now, for final run, get masked assembly and get GeneMark hmm model
        masked_assembly = get_masked_asm(output_dir, est_files)

        # Run gmes or gmsn
        eukgmhmmfile = run_gmes(
            masked_assembly, num_cores, output_dir, log_dir, gmes_fungus
        )

        # Train run3 & run Maker run4
        all_gff_file_run3 = collect_result(
            input_fasta, output_dir, '3', est_prefix,
        )
        logger_time.debug('START training run3 & running maker run4')
        snap_hmm_file_run3 = train_snap(
            output_dir, all_gff_file_run3, '3', est_prefix
        )
        run_flag_run4 = check_maker_finished(
            output_dir, input_fasta, '4', est_prefix
        )
        if run_flag_run4:
            run_maker_trained(
                input_fasta, output_dir, log_dir, augustus_species, num_cores,
                snap_hmm_file_run3, all_gff_file_run3, '4', est_prefix,
                eukgmhmmfile
            )
        else:
            logger_txt.debug('Running Maker has already been finished')
        logger_time.debug('DONE  training run3 & running maker run4')

        # Get final GFF3 & FASTA
        collect_result_final(input_fasta, output_dir, est_prefix)
        all_gff_file = collect_result(input_fasta, output_dir, '4', est_prefix)


def import_file(input_file):
    with open(input_file) as f_in:
        txt = (line.rstrip() for line in f_in)
        txt = list(line for line in txt if line)
    return txt


def replace(fname, srcstr, deststr):
    f = open(fname)
    txt = f.read()
    txt = re.subn(r'\n{}.+'.format(srcstr), '\n{}'.format(deststr), txt)[0]
    f = open(fname, 'w')
    f.write(txt)
    f.close()


def create_dir(output_dir, log_dir):
    if not glob(output_dir):
        os.mkdir(output_dir)

    if not glob(log_dir):
        os.mkdir(log_dir)

    log_program_dir = os.path.join(log_dir, program_name)
    if not glob(log_program_dir):
        os.mkdir(log_program_dir)

    gmes_dir = os.path.join(output_dir, 'genemark_out')
    if not glob(gmes_dir):
        os.mkdir(gmes_dir)


def check_maker_finished(output_dir, input_fasta, version, prefix):
    # For first run
    index_log_file = glob(os.path.join(
        output_dir, prefix,
        'maker_run{}/*output/*master_datastore_index.log'.format(version)
    ))

    if not index_log_file:
        return True

    index_log = import_file(index_log_file[0])
    finished_scaffolds = []
    for line in index_log:
        line_split = line.split('\t')
        finish_tag = line_split[2]
        if finish_tag != 'FINISHED':
            continue

        finished_scaffold = line_split[0]
        finished_scaffolds.append(finished_scaffold)

    fasta = import_file(input_fasta)
    fasta_scaffolds = []
    for line in fasta:
        if not re.search('^>', line):
            continue
        fasta_scaffold = line.split(' ')[0].replace('>', '')
        fasta_scaffolds.append(fasta_scaffold)

    if finished_scaffolds == fasta_scaffolds:
        return False
    else:
        return True


def run_gmes(
    masked_assembly, num_cores, output_dir, log_dir, gmes_fungus
):
    genemark_bin = D_conf['GENEMARK_PATH']
    
    # Run gm_es.pl
    gmes_dir = os.path.join(output_dir, 'genemark_out')
    output_gmes = os.path.join(gmes_dir, 'output/gmhmm.mod')
    log_file = os.path.join(log_dir, program_name, 'gmes.log')

    logger_time.debug('START ruuning gmes to build hmm')
    if not glob(output_gmes):
        os.chdir(gmes_dir)
        command = (
            '{} --ES {} --cores {} --sequence {} --soft_mask 1 > '
            '{}'.format(
                genemark_bin, gmes_fungus, num_cores, masked_assembly,
                log_file
            )
        )
        logger_txt.debug('[Run] {}'.format(command))
        os.system(command)
    else:
        logger_txt.debug('GMES has already been finished')
    logger_time.debug('DONE  running gmes to build hmm')

    return output_gmes


def run_maker_batch(
    input_fasta, output_dir, log_dir, protein_db_fastas,
    num_cores, repeat_model, est_file, all_gff_file
):
    # Get binary
    maker_bin = D_conf['MAKER_PATH']

    est_prefix = os.path.basename(os.path.splitext(est_file)[0])
    est_prefix = est_prefix.replace('Trinity_', '')

    # Change directory
    maker_run1_dir = os.path.join(output_dir, est_prefix, 'maker_run1')
    if not glob(maker_run1_dir):
        os.mkdir(maker_run1_dir)

    # Change directory
    os.chdir(maker_run1_dir)

    # Make CTL files
    os.system('{} -CTL'.format(maker_bin))

    # Editting maker_opts.ctl - general
    replace('maker_opts.ctl', 'genome= ', 'genome={} '.format(input_fasta))
    replace(
        'maker_opts.ctl', 'protein=  ',
        'protein={} '.format(','.join(protein_db_fastas))
    )
    replace('maker_opts.ctl', 'cpus=1', 'cpus={}'.format(num_cores))
    replace('maker_opts.ctl', 'clean_up=0', 'clean_up=1'.format(num_cores))

    # For fungal genome
    replace('maker_opts.ctl', 'split_hit=', 'split_hit=5000')
    replace('maker_opts.ctl', 'single_exon=', 'single_exon=1')
    replace('maker_opts.ctl', 'single_length=', 'single_length=50')
    replace('maker_opts.ctl', 'correct_est_fusion=', 'correct_est_fusion=1')

    # If EST is provided
    if est_file != '':
        replace('maker_opts.ctl', "est= ", "est={} ".format(est_file))
        replace('maker_opts.ctl', "est2genome=0 ", "est2genome=1 ")

    # Set repeat model
    replace('maker_opts.ctl', 'model_org=all', 'model_org=')

    # Run faster feed aligned transcripts, proteins, repeat masking
    if all_gff_file:
        replace(
            'maker_opts.ctl', 'maker_gff= ',
            'maker_gff={} '.format(all_gff_file)
        )
        replace('maker_opts.ctl', 'protein_pass=0', 'protein_pass=1')
        replace('maker_opts.ctl', 'rm_pass=0', 'rm_pass=1')
        replace(
            'maker_opts.ctl', 'repeat_protein=', 'repeat_protein='
        )

    else:
        replace(
            'maker_opts.ctl', 'rmlib= ', 'rmlib={}'.format(repeat_model)
        )

    # Run maker
    maker_log = os.path.join(
        log_dir, program_name, 'maker_{}_run1.log'.format(est_prefix)
    )
    command = '{} -fix_nucleotides > {} 2>&1'.format(maker_bin, maker_log)
    logger_txt.debug('[Run] {}'.format(command))
    os.system(command)


def run_maker_trained(
    input_fasta, output_dir, log_dir, augustus_species, num_cores,
    snap_hmm_file, all_gff_file, version, prefix, eukgmhmmfile=None
):
    # Get binary
    maker_bin = D_conf['MAKER_PATH']

    # Create directory
    maker_run_dir = os.path.join(
        output_dir, prefix, 'maker_run{}'.format(version)
    )

    if not glob(maker_run_dir):
        os.mkdir(maker_run_dir)

    # Change directory
    os.chdir(maker_run_dir)

    # Make CTL files
    os.system('{} -CTL'.format(maker_bin))

    # Editting maker_opts.ctl - general
    replace('maker_opts.ctl', 'genome= ', 'genome={} '.format(input_fasta))
    replace('maker_opts.ctl', 'cpus=1', 'cpus={}'.format(num_cores))

    # For fungal genome
    replace('maker_opts.ctl', 'split_hit=', 'split_hit=5000')
    replace('maker_opts.ctl', 'single_exon=', 'single_exon=1')
    replace('maker_opts.ctl', 'single_length=', 'single_length=50')
    replace('maker_opts.ctl', 'correct_est_fusion=', 'correct_est_fusion=1')

    # Remove repeat org
    replace('maker_opts.ctl', 'model_org=all', 'model_org=')
    replace('maker_opts.ctl', 'repeat_protein=', 'repeat_protein=')

    # Supply SNAP HMM v1
    replace('maker_opts.ctl', 'snaphmm= ', 'snaphmm={} '.format(snap_hmm_file))

    # Run faster feed aligned transcripts, proteins, repeat masking
    replace(
        'maker_opts.ctl', 'maker_gff= ', 'maker_gff={} '.format(all_gff_file)
    )
    replace('maker_opts.ctl', 'est_pass=0', 'est_pass=1')
    replace('maker_opts.ctl', 'protein_pass=0', 'protein_pass=1')
    replace('maker_opts.ctl', 'rm_pass=0', 'rm_pass=1')

    # Last run, keep_preds=1
    if version == '4':
        replace('maker_opts.ctl', 'keep_preds=0', 'keep_preds=1')

        # Set AUGUSTUS species
        replace(
            'maker_opts.ctl', 'augustus_species= ',
            'augustus_species={} '.format(augustus_species)
        )

        # Set gmhmm
        replace('maker_opts.ctl', 'gmhmm= ', 'gmhmm={} '.format(eukgmhmmfile))
        replace(
            'maker_exe.ctl', 'gmhmme3= ',
            'gmhmme3={} '.format(D_conf['GMHMME3_PATH'])
        )
        replace(
            'maker_exe.ctl', 'probuild= ',
            'probuild={} '.format(D_conf['PROBUILD_PATH'])
        )

    # Run maker
    maker_log = os.path.join(
        log_dir, program_name, 'maker_{}_run{}.log'.format(prefix, version)
    )
    command = '{} -fix_nucleotides > {} 2>&1'.format(maker_bin, maker_log)
    logger_txt.debug('[Run] {}'.format(command))
    os.system(command)


def collect_result(
    input_fasta, output_dir, version, prefix
):
    maker_run_dir = os.path.join(
        output_dir, prefix, 'maker_run{}'.format(version)
    )
    input_prefix = (os.path.splitext(os.path.basename(input_fasta))[0])
    index_file = os.path.join(
        maker_run_dir,
        '{}.maker.output/{}_master_datastore_index.log'.format(
            input_prefix, input_prefix
        )
    )

    # Change directory to maker_run_dir
    gff3_merge_bin = D_conf['GFF3_MERGE_PATH']
    os.chdir(maker_run_dir)
    command = '{} -d {}'.format(gff3_merge_bin, index_file)
    logger_txt.debug('[Run] {}'.format(command))
    os.system(command)

    all_gff_file = '{}.all.gff'.format(input_prefix)
    all_gff_file_abs = os.path.abspath(all_gff_file)

    os.chdir(output_dir)

    return all_gff_file_abs


def collect_result_final(input_fasta, output_dir, prefix):
    maker_run_dir = os.path.join(output_dir, prefix, 'maker_run4')
    input_prefix = (os.path.splitext(os.path.basename(input_fasta))[0])
    index_file = os.path.join(
        maker_run_dir,
        '{}.maker.output/{}_master_datastore_index.log'.format(
            input_prefix, input_prefix
        )
    )

    # Change directory to maker_run_dir
    os.chdir(maker_run_dir)
    gff3_merge_bin = D_conf['GFF3_MERGE_PATH']
    command1 = '{} -g -n -d {}'.format(gff3_merge_bin, index_file)
    logger_txt.debug('[Run] {}'.format(command1))
    os.system(command1)

    # Collect FASTA, too
    fasta_merge_bin = D_conf['FASTA_MERGE_PATH']
    command2 = '{} -d {}'.format(fasta_merge_bin, index_file)
    logger_txt.debug('[Run] {}'.format(command2))
    os.system(command2)

    # Copy to maker root directory
    maker_root = os.path.join(output_dir, prefix)
    merged_gff3 = os.path.join(
        maker_root, 'maker_run4', '{}.all.gff'.format(input_prefix)
    )
    merged_faa = os.path.join(
        maker_root, 'maker_run4',
        '{}.all.maker.proteins.fasta'.format(input_prefix)
    )
    output_gff3 = os.path.join(maker_root, 'maker_{}.gff3'.format(prefix))
    output_faa = os.path.join(maker_root, 'maker_{}.faa'.format(prefix))

    copyfile(merged_gff3, output_gff3)
    copyfile(merged_faa, output_faa)

    os.chdir(output_dir)


def train_snap(output_dir, all_gff_file, version, prefix):
    maker_run_dir = os.path.join(
        output_dir, prefix, 'maker_run{}'.format(version)
    )

    maker2zff_bin = D_conf['MAKER2ZFF_PATH']
    fathom_bin = D_conf['FATHOM_PATH']
    forge_bin = D_conf['FORGE_PATH']
    hmm_assembler_bin = D_conf['HMM_ASSEMBLER_PATH']

    # Change directory into Maker run1 directory
    os.chdir(maker_run_dir)
    if not os.path.exists('snp_training'):
        os.makedirs('snp_training')
    os.chdir('snp_training')

    snap_hmm_file = os.path.abspath('snap_hmm_v{}.hmm'.format(version))
    if not os.path.exists(snap_hmm_file):
        # Run maker2zff to select a subset of gene models for training
        command1 = '{} -n {}'.format(maker2zff_bin, all_gff_file)
        logger_txt.debug('[Run] {}'.format(command1))
        os.system(command1)

        # It generates genome.dna and genome.ann
        # split the annotations into four categories: unique genes, warnings,
        # alternative spliced genes, overlapping genes, and errors
        command2 = '{} -categorize 1000 genome.ann genome.dna'.format(
            fathom_bin
        )
        logger_txt.debug('[Run] {}'.format(command2))
        os.system(command2)

        # Export the genes
        command3 = '{} -export 1000 -plus uni.ann uni.dna'.format(fathom_bin)
        logger_txt.debug('[Run] {}'.format(command3))
        os.system(command3)

        # Create directory
        if not os.path.exists('parameters'):
            os.makedirs('parameters')

        # Change directory
        os.chdir('parameters')

        # Generate the new parameters with forge
        command4 = '{} ../export.ann ../export.dna'.format(forge_bin)
        logger_txt.debug('[Run] {}'.format(command4))
        os.system(command4)

        # Generate the new HMM
        os.chdir('..')
        command5 = '{} snap_hmm_v{} parameters > snap_hmm_v{}.hmm'.format(
            hmm_assembler_bin, version, version
        )
        logger_txt.debug('[Run] {}'.format(command5))
        os.system(command5)
    else:
        logger_txt.debug("SNAP training has been alread finished for {}".format(
            os.path.basename(snap_hmm_file)))

    os.chdir(output_dir)

    return snap_hmm_file


def get_masked_asm(output_dir, est_files):
    est_prefix_first = (
        os.path.basename(os.path.splitext(est_files[0])[0])
        .replace('Trinity_', '')
    )

    maker_run_dir = os.path.join(
        output_dir, est_prefix_first, 'maker_run3'
    )
    # masked_asm_files = glob(masked_asm_path)
    masked_asm = os.path.join(output_dir, 'masked_assembly.fasta')
    command = 'find {} -name "query.masked.fasta" | xargs cat > {}'.format(
        maker_run_dir, masked_asm
    )
    logger_txt.debug('[Run] {}'.format(command))
    os.system(command)

    return masked_asm


if __name__ == "__main__":
    main(sys.argv[1:])
