#!/usr/bin/env python2

'''
Check if dependencies are correctly located and installed
 1) GeneMark
 2) RepeatModeler
 3) Hisat2
 4) Trinity
 5) Maker
 6) Braker
 7) BUSCO
 8) Pfam Scan
 9) Blast
 10) Samtools

The locations will be written in 'fungap.conf' file
'''

# Import modules
import os
import sys
import subprocess
from distutils import spawn
from argparse import ArgumentParser

# Get logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)


# Main function
def main(argv):
    argparse_usage = (
        'check_dependencies.py -p <pfam_db> -u <busco_dir> -g <genemark_dir>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-p", "--pfam_db_dir", required=True,
        help="Path of Pfam database directory"
    )
    parser.add_argument(
        "-u", "--busco_db_dir", required=True,
        help="Path of BUSCO database directory"
    )
    parser.add_argument(
        "-g", "--genemark_dir", required=True,
        help="GeneMark installation path (binary directory)"
    )
    parser.add_argument(
        "-r", "--repeat_modeler_dir", required=True,
        help="Repeat Modeler installation path (binary directory)"
    )
    parser.add_argument(
        "-H", "--with_hisat2", nargs='?', default='',
        help="User-defined Hisat2 installation path (binary directory)"
    )
    parser.add_argument(
        "-t", "--with_trinity", nargs='?', default='',
        help="User-defined Trinity installation path (binary directory)"
    )
    parser.add_argument(
        "-m", "--with_maker", nargs='?', default='',
        help="User-defined Maker installation path (binary directory)"
    )
    parser.add_argument(
        "-b", "--with_braker1", nargs='?', default='',
        help="User-defined Braker1 installation path (binary directory)"
    )
    parser.add_argument(
        "-B", "--with_busco", nargs='?', default='',
        help="User-defined BUSCO installation path (binary directory)"
    )
    parser.add_argument(
        "-i", "--with_pfam_scan", nargs='?', default='',
        help="User-defined pfam_scan installation path (binary directory)"
    )

    args = parser.parse_args()
    pfam_db_dir = os.path.abspath(args.pfam_db_dir)
    busco_db_dir = os.path.abspath(args.busco_db_dir)
    genemark_dir = os.path.abspath(args.genemark_dir)
    repeat_modeler_dir = os.path.abspath(args.repeat_modeler_dir)
    with_hisat2 = args.with_hisat2
    with_trinity = args.with_trinity
    with_maker = args.with_maker
    with_braker1 = args.with_braker1
    with_busco = args.with_busco
    with_pfam_scan = args.with_pfam_scan

    pfam_db_path, busco_db_path = check_db(pfam_db_dir, busco_db_dir)
    (
        genemark_path, gmhmme3_path, probuild_path, build_database_path,
        repeat_modeler_path, hisat2_path, trinity_path, maker_path,
        gff3_merge_path, fasta_merge_path, maker2zff_path, fathom_path,
        forge_path, hmm_assembler_path, braker1_path, busco_path,
        pfam_scan_path, blastp_path, blastn_path, blastx_path,
        makeblastdb_path, samtools_path, bamtools_path, augustus_path
    ) = get_path(
        genemark_dir, repeat_modeler_dir, with_hisat2, with_trinity,
        with_maker, with_braker1, with_busco, with_pfam_scan
    )
    check_working(
        genemark_path, gmhmme3_path, probuild_path, build_database_path,
        repeat_modeler_path, hisat2_path, trinity_path, maker_path,
        gff3_merge_path, fasta_merge_path, maker2zff_path, fathom_path,
        forge_path, hmm_assembler_path, braker1_path, busco_path,
        pfam_scan_path, blastp_path, blastn_path, blastx_path, makeblastdb_path,
        samtools_path, bamtools_path, augustus_path
    )

    write_config(
        pfam_db_path, busco_db_path, genemark_path, gmhmme3_path, probuild_path,
        build_database_path, repeat_modeler_path, hisat2_path, trinity_path,
        maker_path, gff3_merge_path, fasta_merge_path, maker2zff_path,
        fathom_path, forge_path, hmm_assembler_path, braker1_path, busco_path,
        pfam_scan_path, blastp_path, blastn_path, blastx_path, makeblastdb_path,
        samtools_path, bamtools_path, augustus_path
    )
    print(
        '\nSetting dependencies is complete. Check fungap.conf file in the '
        'installation directory'
    )


def check_db(pfam_db_dir, busco_db_dir):
    # Check pfam_db_dir
    if not os.path.isdir(pfam_db_dir):
        sys.exit(
            '[ERROR] Pfam directory does not exist: {}'.format(pfam_db_dir)
        )

    file_names = set(os.listdir(pfam_db_dir))
    pfam_files = [
        'Pfam-A.hmm.h3f', 'Pfam-A.hmm.h3i', 'Pfam-A.hmm.h3m', 'Pfam-A.hmm.h3p',
        'Pfam-A.hmm', 'Pfam-A.hmm.dat'
    ]
    pfam_flag = True
    for pfam_file in pfam_files:
        if pfam_file not in file_names:
            pfam_flag = False

    if pfam_flag is False:
        sys.exit(
            '[ERROR] Pfam files are not found. Six files (Pfam-A.hmm  '
            'Pfam-A.hmm.dat Pfam-A.hmm.h3f Pfam-A.hmm.h3i Pfam-A.hmm.h3m '
            'Pfam-A.hmm.h3p) should be located at {} directory. These files can'
            'be generated by "hmmpress Pfam-A.hmm" command'.format(pfam_db_dir)
        )

    # Check BUSCO directory
    if not os.path.isdir(busco_db_dir):
        sys.exit(
            '[ERROR] BUSCO DB directory does not exist: {}'.format(busco_db_dir)
        )

    file_names2 = set(os.listdir(busco_db_dir))
    if 'lengths_cutoff' not in file_names2:
        sys.exit(
            '[ERROR] It looks you did not provide correct BUSCO DB directory'
        )

    return pfam_db_dir, busco_db_dir


def get_path(
    genemark_dir, repeat_modeler_dir, with_hisat2, with_trinity, with_maker,
    with_braker1, with_busco, with_pfam_scan
):
    print('\n** Checking the installed locations of dependencies **\n')

    def check_binary(tool_name, path, binary):
        if path:
            binary_abs_path = os.path.join(path, binary)
            if os.path.exists(binary_abs_path):
                print(
                    '[{} path] {}'.format(tool_name, binary_abs_path)
                )
                return binary_abs_path
            else:
                sys.exit(
                    '[ERROR] No installation of {} found. "{}" should be l'
                    'ocated in the {} directory'.format(tool_name, binary, path)
                )
        elif spawn.find_executable(binary):
            binary_abs_path = spawn.find_executable(binary)
            print(
                '[{} path] {}'.format(tool_name, binary_abs_path)
            )
            return binary_abs_path
        else:
            sys.exit(
                '[ERROR] No installation of {} found. "{}" should be found '
                'in PATH environmental variable'.format(tool_name, binary)
            )

    genemark_path = check_binary('GeneMark', genemark_dir, 'gmes_petap.pl')
    gmhmme3_path = check_binary('GeneMark', genemark_dir, 'gmhmme3')
    probuild_path = check_binary('GeneMark', genemark_dir, 'probuild')
    hisat2_path = check_binary('Hisat2', with_hisat2, 'hisat2')
    trinity_path = check_binary('Trinity', with_trinity, 'Trinity')
    maker_path = check_binary('Maker', with_maker, 'maker')
    gff3_merge_path = check_binary('Maker', with_maker, 'gff3_merge')
    fasta_merge_path = check_binary('Maker', with_maker, 'fasta_merge')
    maker2zff_path = check_binary('Maker', with_maker, 'maker2zff')
    fathom_path = check_binary('Snap', '', 'fathom')
    forge_path = check_binary('Snap', '', 'forge')
    hmm_assembler_path = check_binary('Snap', '', 'hmm-assembler.pl')
    build_database_path = check_binary(
        'RepeatModeler (BuildDatabase)', repeat_modeler_dir, 'BuildDatabase',
    )
    repeat_modeler_path = check_binary(
        'RepeatModeler (RepeatModeler)', repeat_modeler_dir, 'RepeatModeler',
    )
    braker1_path = check_binary('Braker1', with_braker1, 'braker.pl')
    busco_path = check_binary('BUSCO', with_busco, 'run_busco')
    pfam_scan_path = check_binary('Pfam_scan', with_pfam_scan, 'pfam_scan.pl')
    blastp_path = check_binary('BLASTp', '', 'blastp')
    blastn_path = check_binary('BLASTn', '', 'blastn')
    blastx_path = check_binary('BLASTx', '', 'blastx')
    makeblastdb_path = check_binary('MAKEBLASTDB', '', 'makeblastdb')
    samtools_path = check_binary('Samtools', '', 'samtools')
    bamtools_path = check_binary('Bamtools', '', 'bamtools')
    augustus_path = check_binary('Augustus', '', 'augustus')

    return (
        genemark_path, gmhmme3_path, probuild_path, build_database_path,
        repeat_modeler_path, hisat2_path, trinity_path, maker_path,
        gff3_merge_path, fasta_merge_path, maker2zff_path, fathom_path,
        forge_path, hmm_assembler_path, braker1_path, busco_path,
        pfam_scan_path, blastp_path, blastn_path, blastx_path,
        makeblastdb_path, samtools_path, bamtools_path, augustus_path
    )


def check_working(
    genemark_path, gmhmme3_path, probuild_path, build_database_path,
    repeat_modeler_path, hisat2_path, trinity_path, maker_path, gff3_merge_path,
    fasta_merge_path, maker2zff_path, fathom_path, forge_path,
    hmm_assembler_path, braker1_path, busco_path, pfam_scan_path, blastp_path,
    blastn_path, blastx_path, makeblastdb_path, samtools_path, bamtools_path,
    augustus_path
):
    print('\n** Checking the dependencies if they properly work **\n')

    def check_working_internal(binary_path, command_list):
        try:
            subprocess.call(
                command_list, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb')
            )
        except:
            sys.exit(
                '\n[ERROR] {} --help is not working. Please check your '
                'installation'.format(os.path.basename(binary_path))
            )
        print(
            '[{}] Running is OK'.format(os.path.basename(binary_path))
        )

    check_working_internal(genemark_path, [genemark_path])
    check_working_internal(gmhmme3_path, [gmhmme3_path])
    check_working_internal(probuild_path, [probuild_path])
    check_working_internal(build_database_path, [build_database_path, '--help'])
    check_working_internal(repeat_modeler_path, [repeat_modeler_path, '--help'])
    check_working_internal(hisat2_path, [hisat2_path, '--help'])
    check_working_internal(trinity_path, [trinity_path, '--help'])
    check_working_internal(maker_path, [maker_path, '--help'])
    check_working_internal(gff3_merge_path, [gff3_merge_path, '--help'])
    check_working_internal(fasta_merge_path, [fasta_merge_path, '--help'])
    check_working_internal(maker2zff_path, [maker2zff_path, '--help'])
    check_working_internal(fathom_path, [fathom_path, '--help'])
    check_working_internal(forge_path, [forge_path, '--help'])
    check_working_internal(hmm_assembler_path, [hmm_assembler_path, '--help'])
    check_working_internal(braker1_path, [braker1_path, '--help'])
    check_working_internal(busco_path, [busco_path, '--help'])
    check_working_internal(pfam_scan_path, [pfam_scan_path, '-h'])
    check_working_internal(blastp_path, [blastp_path, '--help'])
    check_working_internal(blastn_path, [blastn_path, '--help'])
    check_working_internal(blastx_path, [blastx_path, '--help'])
    check_working_internal(makeblastdb_path, [makeblastdb_path, '--help'])
    check_working_internal(samtools_path, [samtools_path, '--help'])
    check_working_internal(bamtools_path, [bamtools_path, '--help'])
    check_working_internal(augustus_path, [augustus_path, '--help'])

    # For GeneMark, check the .gm_key
    home_dir = os.path.expanduser('~')
    if not os.path.exists(os.path.join(home_dir, '.gm_key')):
        sys.exit(
            '\n[ERROR] You do not have .gm_key in your home directory.\n'
            'Check https://wiki.gacrc.uga.edu/wiki/GeneMark'
        )


def write_config(
    pfam_db_path, busco_db_path, genemark_path, gmhmme3_path, probuild_path,
    build_database_path, repeat_modeler_path, hisat2_path,
    trinity_path, maker_path, gff3_merge_path, fasta_merge_path,
    maker2zff_path, fathom_path, forge_path, hmm_assembler_path,
    braker1_path, busco_path, pfam_scan_path,
    blastp_path, blastn_path, blastx_path, makeblastdb_path, samtools_path,
    bamtools_path, augustus_path
):

    config_file = os.path.join(this_dir, 'fungap.conf')
    outhandle = open(config_file, 'w')

    outhandle.write('PFAM_DB_PATH={}\n'.format(pfam_db_path))
    outhandle.write('BUSCO_DB_PATH={}\n'.format(busco_db_path))
    outhandle.write('GENEMARK_PATH={}\n'.format(genemark_path))
    outhandle.write('GMHMME3_PATH={}\n'.format(gmhmme3_path))
    outhandle.write('PROBUILD_PATH={}\n'.format(probuild_path))
    outhandle.write('BUILDDATABASE_PATH={}\n'.format(build_database_path))
    outhandle.write('REPEATMODELER_PATH={}\n'.format(repeat_modeler_path))
    outhandle.write('HISAT2_PATH={}\n'.format(hisat2_path))
    outhandle.write('TRINITY_PATH={}\n'.format(trinity_path))
    outhandle.write('MAKER_PATH={}\n'.format(maker_path))
    outhandle.write('GFF3_MERGE_PATH={}\n'.format(gff3_merge_path))
    outhandle.write('FASTA_MERGE_PATH={}\n'.format(fasta_merge_path))
    outhandle.write('MAKER2ZFF_PATH={}\n'.format(maker2zff_path))
    outhandle.write('FATHOM_PATH={}\n'.format(fathom_path))
    outhandle.write('FORGE_PATH={}\n'.format(forge_path))
    outhandle.write('HMM_ASSEMBLER_PATH={}\n'.format(hmm_assembler_path))
    outhandle.write('BRAKER1_PATH={}\n'.format(braker1_path))
    outhandle.write('BUSCO_PATH={}\n'.format(busco_path))
    outhandle.write('PFAM_SCAN_PATH={}\n'.format(pfam_scan_path))
    outhandle.write('BLASTP_PATH={}\n'.format(blastp_path))
    outhandle.write('BLASTN_PATH={}\n'.format(blastn_path))
    outhandle.write('BLASTX_PATH={}\n'.format(blastx_path))
    outhandle.write('MAKEBLASTDB_PATH={}\n'.format(makeblastdb_path))
    outhandle.write('SAMTOOLS_PATH={}\n'.format(samtools_path))
    outhandle.write('BAMTOOLS_PATH={}\n'.format(bamtools_path))
    outhandle.write('AUGUSTUS_PATH={}\n'.format(augustus_path))
    outhandle.close()


if __name__ == "__main__":
    main(sys.argv[1:])
