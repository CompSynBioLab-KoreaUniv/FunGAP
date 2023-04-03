#!/usr/bin/env python3

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

Output: the locations will be written in 'fungap.conf' file
Last updated: Aug 12, 2020
'''

import os
import subprocess
import sys
from argparse import ArgumentParser
from distutils import spawn


def main():
    '''Main function'''
    argparse_usage = (
        'check_dependencies.py -p <pfam_db_path> -g <genemark_path> '
        '-m <maker_path>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        '-p', '--pfam_db_path', nargs=1, required=True,
        help='Pfam database path'
    )
    parser.add_argument(
        '-g', '--genemark_path', nargs=1, required=True,
        help='GeneMark bin path'
    )
    parser.add_argument(
        '-m', '--maker_path', nargs=1, required=True,
        help='Maker bin path'
    )
    parser.add_argument(
        '-s', '--snap_path', nargs=1, required=True,
        help='SNAP-HMM bin path'
    )
    parser.add_argument(
        '-r', '--with_repeat_modeler', nargs='?', default='',
        help='User-defined RepeatModeler bin path'
    )
    parser.add_argument(
        '-a', '--with_augustus', nargs='?', default='',
        help='User-defined Augustus bin path (it should be Augustus >=3.4.0)'
    )
    parser.add_argument(
        '-H', '--with_hisat2', nargs='?', default='',
        help='User-defined Hisat2 bin path'
    )
    parser.add_argument(
        '-t', '--with_trinity', nargs='?', default='',
        help='User-defined Trinity bin path'
    )
    parser.add_argument(
        '-b', '--with_braker', nargs='?', default='',
        help='User-defined Braker2 bin path'
    )
    parser.add_argument(
        '-B', '--with_busco', nargs='?', default='',
        help='User-defined BUSCO bin path'
    )
    parser.add_argument(
        '-i', '--with_pfam_scan', nargs='?', default='',
        help='User-defined pfam_scan bin path'
    )

    args = parser.parse_args()
    pfam_db_path = os.path.abspath(args.pfam_db_path[0])
    i_genemark_path = os.path.abspath(args.genemark_path[0])
    i_maker_path = os.path.abspath(args.maker_path[0])
    i_snap_path = os.path.abspath(args.snap_path[0])
    if args.with_repeat_modeler:
        with_repeat_modeler = os.path.abspath(args.with_repeat_modeler)
    else:
        with_repeat_modeler = ''
    if args.with_augustus:
        with_augustus = os.path.abspath(args.with_augustus)
    else:
        with_augustus = ''
    with_hisat2 = os.path.abspath(args.with_hisat2) if args.with_hisat2 else ''
    if args.with_trinity:
        with_trinity = os.path.abspath(args.with_trinity)
    else:
        with_trinity = ''
    with_busco = os.path.abspath(args.with_busco) if args.with_busco else ''
    with_braker = os.path.abspath(args.with_braker) if args.with_braker else ''
    if args.with_pfam_scan:
        with_pfam_scan = os.path.abspath(args.with_pfam_scan)
    else:
        with_pfam_scan = ''
    pfam_db_path = check_db(pfam_db_path)
    (
        genemark_path, gmhmme3_path, probuild_path, build_database_path,
        repeat_modeler_path, hisat2_path, trinity_path, maker_path,
        gff3_merge_path, fasta_merge_path, maker2zff_path, fathom_path,
        forge_path, hmm_assembler_path, braker_path, busco_path,
        pfam_scan_path, blastp_path, blastn_path, blastx_path,
        makeblastdb_path, samtools_path, bamtools_path, augustus_path
    ) = get_path(
        i_genemark_path, i_maker_path, i_snap_path, with_repeat_modeler,
        with_augustus, with_hisat2, with_trinity, with_braker, with_busco,
        with_pfam_scan
    )
    check_working(
        genemark_path, gmhmme3_path, probuild_path, build_database_path,
        repeat_modeler_path, hisat2_path, trinity_path, maker_path,
        gff3_merge_path, fasta_merge_path, maker2zff_path, fathom_path,
        forge_path, hmm_assembler_path, braker_path, busco_path,
        pfam_scan_path, blastp_path, blastn_path, blastx_path, makeblastdb_path,
        samtools_path, bamtools_path, augustus_path
    )

    write_config(
        pfam_db_path, genemark_path, gmhmme3_path, probuild_path,
        build_database_path, repeat_modeler_path, hisat2_path, trinity_path,
        maker_path, gff3_merge_path, fasta_merge_path, maker2zff_path,
        fathom_path, forge_path, hmm_assembler_path, braker_path, busco_path,
        pfam_scan_path, blastp_path, blastn_path, blastx_path, makeblastdb_path,
        samtools_path, bamtools_path, augustus_path
    )
    print(
        '\nSetting dependencies is complete. Check fungap.conf file in the '
        'installation directory'
    )


def check_db(pfam_db_path):
    '''Check pfam_db_path'''
    if not os.path.isdir(pfam_db_path):
        sys.exit(
            '[ERROR] Pfam directory does not exist: {}'.format(pfam_db_path)
        )

    file_names = set(os.listdir(pfam_db_path))
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
            'be generated by "hmmpress Pfam-A.hmm" command'.format(pfam_db_path)
        )

    return pfam_db_path


def get_path(
        i_genemark_path, i_maker_path, i_snap_path, with_repeat_modeler,
        with_augustus, with_hisat2, with_trinity, with_braker, with_busco,
        with_pfam_scan):
    '''Get path'''
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

    genemark_path = check_binary('GeneMark', i_genemark_path, 'gmes_petap.pl')
    gmhmme3_path = check_binary('GeneMark', i_genemark_path, 'gmhmme3')
    probuild_path = check_binary('GeneMark', i_genemark_path, 'probuild')
    hisat2_path = check_binary('Hisat2', with_hisat2, 'hisat2')
    trinity_path = check_binary('Trinity', with_trinity, 'Trinity')
    maker_path = check_binary('Maker', i_maker_path, 'maker')
    gff3_merge_path = check_binary('Maker', i_maker_path, 'gff3_merge')
    fasta_merge_path = check_binary('Maker', i_maker_path, 'fasta_merge')
    maker2zff_path = check_binary('Maker', i_maker_path, 'maker2zff')
    fathom_path = check_binary('Snap', i_snap_path, 'fathom')
    forge_path = check_binary('Snap', i_snap_path, 'forge')
    hmm_assembler_path = check_binary('Snap', i_snap_path, 'hmm-assembler.pl')
    build_database_path = check_binary(
        'RepeatModeler (BuildDatabase)', with_repeat_modeler, 'BuildDatabase',
    )
    repeat_modeler_path = check_binary(
        'RepeatModeler (RepeatModeler)', with_repeat_modeler, 'RepeatModeler',
    )
    augustus_path = check_binary('Augustus', with_augustus, 'augustus')
    braker_path = check_binary('Braker', with_braker, 'braker.pl')
    busco_path = check_binary('BUSCO', with_busco, 'busco')
    pfam_scan_path = check_binary('Pfam_scan', with_pfam_scan, 'pfam_scan.pl')
    blastp_path = check_binary('BLASTp', '', 'blastp')
    blastn_path = check_binary('BLASTn', '', 'blastn')
    blastx_path = check_binary('BLASTx', '', 'blastx')
    makeblastdb_path = check_binary('MAKEBLASTDB', '', 'makeblastdb')
    samtools_path = check_binary('Samtools', '', 'samtools')
    bamtools_path = check_binary('Bamtools', '', 'bamtools')

    return (
        genemark_path, gmhmme3_path, probuild_path, build_database_path,
        repeat_modeler_path, hisat2_path, trinity_path, maker_path,
        gff3_merge_path, fasta_merge_path, maker2zff_path, fathom_path,
        forge_path, hmm_assembler_path, braker_path, busco_path,
        pfam_scan_path, blastp_path, blastn_path, blastx_path,
        makeblastdb_path, samtools_path, bamtools_path, augustus_path
    )


def check_working(
        genemark_path, gmhmme3_path, probuild_path, build_database_path,
        repeat_modeler_path, hisat2_path, trinity_path, maker_path,
        gff3_merge_path, fasta_merge_path, maker2zff_path, fathom_path,
        forge_path, hmm_assembler_path, braker_path, busco_path, pfam_scan_path,
        blastp_path, blastn_path, blastx_path, makeblastdb_path, samtools_path,
        bamtools_path, augustus_path):
    '''Check if programs work properly'''
    print('\n** Checking the dependencies if they properly work **\n')

    def check_working_internal(binary_path, command_list):
        try:
            subprocess.call(
                command_list, stdout=open(os.devnull, 'wb'),
                stderr=open(os.devnull, 'wb')
            )
        except subprocess.CalledProcessError:
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
    check_working_internal(braker_path, [braker_path, '--help'])
    check_working_internal(busco_path, [busco_path, '--help'])
    check_working_internal(pfam_scan_path, [pfam_scan_path, '-h'])
    check_working_internal(blastp_path, [blastp_path, '--help'])
    check_working_internal(blastn_path, [blastn_path, '--help'])
    check_working_internal(blastx_path, [blastx_path, '--help'])
    check_working_internal(makeblastdb_path, [makeblastdb_path, '--help'])
    check_working_internal(samtools_path, [samtools_path, '--help'])
    check_working_internal(bamtools_path, [bamtools_path, '--help'])
    check_working_internal(augustus_path, [augustus_path, '--help'])
    check_augustus_version(augustus_path)

    # For GeneMark, check the .gm_key
    home_path = os.path.expanduser('~')
    if not os.path.exists(os.path.join(home_path, '.gm_key')):
        sys.exit(
            '\n[ERROR] You do not have .gm_key in your home directory.\n'
            'Check https://wiki.gacrc.uga.edu/wiki/GeneMark'
        )

def check_augustus_version(augustus_path):
    '''Check Augustus version 3.4.0'''
    proc = subprocess.Popen(
        [augustus_path, '--version'], stderr=subprocess.PIPE
    )
    output = str(proc.stderr.read().decode('utf-8'))
    if not output.startswith('AUGUSTUS (3.4.0)'):
        sys.exit(
            '\n[ERROR] Augustus version is not 3.4.0. Check with "{} --version"'
            ''.format(augustus_path)
        )


def write_config(
        pfam_db_path, genemark_path, gmhmme3_path, probuild_path,
        build_database_path, repeat_modeler_path, hisat2_path, trinity_path,
        maker_path, gff3_merge_path, fasta_merge_path, maker2zff_path,
        fathom_path, forge_path, hmm_assembler_path, braker_path, busco_path,
        pfam_scan_path, blastp_path, blastn_path, blastx_path, makeblastdb_path,
        samtools_path, bamtools_path, augustus_path):
    '''Write config file'''
    this_path = os.path.realpath(__file__)
    this_dir = os.path.dirname(this_path)
    config_file = os.path.join(this_dir, 'fungap.conf')
    outhandle = open(config_file, 'w')

    outhandle.write('PFAM_DB_PATH={}\n'.format(pfam_db_path))
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
    outhandle.write('BRAKER_PATH={}\n'.format(braker_path))
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


if __name__ == '__main__':
    main()
