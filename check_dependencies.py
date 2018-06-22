#!/usr/bin/python

'''
Check if dependencies are correctly located and installed
    1) Hisat2
    2) Trinity
    3) Maker
    4) RepeatModeler
    5) Braker
    6) BUSCO
    7) InterProScan
    8) GeneMark
'''

# Import modules
import sys
import os
import subprocess
from glob import glob
from argparse import ArgumentParser

# Get logging
this_path = os.path.realpath(__file__)
this_dir = os.path.dirname(this_path)
sys.path.append(this_dir)
from set_logging import set_logging


# Main function
def main(argv):
    argparse_usage = (
        'check_dependencies.py -o <output_dir> -H <with_hisat2>'
        ' -t <with_trinity> -m <with_maker> -r <with_repeat_modeler>'
        ' -b <with_braker1>  -B <with_busco> -i <with_interproscan>'
        ' -g <with_genemark>'
    )
    parser = ArgumentParser(usage=argparse_usage)
    parser.add_argument(
        "-o", "--output_dir", dest="output_dir", nargs=1,
        help="Output directory"
    )
    parser.add_argument(
        "-H", "--with_hisat2", dest="with_hisat2", nargs='?',
        help="User-defined Hisat2 installation path (binary directory)"
    )
    parser.add_argument(
        "-t", "--with_trinity", dest="with_trinity", nargs='?',
        help="User-defined Trinity installation path (binary directory)"
    )
    parser.add_argument(
        "-m", "--with_maker", dest="with_maker", nargs='?',
        help="User-defined Maker installation path (binary directory)"
    )
    parser.add_argument(
        "-r", "--with_repeat_modeler", dest="with_repeat_modeler", nargs='?',
        help="User-defined Repeat Modeler installation path (binary directory)"
    )
    parser.add_argument(
        "-b", "--with_braker1", dest="with_braker1", nargs='?',
        help="User-defined Braker1 installation path (binary directory)"
    )
    parser.add_argument(
        "-B", "--with_busco", dest="with_busco", nargs='?',
        help="User-defined BUSCO installation path (binary directory)"
    )
    parser.add_argument(
        "-i", "--with_interproscan", dest="with_interproscan", nargs='?',
        help="User-defined InterproScan installation path (binary directory)"
    )
    parser.add_argument(
        "-g", "--with_genemark", dest="with_genemark", nargs='?',
        help="User-defined GeneMark installation path (binary directory)"
    )

    args = parser.parse_args()
    if args.output_dir:
        output_dir = os.path.abspath(args.output_dir[0])
    else:
        print '[ERROR] You should provide OUTPUT DIRECTORY'
        sys.exit(2)

    if args.with_hisat2:
        with_hisat2 = os.path.abspath(args.with_hisat2)
    else:
        with_hisat2 = ''

    if args.with_trinity:
        with_trinity = os.path.abspath(args.with_trinity)
    else:
        with_trinity = ''

    if args.with_maker:
        with_maker = os.path.abspath(args.with_maker)
    else:
        with_maker = ''

    if args.with_repeat_modeler:
        with_repeat_modeler = os.path.abspath(args.with_repeat_modeler)
    else:
        with_repeat_modeler = ''

    if args.with_braker1:
        with_braker1 = os.path.abspath(args.with_braker1)
    else:
        with_braker1 = ''

    if args.with_busco:
        with_busco = os.path.abspath(args.with_busco)
    else:
        with_busco = ''

    if args.with_interproscan:
        with_interproscan = os.path.abspath(args.with_interproscan)
    else:
        with_interproscan = ''

    if args.with_genemark:
        with_genemark = os.path.abspath(args.with_genemark)
    else:
        with_genemark = ''

    # Create necessary dirs
    create_dir(output_dir)

    # Set logging
    log_dir = os.path.join(output_dir, 'logs')
    log_file = os.path.join(
        log_dir, 'pipeline', 'check_dependencies.log')
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    # Run functions :) Slow is as good as Fast
    logger_time.debug('Check dependencies: get paths')
    (
        hisat2_path, trinity_path, maker_path, repeat_modeler_path,
        braker1_path, busco_path, interproscan_path, genemark_path
    ) = get_path(
        with_hisat2, with_trinity, with_maker, with_repeat_modeler,
        with_braker1, with_busco, with_interproscan, with_genemark
    )

    logger_txt.debug('')
    logger_time.debug('Check dependencies: check tools working')
    check_working(
        hisat2_path, trinity_path, maker_path, repeat_modeler_path,
        braker1_path, busco_path, interproscan_path, genemark_path
    )

    write_config(
        output_dir, hisat2_path, trinity_path, maker_path,
        repeat_modeler_path, braker1_path, busco_path, interproscan_path,
        genemark_path
    )

    # Check BLAST installation
    check_blast()


def create_dir(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    log_dir = os.path.join(output_dir, 'logs')
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    log_pipeline_dir = os.path.join(log_dir, 'pipeline')
    if not os.path.exists(log_pipeline_dir):
        os.mkdir(log_pipeline_dir)


def get_path(
    with_hisat2, with_trinity, with_maker, with_repeat_modeler,
    with_braker1, with_busco, with_interproscan, with_genemark
):
    def check_binary(tool_name, path, binary, fungap_external):
        if path:
            if os.path.exists(os.path.join(path, binary)):
                logger_txt.debug(
                    '[%s path] %s' % (tool_name, os.path.join(path, binary))
                )
                return os.path.join(path, binary)
            else:
                logger_txt.debug(
                    "\n[ERROR] You provided wrong %s path. Please check" % (
                        path
                    )
                )
                sys.exit(2)
        elif os.path.exists(fungap_external):
            logger_txt.debug(
                '[%s path] %s' % (tool_name, fungap_external)
            )
            return fungap_external
        else:
            logger_txt.debug(
                '\n[ERROR] No installation for %s. We expect %s to exist' % (
                    tool_name, fungap_external
                )
            )
            sys.exit(2)

    fungap_hisat2 = os.path.join(this_dir, 'external/hisat2/hisat2')
    binary_hisat2 = 'hisat2'
    hisat2_path = check_binary(
        'Hisat2', with_hisat2, binary_hisat2, fungap_hisat2
    )

    fungap_trinity = os.path.join(this_dir, 'external/trinityrnaseq/Trinity')
    binary_trinity = 'Trinity'
    trinity_path = check_binary(
        'Trinity', with_trinity, binary_trinity, fungap_trinity
    )

    fungap_maker = os.path.join(this_dir, 'external/maker/bin/maker')
    binary_maker = 'maker'
    maker_path = check_binary(
        'Maker', with_maker, binary_maker, fungap_maker
    )

    fungap_repeat_modeler = os.path.join(
        this_dir, 'external/RepeatModeler/RepeatModeler'
    )
    binary_repeat_modeler = 'RepeatModeler'
    repeat_modeler_path = check_binary(
        'RepeatModeler', with_repeat_modeler, binary_repeat_modeler,
        fungap_repeat_modeler
    )

    fungap_braker1 = os.path.join(
        this_dir, 'external/BRAKER1/braker.pl'
    )
    binary_braker1 = 'braker.pl'
    braker1_path = check_binary(
        'Braker1', with_braker1, binary_braker1, fungap_braker1
    )

    fungap_busco = os.path.join(
        this_dir, 'external/BUSCO_v1.1b1/BUSCO_v1.1b1.py'
    )
    binary_busco = 'BUSCO_v1.1b1.py'
    busco_path = check_binary(
        'Busco', with_busco, binary_busco, fungap_busco
    )

    glob_interproscan = glob(os.path.join(
        this_dir, 'external/interproscan-5.*-*.0/interproscan.sh'
    ))
    if glob_interproscan:
        fungap_interproscan = glob_interproscan[0]
    else:
        fungap_interproscan = os.path.join(
            this_dir, 'external/interproscan-5.18-57.0/interproscan.sh'
        )
    binary_interproscan = 'interproscan.sh'
    interproscan_path = check_binary(
        'InterProScan', with_interproscan, binary_interproscan,
        fungap_interproscan
    )

    fungap_genemark = os.path.join(
        this_dir, 'external/gm_et_linux_64/gmes_petap/gmes_petap.pl'
    )
    binary_genemark = 'gmes_petap.pl'
    genemark_path = check_binary(
        'GeneMark', with_genemark, binary_genemark, fungap_genemark
    )

    return (
        hisat2_path, trinity_path, maker_path, repeat_modeler_path,
        braker1_path, busco_path, interproscan_path, genemark_path
    )


def check_working(
    hisat2_path, trinity_path, maker_path, repeat_modeler_path,
    braker1_path, busco_path, interproscan_path, genemark_path
):
    def check_working_internal(binary_path):
        try:
            subprocess.call(
                [binary_path, "--help"], stdout=open(os.devnull, 'wb')
            )
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                logger_txt.debug(
                    '%s --help is not working. Please check your installation'
                ) % (os.path.basename(binary_path))
                sys.exit(2)
            else:
                logger_txt.debug(
                    '%s --help is not working. Please check your installation'
                ) % (os.path.basename(binary_path))
                sys.exit(2)
                raise
        logger_txt.debug(
            '[%s] Running is OK' % (os.path.basename(binary_path))
        )

    check_working_internal(hisat2_path)
    check_working_internal(trinity_path)
    check_working_internal(maker_path)
    check_working_internal(repeat_modeler_path)
    check_working_internal(braker1_path)
    check_working_internal(busco_path)
    check_working_internal(interproscan_path)

    # For GeneMark, we should check .gm_key
    home_dir = os.path.expanduser('~')
    if not os.path.exists(os.path.join(home_dir, '.gm_key')):
        logger_txt.debug(
            '\n[ERROR] You do not have .gm_key in you home directory.\n'
            'Check https://wiki.gacrc.uga.edu/wiki/GeneMark'
        )
        sys.exit(2)

    try:
        subprocess.call(
            [genemark_path], stdout=open(os.devnull, 'wb')
        )
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            logger_txt.debug(
                '%s --help is not working. Please check your installation'
            ) % (os.path.basename(genemark_path))
            sys.exit(2)
        else:
            logger_txt.debug(
                '%s --help is not working. Please check your installation'
            ) % (os.path.basename(genemark_path))
            sys.exit(2)
            raise
    logger_txt.debug(
        '[%s] Running is OK' % (os.path.basename(genemark_path))
    )


def write_config(
    output_dir, hisat2_path, trinity_path, maker_path,
    repeat_modeler_path, braker1_path, busco_path, interproscan_path,
    genemark_path
):
    config_file = os.path.join(output_dir, 'fungap_exe.config')
    outhandle = open(config_file, 'w')
    outhandle.write('# Program paths\n')
    outhandle.write('HISAT2_PATH=%s\n' % (hisat2_path))
    outhandle.write('TRINITY_PATH=%s\n' % (trinity_path))
    outhandle.write('MAKER_PATH=%s\n' % (maker_path))
    outhandle.write('REPEATMODELER_PATH=%s\n' % (repeat_modeler_path))
    outhandle.write('BRAKER1_PATH=%s\n' % (braker1_path))
    outhandle.write('BUSCO_PATH=%s\n' % (busco_path))
    outhandle.write('INTERPROSCAN_PATH=%s\n' % (interproscan_path))
    outhandle.write('GENEMARK_PATH=%s\n' % (genemark_path))
    outhandle.close()


def check_blast():
    binary_paths = ['blastp', 'blastx', 'blastn', 'makeblastdb']
    for binary_path in binary_paths:
        try:
            subprocess.call(
                [binary_path, "-help"], stdout=open(os.devnull, 'wb')
            )
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                logger_txt.debug(
                    '%s -help is not working. Please check your installation'
                ) % (os.path.basename(binary_paths))
                sys.exit(2)
            else:
                logger_txt.debug(
                    '%s -help is not working. Please check your installation'
                ) % (os.path.basename(binary_paths))
                sys.exit(2)
                raise
        logger_txt.debug(
            '[%s] Running is OK' % (os.path.basename(binary_path))
        )


if __name__ == "__main__":
    main(sys.argv[1:])
