'''Import config function'''

import os
import sys


def import_config():
    '''Import config'''
    this_path = os.path.realpath(__file__)
    this_dir = os.path.dirname(this_path)
    config_file = os.path.join(this_dir, 'fungap.conf')
    with open(config_file) as f_in:
        config_txt = list(line.rstrip() for line in f_in)

    d_conf = {}
    for line in config_txt:
        if line.startswith('#'):
            continue
        line_split = line.split('=')
        if not line_split[1]:
            sys.exit(
                '[ERROR] There is a problem with fungap.conf Please re-run '
                'set_dependencies.py'
            )
        d_conf[line_split[0]] = line_split[1]
    return d_conf
