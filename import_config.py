import os


def import_config(dir):
    config_file = os.path.join(dir, 'fungap.conf')
    with open(config_file) as f_in:
        config_txt = list(line.rstrip('\n') for line in f_in)

    D_conf = {}
    for line in config_txt:
        if line.startswith('#'):
            continue
        line_split = line.split('=')
        if not line_split[1]:
            sys.exit(
                '[ERROR] There is a problem with fungap.conf Please re-run '
                'set_dependencies.py'
            )
        D_conf[line_split[0]] = line_split[1]
    return D_conf