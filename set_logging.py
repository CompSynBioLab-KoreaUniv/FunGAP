import logging


# Define functions
def set_logging(log_file):
    # Set logger_time
    logger_time = logging.getLogger('logger_time')
    logger_time.setLevel(logging.DEBUG)

    fh_time = logging.FileHandler(log_file)
    fh_time.setLevel(logging.DEBUG)

    ch_time = logging.StreamHandler()
    ch_time.setLevel(logging.DEBUG)

    # create formatter and add it to the handlers
    formatter_time = logging.Formatter(
        '[%(asctime)s] %(message)s', datefmt='%m-%d %H:%M'
    )
    fh_time.setFormatter(formatter_time)
    ch_time.setFormatter(formatter_time)

    # add the handlers to the logger
    logger_time.addHandler(fh_time)
    logger_time.addHandler(ch_time)

    # Set logger_txt
    logger_txt = logging.getLogger('logger_txt')
    logger_txt.setLevel(logging.DEBUG)

    fh_txt = logging.FileHandler(log_file)
    fh_txt.setLevel(logging.DEBUG)

    ch_txt = logging.StreamHandler()
    ch_txt.setLevel(logging.DEBUG)

    # create formatter and add it to the handlers
    formatter_txt = logging.Formatter('%(message)s')
    fh_txt.setFormatter(formatter_txt)
    ch_txt.setFormatter(formatter_txt)

    # add the handlers to the logger
    logger_txt.addHandler(fh_txt)
    logger_txt.addHandler(ch_txt)

    return logger_time, logger_txt
