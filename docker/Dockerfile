FROM taniguti/fungap-base:v1.1.1

ENV FUNGAP_DIR=/workspace/FunGAP

# Install GeneMark
COPY gmes_linux_64.tar.gz /workspace/FunGAP/external/
COPY gm_key_64.gz /workspace/FunGAP/external/


WORKDIR /workspace/FunGAP/external/
RUN tar -zxvf gmes_linux_64_4.tar.gz \
    && gunzip gm_key_64.gz \
    && cp gm_key_64 ~/.gm_key \
    && cd $FUNGAP_DIR/external/gmes_linux_64_4/ \
    && cp other/reformat_fasta.pl . \
    && perl change_path_in_perl_scripts.pl "/usr/bin/env perl"

COPY fungap.conf $FUNGAP_DIR/
