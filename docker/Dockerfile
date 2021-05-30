FROM taniguti/fungap-base:v1.1.1

ENV FUNGAP_DIR=/workspace/FunGAP

# Install GeneMark
COPY gmes_linux_64.tar.gz /workspace/FunGAP/external/
COPY gm_key_64.gz /workspace/FunGAP/external/


WORKDIR /workspace/FunGAP/external/
RUN tar -zxvf gmes_linux_64.tar.gz \
    && gunzip gm_key_64.gz \
    && cp gm_key_64 ~/.gm_key \
    && cd $FUNGAP_DIR/external/gmes_linux_64/ \
    && cp other/reformat_fasta.pl . \
    && perl change_path_in_perl_scripts.pl "/usr/bin/env perl"

COPY patch_braker.pl /opt/conda/bin/filterGenesIn_mRNAname.pl
COPY fungap.conf $FUNGAP_DIR/

#$FUNGAP_DIR/fungap.py --genome_assembly GCF_000146045.2_R64_genomic.fna --trans_read_1 SRR1198667_1.fastq --trans_read_2 SRR1198667_2.fastq --augustus_species saccharomyces_cerevisiae_S288C --busco_dataset ascomycota_odb10 --sister_proteome prot_db.faa --num_cores 8
