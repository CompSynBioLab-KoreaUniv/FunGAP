FROM continuumio/miniconda3:4.8.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y postgresql-contrib \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install braker2=2.1.5 trinity=2.11.0 repeatmodeler=2.0.1 hisat2=2.2.0 pfam_scan=1.6 busco=4.1.2 \
    && pip install biopython bcbio-gff networkx markdown2 matplotlib \
    && cpanm Hash::Merge Logger::Simple Parallel::ForkManager YAML

WORKDIR /workspace/FunGAP
ENV PFAM_DB=/workspace/FunGAP/db/pfam/
RUN git clone -b v1.1.0 https://github.com/CompSynBioLab-KoreaUniv/FunGAP.git . \
    && mkdir -p $PFAM_DB \
    && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -P $PFAM_DB \
    && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz -P $PFAM_DB \
    && gunzip $PFAM_DB/Pfam-A.hmm.gz \
    && gunzip $PFAM_DB/Pfam-A.hmm.dat.gz \
    && hmmpress $PFAM_DB/Pfam-A.hmm

ENV TERM=xterm
WORKDIR /opt/conda/share/RepeatMasker
RUN echo "\n2\n/opt/conda/bin\n\n5\n" > tmp \
    && perl ./configure < tmp

RUN conda create -n maker maker=2.31.10

WORKDIR /workspace/FunGAP/db/busco
RUN wget https://busco-archive.ezlab.org/v3/datasets/fungi_odb9.tar.gz \
    && wget https://busco-archive.ezlab.org/v3/datasets/ascomycota_odb9.tar.gz \
    && wget https://busco-archive.ezlab.org/v3/datasets/basidiomycota_odb9.tar.gz \
    && tar -zxvf fungi_odb9.tar.gz \
    && tar -zxvf ascomycota_odb9.tar.gz \
    && tar -zxvf basidiomycota_odb9.tar.gz
