FROM continuumio/miniconda3:4.9.2

ENV DEBIAN_FRONTEND=noninteractive

# Install mamba and dependencies
RUN apt-get update && apt-get install -y postgresql-contrib \
    && conda install mamba -n base -c conda-forge \
    && mamba install -c bioconda -c conda-forge \
        braker2=2.1.5 \
        trinity=2.12.0 \
        repeatmodeler=2.0.1 \
        hisat2=2.2.1 \
        pfam_scan=1.6 \
        busco=5.1.2 \
        augustus=3.4.0 \
    && cpanm Hash::Merge Logger::Simple Parallel::ForkManager YAML \
    && pip install biopython bcbio-gff markdown2 matplotlib

# Install Maker using Mamba (Maker installation is conflict with Busco)
RUN mamba create -c bioconda -c conda-forge -n maker maker=3.01.03

WORKDIR /workspace/FunGAP
ENV PFAM_DB=/workspace/FunGAP/db/pfam/
RUN git clone https://github.com/CompSynBioLab-KoreaUniv/FunGAP.git . \
    && mkdir -p $PFAM_DB \
    && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz -P $PFAM_DB \
    && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz -P $PFAM_DB \
    && gunzip $PFAM_DB/Pfam-A.hmm.gz \
    && gunzip $PFAM_DB/Pfam-A.hmm.dat.gz \
    && hmmpress $PFAM_DB/Pfam-A.hmm

RUN cd $(dirname $(which RepeatMasker))/../share/RepeatMasker \
    && echo "\n2\n/opt/conda/bin\n\n5\n" > tmp \
    && ./configure < tmp
