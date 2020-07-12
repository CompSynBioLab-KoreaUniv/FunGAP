FROM continuumio/miniconda:4.6.14

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y build-essential postgresql postgresql-contrib

RUN conda config --add channels bioconda/label/cf201901 \
    && conda config --add channels conda-forge/label/cf201901 \
    && conda install augustus rmblast maker hisat2 braker busco=3.0.2 blast pfam_scan \
    && pip install biopython bcbio-gff networkx markdown2 matplotlib \
    && cpanm Hash::Merge Logger::Simple Parallel::ForkManager YAML

ENV FUNGAP_DIR=/workspace/FunGAP

WORKDIR /workspace
RUN git clone https://github.com/CompSynBioLab-KoreaUniv/FunGAP.git \
    && cd FunGAP/ \
    && mkdir -p db/pfam \
    && cd db/pfam \
    && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz \
    && wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz \
    && gunzip Pfam-A.hmm.gz \
    && gunzip Pfam-A.hmm.dat.gz \
    && hmmpress Pfam-A.hmm

# BuscoDB download
RUN cd $FUNGAP_DIR \
    && mkdir -p db/busco \
    && cd db/busco \
    && wget https://busco-archive.ezlab.org/v3/datasets/fungi_odb9.tar.gz \
    && wget https://busco-archive.ezlab.org/v3/datasets/ascomycota_odb9.tar.gz \
    && wget https://busco-archive.ezlab.org/v3/datasets/basidiomycota_odb9.tar.gz \
    && tar -zxvf fungi_odb9.tar.gz \
    && tar -zxvf ascomycota_odb9.tar.gz \
    && tar -zxvf basidiomycota_odb9.tar.gz

# Install GeneMark
COPY gm_et_linux_64.tar.gz .
COPY gm_key_64.gz .

RUN mkdir $FUNGAP_DIR/external/ \
    && mv gm_et_linux_64.tar.gz gm_key_64.gz $FUNGAP_DIR/external/  \
    && cd $FUNGAP_DIR/external/ \
    && tar -zxvf gm_et_linux_64.tar.gz \
    && gunzip gm_key_64.gz \
    && cp gm_key_64 ~/.gm_key \
    && cd $FUNGAP_DIR/external/gm_et_linux_64/ \
    && cp other/reformat_fasta.pl . \
    && perl change_path_in_perl_scripts.pl "/usr/bin/env perl"

# Install RECON
RUN cd $FUNGAP_DIR/external/ \
    && wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz \
    && tar -zxvf RECON-1.08.tar.gz \
    && cd RECON-1.08/src/ \
    && make \
    && make install

# Install RepeatScout 1.0.5
RUN cd $FUNGAP_DIR/external/ \
    && wget http://www.repeatmasker.org/RepeatScout-1.0.5.tar.gz \
    && tar -zxvf RepeatScout-1.0.5.tar.gz \
    && cd RepeatScout-1 \
    && make

# Install NSEG
RUN cd $FUNGAP_DIR/external/ \
    && mkdir nseg \
    && cd nseg \
    && wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/genwin.c \
    && wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/genwin.h \
    && wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/lnfac.h \
    && wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/makefile \
    && wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/nmerge.c \
    && wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/nseg.c \
    && wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/runnseg \
    && make

ENV TERM xterm

# # Install RepeatMasker 4.0.8
RUN cd $FUNGAP_DIR/external/ \
  && wget http://www.repeatmasker.org/RepeatMasker-open-4-0-8.tar.gz \
  && tar -zxvf RepeatMasker-open-4-0-8.tar.gz

RUN cd $FUNGAP_DIR/external/RepeatMasker \
  && echo -e "\n/opt/conda/bin/perl\n$FUNGAP_DIR/external/RepeatMasker\n/opt/conda/bin/trf\n2\n/opt/conda/bin\nY\n5\n" > tmp \
  && perl ./configure < tmp

# Install RepeatModeler 1.0.11
RUN cd $FUNGAP_DIR/external/ \
  && wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-open-1.0.11.tar.gz \
  && tar -zxvf RepeatModeler-open-1.0.11.tar.gz \
  && cd RepeatModeler-open-1.0.11 \
  && echo -e "\n/opt/conda/bin/perl\n$FUNGAP_DIR/external/RepeatModeler-open-1.0.11\n$FUNGAP_DIR/external/RepeatMasker\n$FUNGAP_DIR/external/RECON-1.08/bin\n$FUNGAP_DIR/external/RepeatScout-1\n$FUNGAP_DIR/external/nseg\n/opt/conda/bin\n1\n/opt/conda/bin\nY\n3\n" > tmp \
  && perl ./configure < tmp \
  && cd ..

# Add fungap.conf
ADD https://gist.githubusercontent.com/lmtani/d37343a40e143b59336e4606055d1723/raw/fungap.conf \
    $FUNGAP_DIR/

##########
## Trinity
ENV TRINITY_VERSION="2.8.5"
ENV TRINITY_CO="d35f3c1149bab077ca7c83f209627784469c41c6"


RUN apt-get update && apt-get install -y cmake build-essential gcc g++ bowtie2 jellyfish default-jre curl libdb-dev zlib1g-dev bzip2 libncurses5-dev \
  && cd $FUNGAP_DIR/external \
  && git clone https://github.com/trinityrnaseq/trinityrnaseq.git \
  && cd trinityrnaseq \
  && git checkout $TRINITY_CO \
  && make && make plugins

## Jellyfish
RUN cd $FUNGAP_DIR/external \
  && wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.7/jellyfish-2.2.7.tar.gz \
  && tar xvf jellyfish-2.2.7.tar.gz \
  && cd jellyfish-2.2.7/ \
  && ./configure \
  && make

## Salmon
RUN cd $FUNGAP_DIR/external \
  && wget https://github.com/COMBINE-lab/salmon/releases/download/v0.9.1/Salmon-0.9.1_linux_x86_64.tar.gz \
  && tar xvf Salmon-0.9.1_linux_x86_64.tar.gz

ENV PATH=${PATH}:$FUNGAP_DIR/external/trinityrnaseq:$FUNGAP_DIR/external/Salmon-latest_linux_x86_64/bin/:$FUNGAP_DIR/external/jellyfish-2.2.7/bin/


# Need to enter container, configure RepeadMask and RepeatModeler manually.

#python /workspace/FunGAP/fungap.py \
#    --output_dir fungap_out \
#    --trans_read_1 sscita_1.fastq \
#    --trans_read_2 sscita_2.fastq \
#    --genome_assembly genome/pilon.fasta  \
#    --augustus_species ustilago_maydis  \
#    --sister_proteome sister_prot/prot_db.faa  \
#    --num_cores 10
