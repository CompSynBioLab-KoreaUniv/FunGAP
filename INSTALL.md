# Installation of FunGAP v1.1.0

Because FunGAP implements many dependent programs, you may encounter issues during
installation. Please don't hesitate to post on *Issues* or contact me (mbnmbn00@gmail.com) for help.

These steps were tested and confirmed in freshly installed Ubuntu 18.04 LTS.

## 0. FunGAP requirements

### 0.1. Required softwares (and tested versions)

1. [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) v2.1.0
1. [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) v2.6.6
1. [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) v1.0.11
1. [Maker](http://www.yandell-lab.org/software/maker.html) v2.31.10
1. [GeneMark-ES/ET](http://topaz.gatech.edu/GeneMark/license_download.cgi) v4.38
1. [Augustus](https://github.com/Gaius-Augustus/Augustus) v3.3
1. [Braker](http://exon.gatech.edu/braker1.html) v1.9
1. [BUSCO](https://busco.ezlab.org/) v3.0.2
1. [Pfam_scan](https://www.ebi.ac.uk/seqdb/confluence/display/THD/PfamScan) v1.6
1. [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) v2.6.0+
1. [Samtools](http://www.htslib.org/download/) v1.9
1. [Bamtools](https://github.com/pezmaster31/bamtools) v2.4.1

### 0.2. Required databases

1. [BUSCO](https://busco.ezlab.org/) odb9
1. [Pfam](https://pfam.xfam.org/) release 32.0

<br/>

## 1. Setup Anaconda environment

For robust installation, we recommend to use Anaconda environment and install dependent programs and libraries as much as possible in the environment.

### 1.1. Install Anaconda2 (v4.5.12 tested)

Download and install Anaconda2 (We assume that you install it in ```$HOME/anaconda2```)

```
cd $HOME
wget https://repo.continuum.io/archive/Anaconda2-2018.12-Linux-x86_64.sh
bash Anaconda2-2018.12-Linux-x86_64.sh
```

### 1.2. Set conda environment

```
echo ". ~/anaconda2/etc/profile.d/conda.sh" >> ~/.bashrc
source ~/.bashrc
```

### 1.3. Create and activate an environment

```
conda update conda
conda create -n fungap
conda activate fungap
```

### 1.4. Add channels

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### 1.5. Install dependencies

```
conda install -c bioconda augustus rmblast maker trinity hisat2 braker busco blast pfam_scan
pip install biopython bcbio-gff networkx markdown2 matplotlib
cpanm Hash::Merge Logger::Simple Parallel::ForkManager YAML
```

<br />

## 2. Download and install FunGAP

### 2.1. Download FunGAP

Download FunGAP using GitHub clone. Suppose we are installing FunGAP in your `$HOME` directory, but you are free to change the location. `$FUNGAP_DIR` is going to be your FunGAP installation directory.

```
cd $HOME
git clone https://github.com/CompSynBioLab-KoreaUniv/FunGAP.git
cd fungap/
export FUNGAP_DIR=$(pwd)
```

<br />

## 3. Download databases

Download Pfam and BUSCO databases in your `$FUNGAP_DIR/db` directory.

### 3.1. Pfam DB download 

ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release

```
cd $FUNGAP_DIR  # Change directory to FunGAP installation directory
mkdir -p db/pfam
cd db/pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
gunzip Pfam-A.hmm.gz
gunzip Pfam-A.hmm.dat.gz
hmmpress Pfam-A.hmm  # HMMER package (would be automatically installed in the above Anaconda step)
```

### 3.2. BUSCO DB download

There are various databases in BUSCO, so just download one of them fitted to your target genome. Here are example commands.

```
cd $FUNGAP_DIR
mkdir -p db/busco
cd db/busco
wget https://busco.ezlab.org/datasets/fungi_odb9.tar.gz
wget https://busco.ezlab.org/datasets/ascomycota_odb9.tar.gz
wget https://busco.ezlab.org/datasets/basidiomycota_odb9.tar.gz
tar -zxvf fungi_odb9.tar.gz
tar -zxvf ascomycota_odb9.tar.gz
tar -zxvf basidiomycota_odb9.tar.gz
```

<br />

## 4. Install GeneMark

Go to this site and download GeneMark-ES/ET.
http://topaz.gatech.edu/GeneMark/license_download.cgi
Don't forget to download the key, too.

### 4.1. Uncompress downloaded files

```
mkdir $FUNGAP_DIR/external/
mv gm_et_linux_64.tar.gz gm_key_64.gz $FUNGAP_DIR/external/  # Move your downloaded files to this directory
cd $FUNGAP_DIR/external/
tar -zxvf gm_et_linux_64.tar.gz
gunzip gm_key_64.gz
cp gm_key_64 ~/.gm_key
```

### 4.2. Install required perl modules for GeneMark

(if required) You may need to install certain Perl modules. Because GeneMark forces to use `/usr/bin/perl` instead of conda-installed perl, you should install the modules for `/usr/bin/perl` (i.e., not in conda environment). Alternatively, you can modify first lines of GeneMark perl scripts from `#!/usr/bin/perl` to `#!/usr/bin/env perl`

```
conda deactivate
sudo apt-get update
sudo apt-get install build-essential
sudo cpan App::cpanminus  # Install cpanm if you do not have one
sudo cpanm Hash::Merge Logger::Simple Parallel::ForkManager YAML
conda activate fungap
```

### 4.3 Check GeneMark and its dependencies are correctly installed.

```
cd $FUNGAP_DIR/external/gm_et_linux_64/gmes_petap
./gmes_petap.pl
```

<br />

### 5. RepeatModeler installation

**Note: RepeatModerler is available in Anaconda2 (https://anaconda.org/bioconda/repeatmodeler), but conda-installed program did not work at the moment. Installation seemed okay, but no result was produced. I will update this whenever working RepeatModeler is available.**

### 5.1. Check perl version.

```
perl -v
```

It should be > 5.8.8

### 5.2. Install RECON 1.08

```
cd $FUNGAP_DIR/external/
wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
tar -zxvf RECON-1.08.tar.gz
cd RECON-1.08/src/
make
make install
```

### 5.3. Install RepeatScout 1.0.5

```
cd $FUNGAP_DIR/external/
wget http://www.repeatmasker.org/RepeatScout-1.0.5.tar.gz
tar -zxvf RepeatScout-1.0.5.tar.gz 
cd RepeatScout-1
make
```

### 5.4. Install NSEG

```
cd $FUNGAP_DIR/external/
mkdir nseg
cd nseg
wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/genwin.c
wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/genwin.h
wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/lnfac.h
wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/makefile
wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/nmerge.c
wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/nseg.c
wget ftp://ftp.ncbi.nih.gov/pub/seg/nseg/runnseg
make
```

### 5.5. Install RepeatMasker 4.0.8

I could not use conda-installed RepeatMasker for RepeatModeler installation. So I manually installed.

```
cd $FUNGAP_DIR/external/
wget http://www.repeatmasker.org/RepeatMasker-open-4-0-8.tar.gz
tar -zxvf RepeatMasker-open-4-0-8.tar.gz
cd RepeatMasker
perl ./configure
```

- Note: `trf` and `rmblastn` are located at `~/anaconda2/envs/fungap/bin`.

### 5.6. Install RepeatModeler 1.0.11

```
cd $FUNGAP_DIR/external/
wget http://www.repeatmasker.org/RepeatModeler/RepeatModeler-open-1.0.11.tar.gz
tar -zxvf RepeatModeler-open-1.0.11.tar.gz
cd RepeatModeler-open-1.0.11/
perl ./configure
```

 - Note: `trf` and `rmblastn` is located at `~/anaconda2/envs/fungap/bin`

### 5.7. Check the installation

```
cd $FUNGAP_DIR/external/RepeatModeler-open-1.0.11/
./BuildDatabase --help
./RepeatModeler --help
```

<br />

## 6. Configure FunGAP

This script allows users to set and test (by --help command) all the dependencies. If this script runs without any issue, you are ready to run FunGAP!

```
cd $FUNGAP_DIR
python set_dependencies.py \
  --pfam_db_dir db/pfam \
  --busco_db_dir db/busco/basidiomycota_odb9/ \
  --genemark_dir external/gm_et_linux_64/gmes_petap/ \
  --repeat_modeler_dir external/RepeatModeler-open-1.0.11
```
