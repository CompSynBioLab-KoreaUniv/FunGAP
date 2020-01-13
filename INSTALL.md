# Installation of FunGAP v1.0.1

**Last updated: Jan 13, 2020*

**FunGAP is freely available for academic use. For the commerical use or license of FunGAP, please contact In-Geol Choi (email: igchoi (at) korea.ac.kr). Please, cite the following reference**

Reference: Byoungnam Min  Igor V Grigoriev  In-Geol Choi, FunGAP: Fungal Genome Annotation Pipeline using evidence-based gene model evaluation (2017), Bioinformatics, Volume 33, Issue 18, Pages 2936–2937, https://doi.org/10.1093/bioinformatics/btx353

<hr>
Because FunGAP implements many dependent programs, you may encounter issues during
installation. Please don't hesitate to post on *Issues* or contact me (mbnmbn00@gmail.com) for help.

These steps were tested and confirmed in freshly installed Ubuntu 18.04 LTS.

# Insall FunGAP using Docker

Using Docker is the most reliable and robust way to install FunGAP.
https://gist.github.com/lmtani/d37343a40e143b59336e4606055d1723

* This Docker container has been developed and maintained by [Lucas Taniguti](https://gist.github.com/lmtani)

# Install FunGAP using conda

Although we recommend using Docker, some workspaces are not available for Docker (e.g., HPC). Please use the following instruction for conda-based FunGAP installation.

## 0. FunGAP requirements

### 0.1. Required softwares (and tested versions)

1. [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) v2.1.0
1. [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) v2.9.0
1. [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) v1.0.11
1. [Maker](http://www.yandell-lab.org/software/maker.html) v2.31.10
1. [GeneMark-ES/ET](http://topaz.gatech.edu/GeneMark/license_download.cgi) v4.48_3.60_lic
1. [Augustus](https://github.com/Gaius-Augustus/Augustus) v3.3
1. [Braker](http://exon.gatech.edu/braker1.html) v1.9
1. [BUSCO](https://busco.ezlab.org/) v3.0.2
1. [Pfam_scan](https://www.ebi.ac.uk/seqdb/confluence/display/THD/PfamScan) v1.6-2
1. [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) v2.6.0+
1. [Samtools](http://www.htslib.org/download/) v1.9
1. [Bamtools](https://github.com/pezmaster31/bamtools) v2.4.1

### 0.2. Required databases

1. [BUSCO](https://busco.ezlab.org/) odb9
1. [Pfam](https://pfam.xfam.org/) release 32.0

<br/>

## 1. Setup Anaconda environment

For robust installation, we recommend to use Anaconda environment and install dependent programs and libraries as much as possible in the environment.

### 1.1. Install Anaconda2 (v4.8.1 tested)

Download and install Anaconda2 (We assume that you install it in ```$HOME/anaconda2```)

```
cd $HOME
wget https://repo.anaconda.com/archive/Anaconda2-2019.10-Linux-x86_64.sh
bash Anaconda2-2019.10-Linux-x86_64.sh
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

This step is **essential**; otherwise, Maker will stop.

```
conda config --remove channels bioconda
conda config --remove channels conda-forge
conda config --add channels bioconda/label/cf201901
conda config --add channels conda-forge/label/cf201901
```

### 1.5. Install dependencies

```
conda install augustus rmblast maker hisat2 braker busco=3.0.2 blast pfam_scan bowtie2
conda install -c bioconda/label/cf201901 jellyfish  # For Trinity
conda install -c anaconda openjdk  # For Trinity
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
cd FunGAP/
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
wget https://busco-archive.ezlab.org/v3/datasets/fungi_odb9.tar.gz
wget https://busco-archive.ezlab.org/v3/datasets/ascomycota_odb9.tar.gz
wget https://busco-archive.ezlab.org/v3/datasets/basidiomycota_odb9.tar.gz
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

### 4.2. Change perl path

GeneMark forces to use `/usr/bin/perl` instead of conda-installed perl. You can change this by running `change_path_in_perl_scripts.pl` script.

```
cd $FUNGAP_DIR/external/gm_et_linux_64/
perl change_path_in_perl_scripts.pl "/usr/bin/env perl"
```

### 4.3 Check GeneMark and its dependencies are correctly installed.

```
cd $FUNGAP_DIR/external/gm_et_linux_64/
./gmes_petap.pl
```

<br />

## 5. RepeatModeler installation

Note: RepeatModerler is available in Anaconda2 (https://anaconda.org/bioconda/repeatmodeler), but the conda-installed program does not work at the moment. Installation seemed okay, but when I ran, I got no results. I will update this whenever working RepeatModeler is available.

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
sudo apt-get install build-essential  # "make: cc: Command not found" error
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

## 6. Trinity installation

Download and compile Trinity

```
cd $FUNGAP_DIR/external
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.9.0/trinityrnaseq-v2.9.0.FULL.tar.gz
tar -zxvf trinityrnaseq-v2.9.0.FULL.tar.gz
cd trinityrnaseq-v2.9.0/
conda deactivate  # Compile outside conda environment
sudo apt-get install cmake zlib1g-dev  # "zlib.h: No such file or directory" error
make
make plugins
```

Add to `$PATH` variable
```
echo "export PATH=$PATH:$FUNGAP_DIR/external/trinityrnaseq-Trinity-v2.8.5/" >> ~/.bashrc
source ~/.bashrc
```

### 6-1. Salmon installation

```
cd $FUNGAP_DIR/external
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.1.0/salmon-1.1.0_linux_x86_64.tar.gz
tar -zxvf salmon-1.1.0_linux_x86_64.tar.gz 
echo "export PATH=$PATH:$FUNGAP_DIR/external/salmon-latest_linux_x86_64/bin/" >> ~/.bashrc
source ~/.bashrc
```

<br />

## 7. Configure FunGAP

This script allows users to set and test (by --help command) all the dependencies. If this script runs without any issue, you are ready to run FunGAP!

```
cd $FUNGAP_DIR
conda activate fungap
python set_dependencies.py \
  --pfam_db_dir db/pfam \
  --busco_db_dir db/busco/basidiomycota_odb9/ \
  --genemark_dir external/gm_et_linux_64/ \
  --repeat_modeler_dir external/RepeatModeler-open-1.0.11
```

<br />

## 8. Braker1 bug

You have to fix this bug; otherwise, you will encounter this error.

> ERROR: Number of good genes is 0, so the parameters cannot be optimized. Recomended are at least 300 genes <br />
> WARNING: Number of good genes is low (0 <br />
> ). Recomended are at least 300 genes

```
cd $HOME/anaconda2/envs/fungap/bin
vim filterGenesIn_mRNAname.pl
```

Go to line 27, and add "?" character.

From
```
if($_ =~ m/transcript_id \"(.*)\"/) {
```
to
```
if($_ =~ m/transcript_id \"(.*?)\"/) {
```
