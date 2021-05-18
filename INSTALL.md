# Installation of FunGAP v1.1.1

**Last updated: May 18, 2021*

**FunGAP is freely available for academic use. For the commerical use or license of FunGAP, please contact In-Geol Choi (email: igchoi (at) korea.ac.kr). Please, cite the following reference**

Reference: Byoungnam Min  Igor V Grigoriev  In-Geol Choi, FunGAP: Fungal Genome Annotation Pipeline using evidence-based gene model evaluation (2017), Bioinformatics, Volume 33, Issue 18, Pages 2936–2937, https://doi.org/10.1093/bioinformatics/btx353

<hr>

Please don't hesitate to post on *Issues* or contact me (mbnmbn00@gmail.com) for help.
These steps were tested in the freshly installed Ubuntu 20.04.2 LTS.

<br />

# Install FunGAP using Docker

Using Docker is the most reliable and robust way to install FunGAP. [Please follow the instruction](docker/README.md).

<br />

# Install FunGAP using conda

Although we recommend using Docker, some workspaces are not available for Docker (e.g., HPC). Please use the following instruction for conda-based FunGAP installation.

## 0. FunGAP requirements

### 0.1. Required softwares (and tested versions)

1. [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) v2.2.1
1. [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) v2.12.0
1. [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) v2.0.1
1. [Maker](http://www.yandell-lab.org/software/maker.html) v3.01.03
1. [GeneMark-ES/ET](http://topaz.gatech.edu/GeneMark/license_download.cgi) v4.65_lic
1. [Augustus](https://github.com/Gaius-Augustus/Augustus) v3.4.0
1. [Braker](http://exon.gatech.edu/braker1.html) v2.1.5
1. [BUSCO](https://busco.ezlab.org/) v5.1.2
1. [Pfam_scan](https://www.ebi.ac.uk/seqdb/confluence/display/THD/PfamScan) v1.6
1. [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) v2.11.0
1. [Samtools](http://www.htslib.org/download/) v1.10
1. [Bamtools](https://github.com/pezmaster31/bamtools) v2.5.1

### 0.2. Required database

1. [Pfam](https://pfam.xfam.org/) release 34.0

<br/>

## 1. Setup Anaconda environment

### 1.1. Install Anaconda3 (v4.10.1 tested)

Download and install Anaconda3 (We assume that you install it in `$HOME/anaconda3`)

```
# Download and install conda
cd $HOME
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh
bash Anaconda3-2021.05-Linux-x86_64.sh

# Set environment if you select "no" to "Do you wish the installer to initialize Anaconda3?"
echo ". $HOME/anaconda3/etc/profile.d/conda.sh" >> ~/.bashrc
source $HOME/.bashrc
which conda  # It should be $HOME/anaconda3/condabin/conda

# Get up-to-date conda
conda update conda
```

### 1.2. Install dependencies

```
# Install Mamba package manager (faster!)
conda install mamba -n base -c conda-forge

# Create FunGAP environment and install dependencies using Mamba
conda create -y -n fungap
conda activate fungap
mamba install \
  braker2=2.1.5 trinity=2.12.0 repeatmodeler=2.0.1 hisat2=2.2.1 pfam_scan=1.6 busco=5.1.2 \
  -c bioconda -c conda-forge

# Install Maker using Mamba (Maker installation is conflict with Busco)
conda deactivate
conda create -y -n maker
conda activate maker
mamba install maker=3.01.03 -c bioconda -c conda-forge

# Install Python and Perl modules
pip install biopython bcbio-gff markdown2 matplotlib
cpanm YAML Hash::Merge Logger::Simple Parallel::ForkManager MCE::Mutex Thread::Queue threads
```

### 1.6. Install Maker

Because Maker is incompatible with other dependencies (it requires Python2), we will make a new environment and install the Maker in it.

```
conda deactivate
conda create -n maker -c bioconda maker=2.31.10
```

<br />

## 2. Download and install FunGAP

### 2.1. Download FunGAP

Download FunGAP using GitHub clone. Suppose we are installing FunGAP in your `$HOME` directory, but you are free to change the location. `$FUNGAP_DIR` is going to be your FunGAP installation directory.

```
cd $HOME  # or wherever you want
git clone https://github.com/CompSynBioLab-KoreaUniv/FunGAP.git
export FUNGAP_DIR=$(realpath FunGAP/)
# You can put this export command in the your .bashrc file
# so that you don't need to type every time you run the FunGAP
```

<br />

## 3. Download Pfam

Download Pfam databases in your `$FUNGAP_DIR/db` directory.

### 3.1. Pfam DB download 

ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release

```
mkdir -p $FUNGAP_DIR/db/pfam
cd $FUNGAP_DIR/db/pfam
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
gunzip Pfam-A.hmm.gz Pfam-A.hmm.dat.gz
conda activate fungap
hmmpress Pfam-A.hmm  # HMMER package (would be automatically installed in the above Anaconda step)
```

<br />

## 4. Install GeneMark

Go to the below site and download GeneMark-ES/ET.
http://topaz.gatech.edu/GeneMark/license_download.cgi
Don't forget to download the key, too.

### 4.1. Uncompress downloaded files

```
mkdir $FUNGAP_DIR/external/
mv gmes_linux_64.tar.gz gm_key_64.gz $FUNGAP_DIR/external/  # Move your downloaded files to this directory
cd $FUNGAP_DIR/external/
tar -zxvf gmes_linux_64.tar.gz
gunzip gm_key_64.gz
cp gm_key_64 ~/.gm_key
```

### 4.2. Change the perl path

GeneMark forces to use `/usr/bin/perl` instead of conda-installed perl. You can change this by running `change_path_in_perl_scripts.pl` script.

```
cd $FUNGAP_DIR/external/gmes_linux_64/
perl change_path_in_perl_scripts.pl "/usr/bin/env perl"
```

### 4.3 Check GeneMark and its dependencies are correctly installed.

```
cd $FUNGAP_DIR/external/gmes_linux_64/
./gmes_petap.pl
```

<br />

## 5. Download RepeatMasker databases

```
conda activate fungap
cd $(dirname $(which RepeatMasker))/../share/RepeatMasker
# ./configure downloads required databases
echo -e "\n2\n$(dirname $(which rmblastn))\n\n5\n" > tmp && ./configure < tmp

# It should look like this
ls $(dirname $(which RepeatMasker))/../share/RepeatMasker/Libraries
# Artefacts.embl  Dfam.hmm       RepeatAnnotationData.pm  RepeatMasker.lib.nin  RepeatPeps.lib      RepeatPeps.lib.psq
# CONS-Dfam_3.0   README.meta    RepeatMasker.lib         RepeatMasker.lib.nsq  RepeatPeps.lib.phr  RepeatPeps.readme
# Dfam.embl       RMRBMeta.embl  RepeatMasker.lib.nhr     RepeatMaskerLib.embl  RepeatPeps.lib.pin  taxonomy.dat
```

<br />

## 6. Configure FunGAP

This script allows users to set and test (by --help command) all the dependencies. If this script runs without any issue, you are ready to run FunGAP!

```
cd $FUNGAP_DIR
conda activate maker
export MAKER_DIR=$(dirname $(which maker))
echo $MAKER_DIR  # /home/ubuntu/anaconda3/envs/maker/bin
conda activate fungap
./set_dependencies.py \
  --pfam_db_path db/pfam/ \
  --genemark_path external/gmes_linux_64/ \
  --maker_path ${MAKER_DIR}
```

<br />

# Test run

<a name="testdata"></a>

### 1. Download test dataset

You can download yeast (*Saccharomyces cerevisiae*) genome assembly (FASTA) and RNA-seq reads (two FASTQs) from NCBI for testing FunGAP.

```
# Download RNA-seq reads using SRA toolkit (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit)
# Parameter -X indicates the number of read pairs you want to download
fastq-dump -X 1000000 -I --split-files SRR1198667

# Download assembly
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
gunzip GCF_000146045.2_R64_genomic.fna.gz
```

### 2. Download protein sequences of related species

```
$FUNGAP_DIR/download_sister_orgs.py \
  --taxon "Saccharomyces cerevisiae" \
  --email_address <YOUR_EMAIL_ADDRESS> \
  --num_sisters 1
zcat sister_orgs/*faa.gz > prot_db.faa
```

### 3. Get Augustus species

```
$FUNGAP_DIR/get_augustus_species.py \
  --genus_name "Saccharomyces" \
  --email_address byoungnammin@lbl.gov
```

 - saccharomyces_cerevisiae_S288C
 
### 4. Run FunGAP

```
$FUNGAP_DIR/fungap.py \
  --genome_assembly GCF_000146045.2_R64_genomic.fna \
  --trans_read_1 SRR1198667_1.fastq \
  --trans_read_2 SRR1198667_2.fastq \
  --augustus_species saccharomyces_cerevisiae_S288C \
  --busco_dataset ascomycota_odb10 \
  --sister_proteome prot_db.faa \
  --num_cores 8
  ```
  
It took about 8 hours by Intel(R) Xeon(R) CPU E5-2676 v3 @ 2.40GHz with 8 CPU cores.
