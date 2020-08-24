# Installation of FunGAP v1.1.0

**Last updated: Aug 20, 2020*

**FunGAP is freely available for academic use. For the commerical use or license of FunGAP, please contact In-Geol Choi (email: igchoi (at) korea.ac.kr). Please, cite the following reference**

Reference: Byoungnam Min  Igor V Grigoriev  In-Geol Choi, FunGAP: Fungal Genome Annotation Pipeline using evidence-based gene model evaluation (2017), Bioinformatics, Volume 33, Issue 18, Pages 2936–2937, https://doi.org/10.1093/bioinformatics/btx353

<hr>

Please don't hesitate to post on *Issues* or contact me (mbnmbn00@gmail.com) for help.
These steps were tested in the freshly installed Ubuntu 18.04 LTS.

<br />

# Install FunGAP using Docker

Using Docker is the most reliable and robust way to install FunGAP. [Please follow the instruction](docker/README.md).

<br />

# Install FunGAP using conda

You may not have superuser privilege (e.g., HPC) required for Docker. To install dependencies, we recommend using the Anaconda.

## 0. FunGAP requirements

### 0.1. Required softwares (and tested versions)

1. [Hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) v2.2.0
1. [Trinity](https://github.com/trinityrnaseq/trinityrnaseq) v2.11.0
1. [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) v2.0.1
1. [Maker](http://www.yandell-lab.org/software/maker.html) v2.31.10
1. [GeneMark-ES/ET](http://topaz.gatech.edu/GeneMark/license_download.cgi) v4.59_lic
1. [Augustus](https://github.com/Gaius-Augustus/Augustus) v3.3.3
1. [Braker](http://exon.gatech.edu/braker1.html) v2.1.5
1. [BUSCO](https://busco.ezlab.org/) v4.1.2
1. [Pfam_scan](https://www.ebi.ac.uk/seqdb/confluence/display/THD/PfamScan) v1.6
1. [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) v2.9.0+
1. [Samtools](http://www.htslib.org/download/) v1.10
1. [Bamtools](https://github.com/pezmaster31/bamtools) v2.5.1

### 0.2. Required database

1. [Pfam](https://pfam.xfam.org/) release 33.1

<br/>

## 1. Setup Anaconda environment

### 1.1. Install Anaconda3 (v4.8.3 tested)

Download and install Anaconda3 (We assume that you install it in `$HOME/anaconda3`)

```
cd $HOME
wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
bash Anaconda3-2020.07-Linux-x86_64.sh
```

### 1.2. Set conda environment

```
echo ". $HOME/anaconda3/etc/profile.d/conda.sh" >> ~/.bashrc
source $HOME/.bashrc
which conda  # It should be $HOME/anaconda3/condabin/conda
```

### 1.3. Add channels

Set up the channels.

```
# Add two channels
conda config --add channels bioconda
conda config --add channels conda-forge

# Check the channels
conda config --show channels
# channels:
#  - conda-forge
#  - bioconda
#  - defaults
  
# Remove channels if you have unnecessary channels
conda config --remove channels bioconda/label/cf201901
conda config --remove channels conda-forge/label/cf201901
```

### 1.4. Create and activate an environment

```
conda update conda
conda create -n fungap
conda activate fungap
```

### 1.5. Install dependencies

```
conda install braker2=2.1.5 trinity=2.11.0 repeatmodeler=2.0.1 hisat2=2.2.0 pfam_scan=1.6 busco=4.1.2
pip install biopython bcbio-gff networkx markdown2 matplotlib
cpanm Hash::Merge Logger::Simple Parallel::ForkManager YAML
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
gunzip Pfam-A.hmm.gz
gunzip Pfam-A.hmm.dat.gz
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
  --maker_path $MAKER_DIR
```

<br />

## 7. Braker bug

You have to fix this bug; otherwise, you will encounter this error.

> ERROR: Number of good genes is 0, so the parameters cannot be optimized. Recomended are at least 300 genes <br />
> WARNING: Number of good genes is low (0 <br />
> ). Recomended are at least 300 genes

```
conda activate fungap
cd $(dirname $(which braker.pl))
vim filterGenesIn_mRNAname.pl
```

Go to line 38, and add a "?" character.

From
```
    if ( $_ =~ m/transcript_id \"(.*)\"/ ) {
```
to
```
    if ( $_ =~ m/transcript_id \"(.*?)\"/ ) {
```

## 8. Test run

<a name="testdata"></a>

### 8.1. Download test dataset

You can download yeast (*Saccharomyces cerevisiae*) genome assembly (FASTA) and RNA-seq reads (two FASTQs) from NCBI for testing FunGAP.

```
# Download RNA-seq reads using SRA toolkit (https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)
# Parameter -X indicates that we only need <int> pairs from the dataset.
fastq-dump -X 3000000 -I --split-files SRR1198667

# Download assembly
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
gunzip GCF_000146045.2_R64_genomic.fna.gz
```

### 8.2. Download protein sequences of related species

```
$FUNGAP_DIR/download_sister_orgs.py \
  --taxon "Saccharomyces cerevisiae" \
  --email_address <your_email_address>
zcat sister_orgs/*faa.gz > prot_db.faa
```

### 8.3. Get Augustus species

```
$FUNGAP_DIR/get_augustus_species.py \
  --genus_name "Saccharomyces" \
  --email_address byoungnammin@lbl.gov
```

 - saccharomyces_cerevisiae_S288C
 
### 8.4. Run FunGAP

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
  
It took about 9 hours by dual Intel(R) Xeon(R) CPU E5-2670 v3 with 40 CPU cores.
