# Installation of FunGAP v1.1.0

Because FunGAP implements many dependent programs, you may encounter issues during
installation. Please don't hesitate to contact me (mbnmbn00@@gmail.com) for help.

These steps were tested and confirmed in freshly installed Ubuntu 18.04 LTS.

<a name="download"></a>
## 0. FunGAP requirements and tested versions

### 0.1 Required softwares

1. Hisat2 v2.1.0
1. Trinity v2.6.6
1. RepeatModeler v1.0.11
1. Maker v
1. BUSCO v3.0.2
1. Pfam_scan v1.6
1. BLAST v2.6.0+
1. Samtools v1.9
1. Bamtools v2.4.1


Download FunGAP using GitHub clone. Suppose we are installing FunGAP in your `$HOME` directory, but you are free to change the location.

```
cd $HOME
git clone https://github.com/CompSynBioLab-KoreaUniv/FunGAP.git
```

<a name="blast"></a>
## BLAST+ installatioon
**BLAST+** is used in Maker and BUSCO running.

```
sudo apt-get install ncbi-blast+
```

<a name="trinity"></a>
## Trinity installation
**Trinity** performs efficient and robust *de novo* reconstruction of transcriptomes from RNA sequencing data. https://github.com/trinityrnaseq/trinityrnaseq/wiki

Download and install Trinity v2.2.0 using github.

```
cd $HOME/FunGAP/external
git clone https://github.com/trinityrnaseq/trinityrnaseq.git
cd trinityrnaseq
make
```

<a name="maker"></a>
## Maker2 installation
**Maker2** is an easy-to-use annotation pipeline designed for emerging model organism genomes. http://www.gmod.org/wiki/MAKER

**\#\#\# Please note that you need a proper license to use Maker2. \#\#\#**

Download and move `maker-2.31.8.tgz` to your FunGAP directory
```
mv maker-2.31.8.tgz $HOME/FunGAP/external/
```

Unzip maker2
```
cd $HOME/FunGAP/external/
tar -zxvf maker-2.31.8.tgz
```

Install Maker2 pre-requisites
```
cd $HOME/FunGAP/external/maker/src
sudo apt-get install libpq-dev
perl Build.PL
sudo ./Build installdeps
./Build installexes
./Build install
```

Configure RepeatMasker. First, download repbase manually from http://www.girinst.org/server/RepBase/index.php
<br> Then move it to `$HOME/FunGAP/external/maker/exe/RepeatMasker/`

```
cd $HOME/FunGAP/external/maker/exe/RepeatMasker/
tar -zxvf repeatmaskerlibraries-20150807.tar.gz
./configure

# **TRF PROGRAM**
# This is the full path to the TRF program.
# This is now used by RepeatMasker to mask simple repeats.
# Enter path [  ]:
path/to/FunGAP/external/maker/exe/RepeatMasker/trf

# Add a Search Engine:
# 1. CrossMatch: [ Un-configured ]
# 2. RMBlast - NCBI Blast with RepeatMasker extensions: [ Un-configured ]
# 3. WUBlast/ABBlast (required by DupMasker): [ Un-configured ]
# 4. HMMER3.1 & DFAM: [ Un-configured ]

# 5. Done
# Enter Selection:
2

# **RMBlast (rmblastn) INSTALLATION PATH**
# This is the path to the location where
# the rmblastn and makeblastdb programs can be found.
# Enter path [  ]:
path/to/FunGAP/external/maker/exe/RepeatMasker/rmblast/bin
```

<a name="repeatmodeler"></a>
## RepeatModeler installation
**RepeatModeler** is a *de novo* repeat family identification and modeling package.
http://www.repeatmasker.org/RepeatModeler.html

Install RepeatModeler and its dependencies.

Check perl version (ensure version >5.8.8)
```
perl -v
```

Install RepeatModeler
```
cd $HOME/FunGAP/external/
wget http://www.repeatmasker.org/RepeatModeler-open-1-0-8.tar.gz
tar -zxvf RepeatModeler-open-1-0-8.tar.gz
cd RepeatModeler/
perl ./configure

# **REPEATMASKER INSTALLATION PATH**
# This is the path to the location where
# the RepeatMasker program suite can be found.
# Enter path [  ]:
path/to/FunGAP/external/maker/exe/RepeatMasker/

# **RECON INSTALLATION PATH**
# This is the path to the location where
# the RECON program suite can be found.
# Enter path [  ]:
path/to/FunGAP/external/RECON-1.08/bin

# **RepeatScout INSTALLATION PATH**
# This is the path to the location where
# the RepeatScout program suite can be found.
# Enter path [  ]:
path/to/FunGAP/external/RepeatScout-1/

# **TRF INSTALLATION PATH**
# This is the path to the location where
# the TRF program can be found.
# Enter path [  ]:
path/to/FunGAP/external/maker/exe/RepeatMasker

# Add a Search Engine:
# 1. RMBlast - NCBI Blast with RepeatMasker extensions: [ Un-configured ]
# 2. WUBlast/ABBlast: [ Un-configured ]

# 3. Done
# Enter Selection:
1

# **RMBlast (rmblastn) INSTALLATION PATH**
# This is the path to the location where
# the rmblastn and makeblastdb programs can be found.
# Enter path [  ]:
path/to/FunGAP/external/maker/exe/RepeatMasker/rmblast/bin
```

<a name="braker"></a>
## Braker1 installation
Braker1 is an unsupervised RNA-Seq-based genome annotation with GeneMark-ET and AUGUSTUS.
http://exon.gatech.edu/genemark/braker1.html

Install Braker1 and its dependencies.

Copy gm_key to `$HOME`
```
cp $HOME/FunGAP/external/gm_et_linux_64/gm_key ~/.gm_key
```

Install perl modules
```
sudo cpan YAML
sudo cpan App::cpanminus
sudo cpanm File::Spec::Functions
sudo cpanm Hash::Merge
sudo cpanm List::Util
sudo cpanm Logger::Simple
sudo cpanm Module::Load::Conditional
sudo cpanm Parallel::ForkManager
sudo cpanm POSIX
sudo cpanm Scalar::Util::Numeric
sudo cpanm YAML
```

For bamtools,
```
sudo apt-get install zlib1g-dev
```

<a name="busco"></a>
## BUSCO Installation
**Busco**: Benchmarking Universal Single-Copy Orthologs
http://busco.ezlab.org/

Download BUSCO dataset.
```
cd $HOME/FunGAP/data/
wget http://busco.ezlab.org/v1/files/fungi_buscos.tar.gz
tar -zxvf fungi_buscos.tar.gz
```

<a name="interproscan"></a>
## InterProScan installation
**InterProScan** scans a sequence for matches against the InterPro protein signature databases.
https://github.com/ebi-pf-team/interproscan/wiki

Install InterProScan.
```
cd $HOME/FunGAP/external/
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.18-57.0/interproscan-5.18-57.0-64-bit.tar.gz
tar -zxvf interproscan-5.18-57.0-64-bit.tar.gz
```

<a name="pythonmodules"></a>
## Install Python modules
FunGAP requires several python modules and they can be installed by pip.

Install pip
```
sudo apt-get install python-pip
```

Install needed modules
```
sudo pip install biopython  # version 1.65 tested
sudo pip install numpy  # version 1.6.1 tested
sudo pip install networkx  # version 1.1 tested
sudo pip install matplotlib  # version 1.5.x tested
sudo pip install markdown2  # version 2.3.3 tested
```

You can check if FunGAP is correctly installed.
```
python $HOME/FunGAP/check_dependencies.py -o tmp
```
