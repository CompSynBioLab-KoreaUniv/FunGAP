# FunGAP in Docker Container

This gist has instructions about runnig [FunGAP pipeline](https://github.com/CompSynBioLab-KoreaUniv/FunGAP) from inside a Docker Container.

Requirements:

- Docker
- 16Gb of available disk space
- [GeneMark-ES/ET](http://topaz.gatech.edu/GeneMark/license_download.cgi) release and it's key (`gmes_linux_64_4.tar.gz` and `gm_key_64.gz`)


## Steps

### Build FunGAP docker image

Be sure you have the following files in the working directory:

`Dockerfile  fungap.conf  gmes_linux_64_4.tar.gz  gm_key_64.gz`

> GeneMark is not free for everybody, so you need to register in order to have gm_* files. If was not for that I could have push FunGAP docker image ready for use in DockerHub. The Docker image will have about 13Gb.

```bash
# 1. Clone FunGAP repository
git clone https://github.com/CompSynBioLab-KoreaUniv/FunGAP.git
# 2. Go to docker directory
cd FunGAP/docker
# 3. Download gmes_linux_64_4.tar.gz and gm_key_64.gz and put it in same directory
# 4. Build the image
docker build -t fungap .
```

### Enter Docker image and execute FunGAP pipeline

1. Go to the directory you have your rna-seq reads and genome fasta.
1. Enter into a docker container of fungap:

    ```bash
    docker run -it -w /fungap_workspace --rm -v $(pwd):/fungap_workspace fungap bash
    ```

1. Go to `/fungap_workspace` and use helper script to get Augustus species.

    ```bash
    $FUNGAP_DIR/get_augustus_species.py \
      --genus_name "Saccharomyces" \
      --email_address byoungnammin@lbl.gov
    ```

1. Make protein database

    ```bash
    $FUNGAP_DIR/download_sister_orgs.py \
      --taxon "Saccharomyces cerevisiae" \
      --email_address <YOUR_EMAIL_ADDRESS> \
      --num_sisters 1
    zcat sister_orgs/*faa.gz > prot_db.faa
    ```

1. Run FunGAP

    ```bash
    $FUNGAP_DIR/fungap.py \
      --genome_assembly GCF_000146045.2_R64_genomic.fna \
      --trans_read_1 SRR1198667_1.fastq \
      --trans_read_2 SRR1198667_2.fastq \
      --augustus_species saccharomyces_cerevisiae_S288C \
      --busco_dataset ascomycota_odb10 \
      --sister_proteome prot_db.faa \
      --num_cores 8
    ```

Now you can exit docker container. Your current working directory was mounted inside FunGAP container (on /fungap_workspace) so all output files will be available on your system.


# FunGAP in Singularity Container

Singularity is another container platform that, in contrast to Docker, allows you to run data analysis pipelines on HPC clusters. The instructions for running FunGAP pipeline from inside a Singularity container are very similar to the ones described above for Docker container.

Requirements to build the container (Singularity image file, or `.sif`):

- [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps)
- Singularity definition file (`fungap.def`)
- Root priviledges
- ~ 16GB of available disk space
- [GeneMark-ES/ET](http://topaz.gatech.edu/GeneMark/license_download.cgi) release and its key (`gmes_linux_64_4.tar.gz` and `gm_key_64.gz`)

In order to build the container from the def file, you need a machine in which you have root priviledges:

## Steps (tested on Ubuntu 22.04.1 LTS)

### Step 1: Install Singularity

1. Install required libraries

```bash
sudo apt-get update && sudo apt-get install -y \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git \
    cryptsetup \
    libglib2.0-dev
```

1. Install GO

```bash
cd ${HOME}  # or wherever you want to install
wget https://go.dev/dl/go1.19.4.linux-amd64.tar.gz
tar -zxvf go1.19.4.linux-amd64.tar.gz
echo "export PATH=$PATH:${HOME}/go/bin" >> ~/.bashrc
source ~/.bashrc
```

2. Install Singularity

```bash
cd ${HOME}  # or wherever you want to install
wget https://github.com/sylabs/singularity/releases/download/v3.10.4/singularity-ce-3.10.4.tar.gz
tar -zxvf singularity-ce-3.10.4.tar.gz
cd singularity-ce-3.10.4
./mconfig  # You can just save and exit without modifying anything
make -C builddir
sudo make -C builddir install
```

### Step 2: Build the FunGAP Singularity image

```bash
# 1. Clone FunGAP repository
git clone https://github.com/CompSynBioLab-KoreaUniv/FunGAP.git
# 2. Go to docker directory
cd FunGAP/docker
# 3. Download gmes_linux_64_4.tar.gz and gm_key_64.gz and put it in same directory
# 4. Build the Singularity image
sudo singularity build --sandbox sandbox fungap.def
zcat gm_key_64.gz > ~/.gm_key
```

This step will create the image directory `sandbox`, which is about ~ 16GB.

If you don't have the necessary disk space nor have root priviledges, another option is to use the [Remote Builder](https://cloud.sylabs.io/builder).

### Step 3: Test the Singularity image

1. Download a test dataset (_Saccharomyces cerevisiae_ S288C)

```bash
# Download transcriptome FASTQ reads
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.3/sratoolkit.2.11.3-ubuntu64.tar.gz
tar -zxvf sratoolkit.2.11.3-ubuntu64.tar.gz
sratoolkit.2.11.3-ubuntu64/bin/vdb-config --interactive
sratoolkit.2.11.3-ubuntu64/bin/fastq-dump -X 1000000 -I --split-files SRR1198667

# Download genome assembly
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
gunzip GCF_000146045.2_R64_genomic.fna.gz
```

2. Download protein database

```bash
export FUNGAP_DIR=/home/ubuntu/FunGAP  # Where you downloaded FunGAP
export GENUS_NAME=Schizophyllum
export EMAIL_ADDRESS=your_email_address@abc.com

${FUNGAP_DIR}/download_sister_orgs.py \
  --download_dir sister_orgs \
  --taxon ${GENUS_NAME} \
  --num_sisters 3 \
  --email_address ${EMAIL_ADDRESS}
zcat sister_orgs/*faa.gz > prot_db.faa
```

3. Run FunGAP

```
screen  # I recommend it because FunGAP can take long time
singularity exec --writable ${FUNGAP_DIR}/docker/sandbox /workspace/FunGAP/fungap.py \
  --output_dir fungap_out \
  --trans_read_1 SRR1198667_1.fastq \
  --trans_read_2 SRR1198667_2.fastq \
  --genome_assembly GCF_000146045.2_R64_genomic.fna  \
  --augustus_species saccharomyces_cerevisiae_S288C  \
  --sister_proteome prot_db.faa  \
  --num_cores 8 \
  --busco_dataset saccharomycetes_odb10
# You can escape the screen by "ctrl + a + d"
# You can check the list of screens by a command "screen -ls"
# You can re-enter a certain screen by a command "screen -r <screen_id>"
```
