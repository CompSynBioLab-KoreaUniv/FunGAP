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
- ~ 5Gb of available disk space
- [GeneMark-ES/ET](http://topaz.gatech.edu/GeneMark/license_download.cgi) release and it's key (`gmes_linux_64_4.tar.gz` and `gm_key_64.gz`)

In order to build the container from the def file, you need a machine in which you have root priviledges:

```bash
# 1. Clone FunGAP repository
git clone https://github.com/CompSynBioLab-KoreaUniv/FunGAP.git
# 2. Go to docker directory
cd FunGAP/docker
# 3. Download gmes_linux_64_4.tar.gz and gm_key_64.gz and put it in same directory
# 4. Build the Singularity image
sudo singularity build fungap.sif fungap.def
```

This step will create the image file `fungap.sif`, which is about ~ 4.4Gb.

If you don't have the necessary disk space nor have root priviledges, another option is to use the [Remote Builder](https://cloud.sylabs.io/builder).

### Enter Singularity image and execute FunGAP pipeline

1. Transfer the `fungap.sif` image file to the directory you have your sequencing data.

1. Enter the Singularity image, and execute the following:

    ```bash
    singularity shell fungap.sif
    ```

1. Go to `/fungap_workspace` and use helper script to get Augustus species.

    ```bash
    python /workspace/FunGAP/get_augustus_species.py \
      --genus_name "Saccharomyces" \
      --email_address your-email-address
    ```

1. Make protein database

    ```bash
    python /workspace/FunGAP/download_sister_orgs.py \
      --taxon "Saccharomyces" \
      --email_address byoungnammin@lbl.gov

    zcat sister_orgs/*faa.gz > prot_db.faa
    ```

1. Run FunGAP!

    ```bash
    python /workspace/FunGAP/fungap.py \
      --output_dir fungap_out \
      --trans_read_1 SRR1198667_1.fastq \
      --trans_read_2 SRR1198667_2.fastq \
      --genome_assembly GCF_000146045.2_R64_genomic.fna  \
      --augustus_species saccharomyces_cerevisiae_S288C  \
      --sister_proteome prot_db.faa  \
      --num_cores 8 \
      --busco_dataset saccharomycetes_odb10
    ```

1. Exit the image:

    ```bash
    exit
    ```
