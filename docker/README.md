# FunGAP in Docker Container

This gist has instructions about runnig [FunGAP pipeline](https://github.com/CompSynBioLab-KoreaUniv/FunGAP) from inside a Docker Container.

Requirements:

- Docker
- 16Gb of available disk space
- [GeneMark-ES/ET](http://topaz.gatech.edu/GeneMark/license_download.cgi) release and it's key (`gmes_linux_64.tar.gz` and `gm_key_64.gz`)


## Steps

### Build FunGAP docker image

Be sure you have the following files in the working directory:

`Dockerfile  fungap.conf  gmes_linux_64.tar.gz  gm_key_64.gz  patch_braker.pl`

> GeneMark is not free for everybody, so you need to register in order to have gm_* files. If was not for that I could have push FunGAP docker image ready for use in DockerHub. The Dcoker image will have about 13Gb.

```bash
# 1. Clone FunGAP repository
git clone https://github.com/CompSynBioLab-KoreaUniv/FunGAP.git
# 2. Go to docker directory
cd FunGAP/docker
# 3. Download gmes_linux_64.tar.gz and gm_key_64.gz and put it in same directory
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
    python /workspace/FunGAP/get_augustus_species.py \
      --genus_name "Saccharomyces" \
      --email_address byoungnammin@lbl.gov
    ```

1. Make protein database

    ```bash
    python /workspace/FunGAP/download_sister_orgs.py \
      --taxon "Saccharomyces" \
      --email_address byoungnammin@lbl.gov
    zcat sister_orgs/*faa.gz > prot_db.faa
    ```

1. Run FunGAP

    ```bash
    python /workspace/FunGAP/fungap.py \
      --output_dir fungap_out \
      --trans_read_1 SRR1198667_1.fastq \
      --trans_read_2 SRR1198667_2.fastq \
      --genome_assembly GCF_000146045.2_R64_genomic.fna  \
      --augustus_species saccharomyces_cerevisiae_S288C  \
      --sister_proteome prot_db.faa  \
      --num_cores 8
    ```

Now you can exit docker container. Your current working directory was mounted inside FunGAP container (on /fungap_workspace) so all output files will be available on your system.
