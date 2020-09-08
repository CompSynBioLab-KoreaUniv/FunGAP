# FunGAP: Fungal Genome Annotation Pipeline v1.0.1

**Last updated: Jan 7, 2019**

**FunGAP is freely available for academic use. For the commerical use or license of FunGAP, please contact In-Geol Choi (igchoi (at) korea.ac.kr). Please, cite the following reference**

Reference: Byoungnam Min,  Igor V Grigoriev, and In-Geol Choi, **FunGAP: Fungal Genome Annotation Pipeline using evidence-based gene model evaluation (2017), Bioinformatics**, Volume 33, Issue 18, Pages 2936–2937, https://doi.org/10.1093/bioinformatics/btx353

<hr>

# Usage of FunGAP

* [0. Prerequisites](#prerequisites)
* [1. Preparing protein database](#protdb)
* [2. Augustus species model](#augustusgenemodel)
* [3. Running FunGAP](#runningfungap)
* [4. FunGAP output](#output)
* [5. Test dataset](#testdata)
* [6. After FunGAP](#afterfungap)

<a name="prerequisites"></a>
## 0. Prerequisites

To run FunGAP, users are required to prepare three main arguments

 - Genome assembly (FASTA)
 - Transcriptomic reads (FASTQ)
 - Protein database (FASTA)

Currently, FunGAP takes only Illumina-sequenced reads (paired-end or single-read). Paired-end read FASTQ files should have formatted file names such as `XXXX_1.fastq` and `XXXX_2.fastq`. An example would be `hyphae_1.fastq` and `hyphae_2.fastq`. For single-read, it should be like `XXXX_s.fastq`. Also BAM file is acceptable with `--trans_bam` option.

<a name="protdb"></a>
## 1. Preparing protein database

FunGAP requires `PROTEIN DATABASE` in FASTA file. We recommend three or four relatives' proteome to reduce computing time. For convenience, we provide a script `download_sister_orgs.py` to build your own protein database for your genome using NCBI API.

Example command (`$FUNGAP_DIR` is your FunGAP installation directory):
```
python $FUNGAP_DIR/download_sister_orgs.py \
  --download_dir sister_orgs \
  --taxon "Schizophyllum" \
  --num_sisters 3 \
  --email_address mbnmbn00@gmail.com
```

E-mail address is needed for NCBI Entrez. All taxon levels are allowed for `--taxon` argument, but *genus* level is appropriate. Now make a protein database.

```
cd sister_orgs/
zcat ./*faa.gz > prot_db.faa
```

You can now input `prot_db.faa` in the `--sister_proteome` argument. 

<a name="augustusgenemodel"></a>
## 2. Augustus species model

Augustus gene predictor requires to select pre-defined species model.
Run `augustus --species=help` to print out the model list. FunGAP provides a script `get_augustus_species.py` to help choose proper species model.

Example command (`$FUNGAP_DIR` is your FunGAP installation directory):
```
python $FUNGAP_DIR/get_augustus_species.py \
  --genus_name "Schizophyllum" \
  --email_address mbnmbn00@gmail.com
```

This will suggest *coprinus_cinereus* and *laccaria_bicolor*.

<a name="runningfungap"></a>
## 3. Running FunGAP

Usage (`$FUNGAP_DIR` is your FunGAP installation directory):
```
python $FUNGAP_DIR/fungap.py \
  --output_dir <output_directory> \
  --trans_read_1 <transcriptome_reads_fastq_1> \
  --trans_read_2 <transcriptome_reads_fastq_2> \
  --genome_assembly <genome_assembly_fasta> \
  --augustus_species <augustus_species> \
  --sister_proteome <sister_proteome> \
  --num_cores <number_of_cpus_to_be_used>
```

<a name="output"></a>
## 4. Output
Final output will be located in `fungap_out` directory

- fungap_out_prot.faa
- fungap_out.gff3
- fungap_out_stats.html

<a name="testdata"></a>
## 5. Test dataset
You can download yeast (*S. cerevisiae*) genome assembly (FASTA) and RNA-seq reads (two FASTQs) from NCBI for testing FunGAP.

```
# Download RNA-seq reads using SRA toolkit (https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)
# Parameter -X indicates that we only need <int> pairs from the dataset.
fastq-dump -X 3000000 -I --split-files SRR1198667

# Download assembly
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz
```

It took about 9 hours by dual Intel(R) Xeon(R) CPU E5-2670 v3 with 40 CPU cores.

<a name="afterfungap"></a>
## 6. Add Pfam annotation to FunGAP GFF3 output
[Interproscan](https://www.ebi.ac.uk/interpro/search/sequence-search) can infer the functions of predicted genes.
The ```gff3_add_pfam.py``` script adds the annotation to the GFF3 file.

Example commands:
```
$ interproscan.sh -i <protein.fasta> -f tsv -appl Pfam --goterms -pa --iprlookup -b <base_name> --tempdir <TEMP-DIR>
$ gff3_add_pfam.py --input_gff3 <fungap_out.gff3> --pfam_file <interproscan_output.tsv>
```
