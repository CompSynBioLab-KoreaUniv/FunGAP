# FunGAP: Fungal Genome Annotation Pipeline

FunGAP performs gene prediction on given genome assembly and RNA-seq reads. See **INSTALL.md** and **USAGE.md** for installation and usage, or you can go wiki tab for the same.

* [Download FunGAP](#download)
* [BLAST+ installation](#blast)
* [Trinity installation](#trinity)
 * [Maker2 installation](#maker)

## FunGAP INPUT & OUTPUT

FunGAP inputs:
```
--output_dir                      | Output directory
--trans_read_files                | Illumina paired-end mRNA reads files (FASTQ)
--project_name                    | Project name without white space
--genome_assembly                 | Genome assembly file (FASTA)
--augustus_species                | Augustus --species argument
--org_id                          | Organism ID for assigning gene IDs
--sister_proteome                 | Protein database (FASTA)
--num_cores                       | Number of CPU cores to be used
```
FunGAP outputs:
```
fungap_output.gff3                | Tab-delimited genomic feature format
fungap_output_prot.faa            | Translated protein sequences
fungap_output_result_summary.html | Summary of FunGAP results
```

## Pipeline description



## Contact

* Project principal investigator: In-Geol Choi at Korea University
* Contact: mbnmbn00@korea.ac.kr (or mbnmbn00@gmail.com)

If you have any problem to install or run, please don't hesistate to contact me. I will help you as much as I can and that helps build more robust program.



![](http://compbio.korea.ac.kr/bnmin/fungap/fungap_logo-01.png)

