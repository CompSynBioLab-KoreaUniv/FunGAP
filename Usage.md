## Preparing protein database

FunGAP requires `PROTEIN DATABASE` in FASTA file. We recommend three or four relatives' proteome to save running time. For convenience, we provide a script `download_sister_orgs.py` to build your own protein database for your genome using NCBI API.

Assume you downloaded FunGAP in `$HOME/FunGAP`

Usage of `download_sister_orgs.py`
```
python $HOME/FunGAP/download_sister_orgs.py\
 --download_dir <download_directory>\
 --taxon <taxon>\
 --num_sisters <number_of_sisters>\
 --email_address <your_email_address>
```

E-mail address is needed for NCBI Entrez.

Example command.
```
python $HOME/FunGAP/download_sister_orgs.py\
 --download_dir sister_orgs\
 --taxon "Schizophyllum"\
 --num_sisters 3\
 --email_address mbnmbn00@gmail.com
```

All taxon levels are allowed for `--taxon` argument, but `genus` level is appropriate.<br>
Now make protein database.

```
cd sister_orgs/
zcat ./*faa > prot_db.faa
```

You can now input `prot_db.faa` to FunGAP as `PROTEIN DATABASE`. 


## Running FunGAP

To run FunGAP, you need three main arguments

 - Genome assembly (FASTA)
 - Transcriptomic reads (FASTQ)
 - Protein database (FASTA)

Currently, FunGAP gets only **Illumina paired-end** reads files. The file names of two paired-end reads should have identical prefixes like `XX_1.fastq` and `XX_2.fastq`. For example, File names will be like `hyphae_1.fastq` and `hyphae_2.fastq`.

Usage (Assume you downloaded FunGAP in $HOME/FunGAP)
```
python $HOME/FunGAP/fungap.py\
 --output_dir <output_directory>\
 --trans_read_files <transcriptome_reads_fastqs>\
 --project_name <project_name_without_space>\
 --genome_assembly <genome_assembly_fasta>\
 --augustus_species <augustus_species>\
 --org_id <organism_id>\
 --sister_proteome <sister_proteome>\
 --num_cores <number_of_cpus_to_be_used>\
```

- Organism ID will be used in naming gene ID. For example, gene IDs will look like `Hypma_00001.t1`, `Hypma_00002.t1`, and so on.

### Augustus species
You should provide one `augustus_species` used in Augustus.<br>
This is the list what Augustus provides.

Phylum | Class | Species | augustus_species |
--- | --- | --- | ---
Ascomycota | Eurotiomycetes | *Aspergillus fumigatus* | aspergillus_fumigatus |
Ascomycota | Eurotiomycetes | *Aspergillus nidulans* | aspergillus_nidulans |
Ascomycota | Eurotiomycetes | *Aspergillus oryzae* | aspergillus_oryzae |
Ascomycota | Eurotiomycetes | *Aspergillus terreus* | aspergillus_terreus |
Ascomycota | Leotiomycetes | *Botrytis cinerea* | botrytis_cinerea |
Ascomycota | Saccharomycetes | *Candida albicans* | candida_albicans |
Ascomycota | Saccharomycetes | *Candida guilliermondii* | candida_guilliermondii |
Ascomycota | Saccharomycetes | *Candida tropicalis* | candida_tropicalis |
Ascomycota | Sordariomycetes | *Chaetomium globosum* | chaetomium_globosum |
Ascomycota | Eurotiomycetes | *Coccidioides immitis* | coccidioides_immitis |
Basidiomycota | Agaricomycetes | *Coprinus cinereus* | coprinus |
Basidiomycota | Agaricomycetes | *Coprinus cinereus* | coprinus_cinereus |
Basidiomycota | Agaricomycetes | *Cryptococcus neoformans gattii* | cryptococcus_neoformans_gattii |
Basidiomycota | Agaricomycetes | *Cryptococcus neoformans gattii* | cryptococcus_neoformans_neoformans_B |
Basidiomycota | Agaricomycetes | *Cryptococcus neoformans gattii* | cryptococcus_neoformans_neoformans_JEC21 |
Ascomycota | Saccharomycetes | *Debaryomyces hansenii* | debaryomyces_hansenii |
Microsporidia | | *Encephalitozoon cuniculi* | encephalitozoon_cuniculi_GB |
Ascomycota | Saccharomycetes | *Eremothecium gossypii* | eremothecium_gossypii |
Ascomycota | Sordariomycetes | *Fusarium graminearum* | fusarium_graminearum |
Ascomycota | Eurotiomycetes | *Histoplasma capsulatum* | histoplasma_capsulatum |
Ascomycota | Saccharomycetes | *Kluyveromyces lactis* | kluyveromyces_lactis |
Basidiomycota | Agaricomycetes | *Laccaria bicolor* | laccaria_bicolor |
Ascomycota | Saccharomycetes | *Lodderomyces elongisporus* | lodderomyces_elongisporus |
Ascomycota | Sordariomycetes | *Magnaporthe grisea* | magnaporthe_grisea |
Ascomycota | Sordariomycetes | *Neurospora crassa* | neurospora_crassa |
Basidiomycota | Agaricomycetes | *Phanerochaete chrysosporium* | phanerochaete_chrysosporium |
Ascomycota | Saccharomycetes | *Pichia stipitis* | pichia_stipitis |
Mucoromycotina | Mucorales | *Rhizopus oryzae* | rhizopus_oryzae |
Ascomycota | Saccharomycetes | Saccharomyces cerevisiae | saccharomyces_cerevisiae_S288C |
Ascomycota | Saccharomycetes | Saccharomyces cerevisiae | saccharomyces_cerevisiae_rm11-1a_1 |
Ascomycota | Schizosaccharomycetes | Schizosaccharomyces pombe | schizosaccharomyces_pombe |
Basidiomycota | Ustilaginomycetes | Ustilago maydis | ustilago_maydis |
Ascomycota | Saccharomycetes | Yarrowia lipolytica | yarrowia_lipolytica |

## Output
Final output will be located in output directory you gave in the arguments

- fungap_output_prot.faa
- fungap_output.gff3
- fungap_output_stats.html

## Trouble-shootings
This is very beta version of the software, so please don't hesitate reporting any bug or error you have encountered at mbnmbn00@korea.ac.kr or mbnmbn00@gmail.com.
