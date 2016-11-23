******************** BUSCO README File ********************

A. Introduction
B. BUSCO setup
C. Running BUSCO assessments
D. BUSCO options
E. Output files
F. BUSCO setup test with sample data

BUSCO v1.1b May 2015


***********************************************************

A. Introduction
------------------

BUSCO completeness assessment employs sets of Benchmarking Universal Single-Copy
Orthologs from OrthoDB (www.orthodb.org) to provide quantitative measures of the
completeness of genome assemblies, annotated gene sets, and transcriptomes in terms of
expected gene content. Genes that make up the BUSCO sets for each major lineage are selected
from orthologous groups with genes present as single-copy orthologs in at least 90% of the
species. While allowing for rare gene duplications or losses, this establishes an evolutionary
informed expectation that these genes should be found as single-copy orthologs in the genome
of any newly-sequenced species.

Usage of the BUSCO software requires a working installation of Python 3, HMMER 3.1, Blast+,
Augustus (genome assessment only) and EMBOSS transeq (transcriptome assessment only).
BUSCO genome assembly assessment first identifies candidate regions from the genome to be
assessed with tBLASTn searches using BUSCO consensus sequences. Gene structures are then
predicted using Augustus with BUSCO block profiles. Finally, these predicted genes, or all genes
from an annotated gene set or transcriptome, are assessed using HMMER and lineage- specific
BUSCO profiles to classify matc


***********************************************************

B. BUSCO setup
------------------

The BUSCO distribution is released as a compressed archive file (BUSCO_v1.1b.tar.gz) for
download. 

Extracting the files to your current directory "tar -zxvf BUSCO_v1.1b.tar.gz will"
create the directory BUSCO, containing the required files.

Depending on the species you wish to assess, you should now download the appropriate lineage-
specific profile libraries: Metazoa (M), Eukaryota (E), Arthropoda (A), Vertebrata (V), Fungi (F),
or Bacteria (B) from http://busco.ezlab.org to your BUSCO directory.

Before you begin, you will need to make sure that the following required software (some only
required for genome or transcriptome assessments) are installed and accessible from the
command-line, e.g. set environment variable PATH=$PATH:/path/to/software/bin

- Python 3 
Present in most linux repositories

- NCBI BLAST+ 
http://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

- HMMER (HMMER 3.1b2) 
http://hmmer.janelia.org/

- Augustus 3.0.x (genome only)
http://bioinf.uni-greifswald.de/augustus/

**Make sure that the environmental variable AUGUSTUS_CONFIG_PATH was set during installion
(e.g: export AUGUSTUS_CONFIG_PATH=/my_path_to_AUGUSTUS/augustus/config/)**

- EMBOSS tools 6.x.x (transcriptome only)
ftp://emboss.open-bio.org/pub/EMBOSS/

Please check that you have the correct version of the software above.


***********************************************************

C. Running BUSCO assessments 
------------------

Before starting be sure that the necessary software is installed and accessible from the command-line. 
This can be accomplished by adding the necessary paths to your environmental variables (e.g. PATH=$PATH:/path/to/software/bin)

One last thing to do before running BUSCO

Choose and download library of lineage-specific BUSCO data (http://busco.ezlab.org), we recommend using
the largest library possible for your species. Decompress it (e.g. tar -zxvf library.tar.gz) and make 
sure this directory accessible (i.e. move or symlink) from the directory where the analysis is being run. 


1- Genome assembly assessment:

python BUSCO_v1.1b.py -o NAME -in ASSEMBLY -l LINEAGE ?m genome

NAME	name to use for the run and all temporary files
ASSEMBLY	genome assembly file in fasta format
LINEAGE	 path to the lineage to be used (-l /path/to/lineage)

2- Gene set assessment:

python BUSCO_v1.1b.py -o NAME -in GENE_SET -l LINEAGE -m OGS

NAME	name to use for the run and temporary files
GENE_SET	gene set protein sequence file in fasta format
LINEAGE	path to the lineage to be used (-l /path/to/lineage)


3- Transcriptome assessment:

python BUSCO_v1.1b.py -o NAME -in TRANSCRIPTOME -l LINEAGE -m trans

NAME	name to use for the run and temporary files
TRANSCRIPTOME	transcript set sequence file in fasta format
LINEAGE	path to the lineage to be used (-l /path/to/lineage)


***********************************************************

D. BUSCO options:
------------------

1- Mandatory arguments

-o name	Name used for naming output files

-in input_file	Genome assembly / gene set / transcriptome in fasta format

-l lineage	Name of the BUSCO lineage data to be used 
Valid options: metazoa, eukaryota, arthropoda, vertebrata, fungi, and bacteria

-m mode	Mode of analysis 
Valid options: genome, ogs, trans
Default: genome

2- Optional arguments

-h ?help	Print help

-c integer	Number of CPU threads to be used
Default: 1

-sp species	Select from the pre-computed Augustus metaparameters
Selecting a closely-related species usually produces better results
Valid options: see Augustus help for list of options
Default: generic

-e evalue	Use a custom blast e-value cutoff
Default: 0.01

-f	Force overwriting of results files from a previous run with the same name

--flank N	Custom flanking genomic regions in base pairs (bp)
Used when extending selected candidate regions before gene prediction
Default: Automatically calculated flank sizes based on genome size

--long	Performs full optimization for Augustus gene finding training
Default: Off


***********************************************************

E. BUSCO Output
------------------

Successful execution of the BUSCO assessment pipeline will create a directory named name_OUTPUT 
where ?name? is your assigned name for the assessment run. The directory will contain several 
files and directories:

1- Files

short_summary_	Contains summary results in BUSCO notation
and a brief breakdown of the metrics

full_table_	Complete results in tabular format with
coordinates, scores and lengths of BUSCO matches

training_set_	Set of complete BUSCO matches used for training Augustus
Only created during genome assessment

_tblastn	Results in tabular format of tBLASTn searches
with BUSCO consensus sequences

2- Directories

augustus_	Augustus-predicted genes
Only created during genome assessment

augutus_proteins	Corresponding Augustus-predicted proteins
Only created during genome assessment

Selected	Complete BUSCO matches, used for training Augustus

gb	Complete BUSCO matches, GenBank format

gffs	Complete BUSCO matches, GFF format

hmmer_output	Tabular format HMMER output of searches with BUSCO HMMs


***********************************************************

F. BUSCO setup test with sample data
------------------

Sample data are provided to test your BUSCO setup. Execute the following commands 
and compare the final output ?run_SAMPLE? with the provided files in ?run_TEST?.

1.	Change directory to ?sample_data?

cd sample_data/

2.	Run BUSCO assessment on sequence file ?target.fa? in genome mode.

python BUSCO_v1.1b.py -in target.fa -o SAMPLE -l ./example -m genome

3.	Compare the final output ?run_SAMPLE? with the provided files in ?run_TEST?.


***********************************************************

http://busco.ezlab.org