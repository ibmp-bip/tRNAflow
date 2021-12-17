# tRNAflow: Snakemake Workflow for genome annotation of tRNAs

## Quick start

### Installation
Clone the repository:

`git clone https://github.com/ibmp-bip/tRNAflow.git`

`cd tRNAflow/`  

Create the environment:

`conda env create -f trnaflow-env.yml`

Activate the environment:

`conda activate trnaflow`

# install the R package to split tRNA structure
`R --vanilla -e 'devtools::install("workflow/trnascanimport-patch-david", dependencies = TRUE)'`

## Data
### Download genome file(s) of interest and split nuclear genome & organellar genomes in different files (different parameters used in the tRNAscanSE step)
Store it in genomes/ dir.  
*Genome*-nucl.fa  
*Genome*-orga.fa  

### Create a file with link between accession number and chromosome number for each genome to analyze.
Add chromosome size for gff creation and chr_type to import in plantRNA DB. The chromosome type should be chromosome or unplaced (scaffold, contig).  
Look at example *Genome*-alias.txt in genomes/ dir.  

### Create or modify the sample table describing genome information.
Create or modify the sample table (Genome_description.txt) describing genome information:

| Genome                    | Fasta                  | Model         | Alias                       |
|---------------------------|------------------------|---------------|-----------------------------|
| Arabidopsis_thaliana-nucl | genomes/TAIR10_nucl.fa | eukaryotic    | genome/TAIR10_chr-alias.txt |
| Arabidopsis_thaliana-orga | genomes/TAIR10_orga.fa | organellar    | genome/TAIR10_chr-alias.txt |
| genome_name               | genome/file/path       | tRNAscan_model | chromosome_alias/file/path  |

The genome column have to contains unique names !  

The tRNAscanSE model should be :
* eukaryotic  : search for eukaryotic tRNAs
* general     : use general tRNA model (cytoslic tRNAs from all 3 domains included)
* archae      : search for archaeal tRNAs
* bacteria    : search for bacterial tRNAs
* mito_vert   : search for mitochondrial tRNAs in vertebrate
* mito_mamm   : search for mitochondrial tRNAs in mammal
* organellar  : search for other organellar tRNAs

If you need more tRNAscanSE options, you have to modify the tRNAscanSE rule in workflow/Snakefile.

### Configure the yaml file for Snakemake
Change workdir and  genome description file path in config.yaml

## Run snake workflow
snakemake -j 2 -k -s workflow/Snakefile
NB: the -k parameter (*i.e* --keep-going) is important if you have several genomes to process. If one failed, the workflow will keep going.   

The steps of the workflow are:  
- run tRNAscanSE on each genome (with specified parameters)
- extract up and down sequence of the tDNA genes with seqkit
- parse the tRNA structure with a R package
- concat all the data in a csv file

NB : be careful, tRNAscan does not predict chloroplastic tDNA with intron of type II.
The parameters to predict tDNAs are set with very low stringency to filter it by hand.
