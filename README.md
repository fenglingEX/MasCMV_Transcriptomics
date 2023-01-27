# MasCMV_Transcriptomics
A collection of scripts and Snakemake pipelines used in the characterisation of the Mastomys natalaensis cytomegalovirus (MnatCMV) transcriptome.

## Dependencies

### Conda Environment
mnatcmv_transcriptomics.txt contains most of the packages needed to run the Snakemake pipelines and associated python scripts. 

Conda can be installed from: https://docs.conda.io/en/latest/miniconda.html

Once conda has been installed, you can create the environment with:
 
 `conda create --name mnatcmv_transcriptomics --file mnatcmv_transcriptomics.txt`

### ContextMap2

**ContextMap2** is a tool for accurately mapping RNASeq reads. Here it is used to identify poly-A tail signatures in genes. It is not available on any conda channel and must be downloaded separately from:

https://www.bio.ifi.lmu.de//files/Software/ContextMap/manual/ContextMap-manual.html

It is available as a .jar file and requires Java 1.6 or newer to run. Once you have downloaded it, modify the contextmap_command entry in the config.yml to point to the .jar, e.g.:

`contextmap_jar:
   /home/tools/contextmap2.jar`
   
## Usage

The pipeline consists of 2 Snakemake files and accessory python scripts to be run sequentially. A config.yml supplies supporting information. 
The Snakemake files should be run from the project directory using:

`snakemake --snakefile {filename} -j {allocated threads}.`

### Config File

config.yml stores the following:

#### Sample/Animal Information:

- Sample IDs under `all_samples`, starting with animal ID followed by Illumina run identifiers.
- Reference genome names, for mapping reads onto.
- Sample-reference map, denoting which samples ought to be mapped to which reference genomes

#### Mundane processing information:

- Expected file extensions (e.g. `fastq`, `fq.gz`)
- ContextMap2 .jar filepath

### Snakemake Pipeline

 The pipeline consists of 2 separate Snakemake files, to be run sequentially. 

#### Snakemap

Maps stranded RNASeq reads onto reference genomes and then separates them into reads representing sense and antisense transcripts. 
STAR then re-maps these segregated reads onto the reference genome to identify potential splice junctions. A custom script, getIntronCoverage_directional.py, is then run to quantify how many reads span splice junctions.

Run with `snakemake --snakefille Snakemap -j {threads}

#### Snaketail

Maps RNASeq reads onto reference genomes and then identifies potential poly-A sites using a custom script, get_polyA_siganture.py

Run with `snakemake --snakefille Snaketail -j {threads}`
