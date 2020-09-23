# The TaxIt Pipeline
[![Snakemake](https://img.shields.io/badge/snakemake-=3.13.3-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

*TaxIt: An iterative and automated computational pipeline for untargeted strain
level identification of microbial tandem MS spectra*  
Mathias Kuhring, Thilo Muth, Bernhard Y. Renard (in submission)

*Iterative and semi-automatic strain identification of microbial tandem MS
spectra*  
Mathias Kuhring, Bernhard Y. Renard, 2016, Poster at Metaproteomics Symposium
Magdeburg

## Installation

### Download

Download the pipeline from the repository using git and change to the directory:

    git clone --recursive --depth=1 https://gitlab.com/mkuhring/TaxIt.git
    cd TaxIt

### Dependencies

Dependencies are easily installed via conda. If conda is not yet available on
your system install Miniconda as follows (recommended version 4.5.12):

    wget https://repo.continuum.io/miniconda/Miniconda3-4.5.12-Linux-x86_64.sh; chmod +x Miniconda3-4.5.12-Linux-x86_64.sh
    ./Miniconda3-4.5.12-Linux-x86_64.sh

Then, create an environment and install the dependencies using conda:

    conda create --name taxit
    source activate taxit
    conda install -c defaults --override-channels openjdk=8.0 maven=3.5.3 perl=5 r-ggplot2=2.2.1 rpy2=2 bioconda::snakemake=3.13.3 bioconda::msgf_plus=2016.10.26 bioconda::xtandem=15.12.15.2 bioconda::seqkit=0.7.0

After using the pipeline, you may deactivate the conda environment via:

    source deactivate taxit

### Compiling

The pipeline makes use of several java modules which need to be compiled before
execution. Use the Maven build system (available in the active conda
environment) to automatically resolve dependencies and create required jar
files:

    cd modules_java
    mvn clean package

If internet access is restricted, a proxy may be configured as described in the
[Maven manual](https://maven.apache.org/guides/mini/guide-proxies.html).

## Preparation

Besides an input tandem MS spectra file in MGF format, the TaxIt pipeline needs
a comprehensive reference database for initial species identification (NCBI
RefSeq is recommended, UniProt UniRef is experimental) as well as taxonomic
mapping files and NCBI Taxonomy database dump files.

Note: A small example of spectra and databases is included for easy testing as
described in the next section.

Latest NCBI RefSeq proteins can be found here (choose according to expected
kingdom):
* [Bacteria](ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/)
* [Viruses](ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/)

Download all protein files (*\*.protein.faa.gz*), unpack them and merge them to one fasta file.  
E.g. for bactaria run:

    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/\*.protein.faa.gz &&
    gzip -d bacteria.*.protein.faa.gz &&
    cat bacteria.*.protein.faa > refseq_bacteria.fasta &&
    rm bacteria.*.protein.faa

Or for viruses run:

    wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/\*.protein.faa.gz &&
    gzip -d viral.*.protein.faa.gz &&
    cat viral.*.protein.faa > refseq_viral.fasta &&
    rm viral.*.protein.faa

Id mapping files (i.e. protein accession to NCBI taxid) has to be selected
according to the initial database. However, the NCBI Taxonomy mapping file is
also required independently from the initial database, since strain level
identification is based on automatically downloaded NCBI Protein sequences.

* For NCBI RefSeq, the NCBI Taxonomy mapping file:
[prot.accession2taxid](ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz)
* For UniRef, the UniProt ID Mapping file:
[idmapping_selected.tab](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz)

(Unpack with `gzip -d prot.accession2taxid.gz` resp.
`gzip -d idmapping_selected.tab.gz`)

Finally, the latest NCBI Taxonomy database dump files *nodes.dmp* and
*names.dmp* can be obtained [here](ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)
(unpack with `tar xfzv taxdump.tar.gz`).

NCBI Taxonomy mapping and database dump files can be downloaded and extracted
as follows

    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz &&
    gzip -d prot.accession2taxid.gz
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz &&
    tar xfzv taxdump.tar.gz &&
    rm taxdump.tar.gz

## Execution

The TaxIt pipeline is implemented and therefore executed using the Snakemake
workflow management system. Since the pipeline operates on the directory it
has been executed in, proceed as follows:

1. Create a directory for processing your tandem MS spectra (MGF format) with
TaxIt (working directory).
2. Switch to this directory and create a configuration file (as described 
below).
3. Let Snakemake verify configuration and required files by indicating the TaxIt
Snakefile as well as the configuration as follows (adapt path to TaxIt as necessary):

    snakemake -s ~/your/path/to/TaxIt/Snakefile --configfile snakemake.config.json --cores 4 -pkn

4. Then execute the TaxIt pipeline simple by removing the Snakemake dry-run
parameter "-n":

    snakemake -s ~/your/path/to/TaxIt/Snakefile --configfile snakemake.config.json --cores 4 -pk

After executing the pipeline, you will find the identified strain (NCBI taxon
id and scientific name) in *results/identifid.names.txt* in addition to 
plots illustrating strain ratios.

For simple technical testing of the pipeline, as small minimal mockup
set of sample, database and taxonomy files is provided with the code. Switch
to the directory *example* (it will be used as working directory) and execute
the pipeline:

    cd example
    snakemake -s ../Snakefile --configfile snakemake.config.json --cores 4 -pk

## Configuration

The TaxIt pipeline execution via Snakemake is configured with a json file
(see example below).

Following minimum set of parameters is required:

* **sample** - Path to the tandem MS spectra (*.mgf)
* **database_reference** - Path to the initial search database (*.fasta)
* **taxonomy_mapping_ncbi** - Path to the NCBI id mapping file
(prot.accession2taxid)
* **taxonomy_mapping** - Path to id mapping file for the initial search database
* **taxonomy_nodes** - Path to the NCBI Taxonomy nodes file (nodes.dmp)
* **taxonomy_names** - Path to the NCBI Taxonomy names file (names.dmp)

Following parameters are optional:

* **threads** - Maximum number of threads used (default=4)
* **memory** - Maximum amount of memory used in Java subroutines (default="4G")
* **database_contaminants** Path to a protein contaminant database (*.fasta)
(default="*TaxIt/resources/crap.fasta*")
* **database_host** Path to a host database used to filter host spectra
(*.fasta) (default=None)
* **fdr_cutoff** Target-decoy FDR threshold (default=0.01)
* **stacked_legend_max** Maximum number of labeled taxa in results plots
(default=15)
* **https_proxy_url** Proxy URL used for NCBI Protein downloads (default=None)
* **https_proxy_port** Proxy port used for NCBI Protein downloads (default=None)
* **search_engine** Peptide search engine applied, either xtandem or msgfplus
(default="xtandem")

Furthermore, peptide search engines may be additionally configured with
following parameters which are passed either to the according configuration
file (xtandem_input.xml for X!Tandem) or via command-line parameters (MS-GF+).
See the corresponding documentation
([X!Tandem](http://www.thegpm.org/TANDEM/api/index.html)) or command-line help
(MS-GF+) for more valid parameter options.

* **xtandem_default** X!Tandem default configuration file
(default="*TaxIt/resources/xtandem_default_input.xml*")
* **xtandem_fmme** - Fragment mass tolerance (default=0.4)
* **xtandem_fmmeu** - Fragment mass tolerance unit (default="DA")
* **xtandem_pmmep** - Precursor mass tolerance plus (default=100)
* **xtandem_pmmem** - Precursor mass tolerance minus (default=100)
* **xtandem_pmmeu** - Precursor mass tolerance unit (default="ppm")
* **xtandem_mods_fixed** - Fixed modifications, comma separated (default="57@C")
* **xtandem_mods_variable** - Variable modifications, e.g. "16@M", comma
separated (default=None)
* **xtandem_mods_variable_nterm** - Variable N-terminal modifications, e.g.
"+42.0@[", comma separated (default=None)
* **xtandem_add_params** - Additional parameters for X!Tandem as json dictionary
{"param1" : "value", "param2" : "value", ...} which will be added to the
xtandem_input.xml
* **msgfplus_jar** - Path to MS-GF+ jar file, if not installed via Conda
(default=None)
* **msgfplus_mods** - Path to MS-GF+ modification file (default=None)
* **msgfplus_PrecursorMassTolerance** - Precursor mass tolerance and unit
(default=10ppm)
* **msgfplus_MS2DetectorID** - MS2 mass detector used (default=1, i.e. Orbitrap)
* **msgfplus_MinCharge** - Minimum precursor charge to consider (default=2)
* **msgfplus_MaxCharge** - Maximum precursor charge to consider (default=4)
* **msgfplus_add_params** - Additional parameters for MS-GF+ as you would
indicate on command-line

Example json configuration file:
```json
{
  "threads" : 4,
  "memory" : "4G",

  "sample" : "sample/cowpox.mgf",

  "database_reference"    : "database/refseq.fasta",

  "taxonomy_mapping_ncbi" : "taxonomy/accession2taxid/prot.accession2taxid",
  "taxonomy_mapping"  : "taxonomy/accession2taxid/prot.accession2taxid",
  "taxonomy_nodes"    : "taxonomy/taxdump/nodes.dmp",
  "taxonomy_names"    : "taxonomy/taxdump/names.dmp",

  "search_engine" : "xtandem",
  "xtandem_mods_variable" : "16@M",
  "xtandem_pmmep" : 10,
  "xtandem_pmmem" : 10
}
```
