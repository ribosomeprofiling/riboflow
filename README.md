[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3376949.svg)](https://doi.org/10.5281/zenodo.3376949)

![RiboFlow](/docs/figures/riboflow_logo.jpg "RiboFlow Logo")

# RiboFlow #

RiboFlow is a [Nextflow](https://www.nextflow.io/) based pipeline
for processing ribosome profiling data. As output, it generates [ribo files](https://ribopy.readthedocs.io/en/latest/ribo_file_format.html) that can be analyzed using [RiboR](https://github.com/ribosomeprofiling/ribor) or [RiboPy](https://github.com/ribosomeprofiling/ribopy).
RiboFlow belongs to a [software ecosystem](https://ribosomeprofiling.github.io/) desgined to work with ribosome profiling data.


![Overview](/docs/figures/ecosystem_overview.jpg "Ribo Ecosystem Overview")

## Contents

* [Installation](#installation) 
* [Test Run](#test-run)  
* [Output](#output)  
* [RiboFlow on Your Data](#riboflow-on-your-data)
* [UMIs](#working-with-unique-molecular-identifiers)  
* [A Note on References](#a-note-on-references)  
* [Advanced Features](#advanced-features)  
* [Frequently Asked Questions](https://github.com/ribosomeprofiling/riboflow/blob/master/FAQ.md)  
* [Release Notes](https://github.com/ribosomeprofiling/riboflow/blob/master/CHANGELOG.md)  

## Installation

### Requirements

* [Nextflow](https://www.nextflow.io/)
* [Docker](https://docs.docker.com/install/) (Optional)
* [Conda](https://conda.io/en/latest/miniconda.html) (Optional)

First, follow the instructions in [Nextflow website](https://www.nextflow.io/) and install Nextflow.

### Docker Option

Install [Docker](https://docs.docker.com/install/).
Here is a [tutorial for Ubuntu.](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04)

All remaining dependencies come in the Docker image [hakanozadam/riboflow](https://hub.docker.com/r/hakanozadam/riboflow).
This image is automatically pulled by RiboFlow when run with Docker (see test runs below).

### Conda Option

This option has been tested on Linux systems only.

Install  [Conda](https://conda.io/en/latest/miniconda.html).

All other dependencies can be installed using the environment file,
environment.yaml, in this repository.
```
git clone https://github.com/ribosomeprofiling/riboflow.git
conda env create -f riboflow/environment.yaml
```

The above command will create a conda environment called _ribo_
and install dependencies in it.
To start using RiboFlow, you need to activate the _ribo_ environment.

`conda activate ribo`

## Test Run

For fresh installations, before running RiboFlow on actual data,
it is recommended to do a test run.


### Run Using Docker

```
# Clone this repository in a new folder and change your working directory to the RiboFlow folder.
mkdir rf_test_run && cd rf_test_run
git clone https://github.com/ribosomeprofiling/riboflow.git
cd riboflow

# Obtain a copy of the sample data in the working directory.
git clone https://github.com/ribosomeprofiling/rf_sample_data.git
nextflow RiboFlow.groovy -params-file project.yaml -profile docker_local
```

Note that we provided the argument `-profile docker_local` to Nextflow to indicate that RiboFlow will be run via Docker containers. In other words, the steps of RiboFlow will be executed inside Docker containers by Nextflow. 
Hence, no locally installed software (other than Java and Nextflow) is needed by RiboFlow.  


### Run Using Conda Environment

In Conda option, the steps of RiboFlow are run locally. So, we need to install the dependencies first. This can easily be done via conda. The default profile directs RiboFlow to run locally, so we can simply skip the `-profile` argument. Also note that the conda environment has to be activated before running RiboFlow. 

Before running the commands below, make sure that you have created the conda environment, called _ribo_,
using the instructions above. 

```
# List the environments to make sure that ribo environment exists
conda env list

# Activate the ribo environment
conda activate ribo

# Get RiboFlow repository
mkdir rf_test_run && cd rf_test_run
git clone https://github.com/ribosomeprofiling/riboflow.git
cd riboflow

# Obtain a copy of the sample data in the working directory.
git clone https://github.com/ribosomeprofiling/rf_sample_data.git

# Finally run RiboFlow
nextflow RiboFlow.groovy -params-file project.yaml

```

## Output

Pipeline run may take several minutes.
When finished, the resulting files are in the `./output` folder.

Mapping statistics are compiled in a csv file called `stats.csv`

```
ls output/stats/stats.csv
```

Ribosome occupancy data is in a single
[ribo file](https://ribopy.readthedocs.io/en/latest/ribo_file_format.html) called `all.ribo`.

`ls output/ribo/all.ribo`

You can use
[RiboR](https://github.com/ribosomeprofiling/ribor) or
[RiboPy](https://github.com/ribosomeprofiling/ribopy) to work with ribo files.


## RiboFlow on Your Data

For running RiboFlow on actual data, files must be organized and a parameters file must be prepared.
You can examine the sample run above to see an example.

1. Organize your data. The following files are required for RiboFlow
   * **Ribosome profiling sequencing data:** in gzipped fastq files
   * **Transcriptome Reference:** Bowtie2 index files
   * **Filter Reference:** Bowtie2 index files (typically for rRNA sequences)
   * **Annotation:** A bed file defining CDS, UTR5 and UTR3 regions.
   * **Transcript Lengths:** A two column tsv file containing transcript lengths

2. Prepare a custom `project.yaml` file.
You can use the sample file `project.yaml`, provided in this repository,
as template.

3. In `project.yaml`, provide RiboFlow parameters such as `clip_arguments`, alignment arguments etc.
You can simply modify the arguments in the sample file `project.yaml` in this repository.

4. You can adjust the hardware and computing environment settings in Nextflow configuration file(s).
For Docker option, see `configs/docker_local.config`. If you are not using Docker,
see `configs/local.config`.

5. RNA-Seq data is optional for RiboFlow. So, if you do NOT have RNA-Seq data, in the project file, set

`do_rnaseq: false`

If you have RNA-Seq data to be paired with ribosome profiling data, see the [Advanced Features](#advanced-features) below.


6. Metadata is optional for RiboFlow. If you do NOT have metadata, in the project file, set

`do_metadata: false`

If you have metadata, see [Advanced Features](#advanced-features) below.

7. Run RiboFlow using the new parameters file `project.yaml`.

Using Docker:

`nextflow RiboFlow.groovy -params-file project.yaml -profile docker_local`

Without Docker:

`nextflow RiboFlow.groovy -params-file project.yaml`

## Working with Unique Molecular Identifiers
Unique Molecular Identifiers (UMIs) can be ligated to either side of the molecules and 
they allow labeling molecules uniqely. This way UMIs can be used to deduplicate mapped reads
for more accurate quantification.

If there are UMIs in your ribosome profiling data, Riboflow can trim them and deduplicate reads based on UMIs. 

RiboFlow extracts UMIs and stores them in the Fastq headers and uses the UMIs in deduplication 
(instead of position based read collapsing). For this purpose RiboFlow uses 
[umi_tools](https://github.com/CGATOxford/UMI-tools).


### Project File

Here we explain the related parts of the project file to be able to use UMIs feature of Riboflow.

Also, we provide a working example of project file in this repository: *project_umi.yaml*.

The following parameter must be set:
```
dedup_method: "umi_tools"
```

Also, users must set the following two parameters: `umi_tools_extract_arguments` and `umi_tools_dedup_arguments`.

For example: 
```
umi_tools_extract_arguments: "-p \"^(?P<umi_1>.{12})(?P<discard_1>.{4}).+$\" --extract-method=regex"
umi_tools_dedup_arguments:   "--read-length"
```

The above example takes the first 12 nucleotides from the 5' end, discards the 4 nucleotides downstream and writes the 12 nt UMI sequence to the header.
The second parameter tells umi_tools to use read lengths IN ADDITION to UMI sequencing in collapsing reads. Note that these two arguments are directly provided to umi_tools. So users are encouraged to familirize themselves with [umi_tools](https://umi-tools.readthedocs.io/en/latest/).

### Test Run with UMIs

We provide a mini dataset, with two samples, to try Riboflow with sequencing reads having UMIs.
In this sample dataset, the first 12 nucleotides on the 5' end of the reads are UMIs.
Four nucleotides following the UMIs need to be discarded.
On the 3' end of the reads, there are adapters having the sequence `AAAAAAAAAACAAAAAAAAAA`.
The parameters of this sample run are provided in the file *project_umi.yaml*.
Below are the steps to process this data.


```
# List the environments to make sure that ribo environment exists
conda env list

# Activate the ribo environment
conda activate ribo

# Get RiboFlow repository
mkdir rf_test_run && cd rf_test_run
git clone https://github.com/ribosomeprofiling/riboflow.git
cd riboflow

# Obtain a copy of the sample data in the working directory.
git clone https://github.com/ribosomeprofiling/rf_sample_data.git

# Finally run RiboFlow
nextflow RiboFlow.groovy -params-file project_umi.yaml

# At the end of the run
# checkk the ribo file
ribopy info output_umi/ribo/all.ribo
```

 ### UMI support for RNA-Seq
 
 In the current version, UMIs are supported for ribosome profiling data only. So RNA-Seq libraries can either be used without deduplication or the reads can be collapsed based on position.

## A Note on References

RiboFlow is designed to work with transcriptomic references. RiboFlow does **NOT** work with genomic references.
The users need to provide a transcriptome reference and annotation to run this software.
There is a curated set of RiboFlow references, that users can download and use, in
[this GitHub repository](https://github.com/ribosomeprofiling/references_for_riboflow)

## Advanced Features

### RNA-Seq Data

If you have RNA-Seq data that you want to pair with ribosome profiling experiments,
provide the paths of the RNA-Seq (gzipped) fastq files  in the configuration file in
_input -> metadata_. See the file `project.yaml` in this repository for an example.
Note that the names in defining RNA-Seq files must match the names in definig ribosome profiling data.
Also turn set the do_rnaseq flag to true, in the project file:

`do_rnaseq: true`

Transcript abundance data will be stored in the output ribo file.

### Metadata

If you have metadata files for the ribosome profiling experiments,
provide the paths of the metadata files (in yaml format) in the configuration file in
_input -> metadata_. See the file `project.yaml` in this repository for an example.
Note that the names in defining metadata files must match the names in definig ribosome profiling data.
Also turn set the metadata flag to true, in the project file:

`do_metadata: true`

Metadata will be stored in the output ribo file.

## Citing

[RiboFlow, RiboR and RiboPy: an ecosystem for analyzing ribosome profiling data at read length resolution, H. Ozadam, M. Geng, C. Cenik
Bioinformatics 36 (9), 2929-2931](https://academic.oup.com/bioinformatics/article/36/9/2929/5701654)

```bibtex
@article{ozadam2020riboflow,
  title={RiboFlow, RiboR and RiboPy: an ecosystem for analyzing ribosome profiling data at read length resolution},
  author={Ozadam, Hakan and Geng, Michael and Cenik, Can},
  journal={Bioinformatics},
  volume={36},
  number={9},
  pages={2929--2931},
  year={2020},
  publisher={Oxford University Press}
}
```

## [Frequently Asked Questions](https://github.com/ribosomeprofiling/riboflow/blob/master/FAQ.md)  

  
## [Release Notes](https://github.com/ribosomeprofiling/riboflow/blob/master/CHANGELOG.md)  
