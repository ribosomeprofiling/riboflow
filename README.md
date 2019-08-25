# RiboFlow

RiboFlow is a [Nextflow](https://www.nextflow.io/) based pipeline 
for processing ribosome profiling data.

## Installation

### Requirements

* [Nextflow](https://www.nextflow.io/)
* [Docker](https://docs.docker.com/install/) (Optional) 
* [Conda](https://conda.io/en/latest/miniconda.html) (Optional)

First, follow the instructions in [Nextflow website](https://www.nextflow.io/) and install Nextflow. 

The easiest way of using RiboFLow is using Docker.
If using Docker is not an option, you can install the dependencies using Conda
and run RiboFlow without Docker. 

### Docker Option

Install [Docker](https://docs.docker.com/install/). 
Here is a [tutorial for Ubuntu.](https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04)

All remaining dependencies come in the Docker image [ceniklab/riboflow](https://hub.docker.com/r/ceniklab/riboflow).
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

Clone this repository in a new folder and change your working directory to the RiboFlow folder. 
```
mkdir rf_test_run && cd rf_test_run
git clone https://github.com/ribosomeprofiling/riboflow.git
cd riboflow
```

Obtain a copy of the sample data in the working directory.
```
git clone https://github.com/ribosomeprofiling/rf_sample_data.git
```

### Run Using Docker

Provide the argument `-profile docker_local` to Nextflow to indicate Docker use. 

`nextflow RiboFlow.groovy -params-file project.yaml -profile docker_local`

### Run Using Conda Environment

Make sure that you have created the conda environment, called _ribo_,
using the instructions above. Then activate the conda environment.

`conda activate ribo` 

If the above command fails to activate the ribo environment, try
`source activate ribo`
 
Now RiboFlow is ready to run.

`nextflow RiboFlow.groovy -params-file project.yaml`

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


## Actual Run

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

If you have RNA-Seq data to be paired with ribosome profiling data, see the __Advanced Features__ below.


6. Metadata is optional for RiboFlow.. If you do NOT have metadata, in the project file, set

`do_metadata: false`

If you have metadata, see Advanced feature below.

7. Run RiboFlow using the new parameters file `project.yaml`.

Using Docker:
 
`nextflow RiboFlow.groovy -params-file project.yaml -profile docker_local`

Without Docker:

`nextflow RiboFlow.groovy -params-file project.yaml`

## Advanced Features

### RNA-Seq Data

If you have RNA-Seq data that you want to pair with ribosome profiling experiments,
provide the paths of the RNA-Seq (gzipped) fastq files  in the configuration file in
_input -> metadata_. See the file `project.yaml` in this repository for an example.
Note that the names in defining RNA-Seq files must match the names in definig ribosome profiling data.
Also turn set the do_rnaseq flag to true, in the project file:

`do_rnaseq: true`

Transcript abundance data will be stored in 

### Metadata

If you have metadata files for the ribosome profiling experiments,
provide the paths of the metadata files (in yaml format) in the configuration file in
_input -> metadata_. See the file `project.yaml` in this repository for an example.
Note that the names in defining metadata files must match the names in definig ribosome profiling data.
Also turn set the metadata flag to true, in the project file:

`do_metadata: true`

Metadata will be stored in the output ribo file.

