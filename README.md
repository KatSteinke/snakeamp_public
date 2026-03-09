# SnakeAmp: Snakemake-based Nanopore amplicon analysis
<!--- # (Name TBD; may find a punny acronym so long as it starts with S. Snatched? taxonomy, clinical, EMU) --->

This is a pipeline for amplicon-based taxonomic assignment on the basis of Nanopore sequences.
The core functionality is based around cleaning input reads and assigning taxonomy using Emu,
then creating a report on the basis of Emu results. Other features such as improved contaminant
depletion and de novo clustering are under development.

The current maintainer is Kat Steinke at KMA Odense;
the structure of the pipeline is based on a script created by Flemming Damgaard.

## Prerequisites
The pipeline can run in two setups:
* Nomad mode: accessed by setting `run_on: nomad` in the config file. 
The pipeline runs as a Nomad job (to be defined in `nomad.hcl`) in a Docker container with
prebuilt conda environments. \
Implementing this mode requires Nomad infrastructure.
* Local mode: accessed by setting `run_on: local` in the config file.
The pipeline runs in Snakemake, by default locally on the computer/server/node the `run_pipeline.py`
script is executed on. Additional configuration is possible, both through the config file and by 
passing flags directly to Snakemake.\
This mode requires Snakemake to be installed in the environment it is run in. \
**Note**: running the pipeline on a cluster is currently only compatible with Snakemake version 7.*
or lower. \
It is possible (and recommended) to use conda or mamba for managing environments.

## Setting up the pipeline

### Installation
For now the pipeline can only be installed by cloning the repository. 

### Installing the base environment - TODO
The base environment for launching the pipeline is defined as `base_env.yml`. This file can be used
to create the base environment using mamba or conda. Navigate to the directory containing this 
repository for the next step. 

#### Using mamba
```
mamba env create --file envs/base_env.yml
```

#### Using conda
```
conda env create --file envs/base_env.yml
```
### Configuring the pipeline
A default config file with placeholders is provided as `pipeline_routine.yaml`; replace 
the placeholders and you're ready to go. For testing etc., it is recommended to create separate 
config files. 

TODO: should we have input config as with RSYD-BASIC?

#### Laboratory information system settings
This pipeline may incorporate data from a laboratory information system (LIS) in the final report. LIS data is read from the file given under
`lis_report`.

If LIS data is given, the final report will contain additional data on each sample:
* sample date ("prøvetagningsdato")
* patient ID - note that this is *not* the patient's CPR number but a per-run patient ID meant to 
aid in identifying samples from the same patient in one run; while the pipeline bases this on 
the patient's CPR number, no mapping of patient ID to CPR number is stored by the pipeline
* sample material ("prøvemateriale")
* anatomical location ("anatomi")
* indication for sampling ("indikation")

The pipeline was designed with OUH's MADS setup in mind. Report files therefore are 
assumed to
* be semicolon-separated
* be encoded using `latin-1` encoding
* contain fields for:
  * date received ("modtaget")
  * patient CPR number ("cprnr.")
  * sample material/category ("prøvekategori")
  * anatomical location ("anatomi")
  * indication for sampling ("Indikation")

An example of a LIS report file can be found [here](docs/sample_mads.csv). TODO  

If your laboratory information system differs, get in touch with this pipeline's maintainer(s) for help with
implementation - the idea is to make things as flexible as needed over time.

## Running the pipeline
The pipeline can be run in two modes: "classic" mode, primarily meant for users not 
comfortable with running commands on the commandline, and "commandline" mode, which 
offers greater control over the pipeline.
### Activating the environment
Before running the pipeline, the base environment needs to be activated. This is done by running
```shell
conda activate snake_amp
```

### Classic mode
Classic mode can be started by running 
```shell
python3 run_pipeline.py
```
in the virtual environment containing all prerequisites. \
Optionally, a configuration file may be specified:
```shell
python3 run_pipeline.py --workflow_config_file config/pipeline_18s.yml
```
**Note** that this is the only flag that may be given in classic mode. 


Classic mode guides the user through the process in a "questionnaire" style. \
As it is meant for routine use by non-expert users, some options cannot be changed in classic mode:
* **continuing a stopped run**: a stopped run cannot be continued in classic mode, only in commandline mode
* **additional flags for Snakemake**: passing specific flags to Snakemake is only possible in commandline mode
* **logfile location**: logfiles for starting the pipeline will always be saved in the default
location (`logs/start_pipeline.log` in the output directory) in classic mode
* **test vs routine Nomad jobs**: in Nomad mode, starting the pipeline in classic mode will always 
follow the config file's settings for what job to dispatch, so the `16s-snake-emu-staging` job
will only be dispatched if `debug` is set to `True` in the config file given to the pipeline.
### Commandline mode
The minimal input for commandline mode is
```shell
python3 run_pipeline.py --rundir path/to/rundir --runsheet path/to/runsheet
```
where `---rundir` is the directory containing sequencing results and
`--runsheet` is the path to the runsheet (see below). \
This runs the pipeline with the same default settings as in classic mode. Additionally, the 
following options can be specified:
* `--workflow_config_file`: the config file from which to read settings for the run
* `--outdir`: output directory (also possible to specify in classic mode when prompted)
* `--run_time`: the expected maximum time for sequencing to finish, in hours
(also possible to specify in classic mode when prompted)
* `--continue_pipeline`: flag to specify continuing a previously stopped run
* `--snake_flags`: flags to pass to Snakemake, enclosed in quotes; 
  note that if only a single flag is given, a space must be appended, e.g. `--snake_flags "-n "`
* `--logfile`: the path where to store logs for pipeline start
* `--dry_run`: instead of running the pipeline, determines the pipeline start command,
sets up required dirs and exits
* `--test_run`: only relevant when running on nomad - always dispatches `16s-snake-emu-staging`
  regardless of config settings

### Input formats
The pipeline requires the following inputs:
*  sequencing data in .fastq format (optionally compressed). One can either specify the `fastq_pass` directory containing all barcode directories, 
or a "base" directory with the default Nanopore output structure (`rawdata/*/fastq_pass/barcode*`).
* a "runsheet" with information about each sample in .xlsx format.
* and optionally:
  * a report from your laboratory information system, containing at least
    * sample numbers (TODO) 
    * date received ("modtaget")
    * patient CPR number ("cprnr.")
    * sample material/category ("prøvekategori")
    * anatomical location ("anatomi")
    * indication for sampling ("Indikation")

#### Minimal runsheet
As the pipeline was developed with the runsheets in use at KMA Odense in mind, the runsheet 
is expected to be an Excel sheet. 
TODO description

<!--- TODO runsheet An example for an Illumina runsheet can be found [here](docs/minimal_runsheet.xlsx). --->

Config: use_conda, frontend, prefix
snake_flags to be passed through to snakemake - currently used e.g. to specify a snakemake `--profile`