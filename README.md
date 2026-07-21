# SnakeAmp: Snakemake-based Nanopore amplicon analysis

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
passing flags directly to Snakemake. For easier execution on a cluster, a profile may be specified 
in the `profile` section of the config file.\
This mode requires Snakemake to be installed in the environment it is run in. \
**Note**: convenience features for running the pipeline on a cluster require Snakemake version 8 or higher. \
It is possible (and recommended) to use conda or mamba for managing environments.

## Setting up the pipeline

### Installation
For now the pipeline can only be installed by cloning the repository. 

### Installing the base environment
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

#### Language settings
The pipeline can produce output in Danish or English. For Danish output, set the `"language"` parameter in the config
file to `"da"`. For English output, set the `"language"` parameter in the config  file to `"en"`.

##### Input settings
The expected names of columns in the metadata "runsheet" and the output from the LIS can be configured in 
a YAML file describing input settings. A template for this can be found in `config/input_config/input_template.yaml`,
and suggested Danish and English language settings can be found in this directory as well. 

<details> <summary> Overview of available input settings </summary>
<h6>Runsheet</h6>
The runsheet is defined in the <code>run_sheet</code> section of the input config. Column names that can be set are
<ul>
<li><code>sample_number</code>: column for the sample number(s) </li>
<li><code>barcode</code>: column for the barcode used for the samples</li>
<li><code>experiment_name</code>: column in the header containing the experiment name</li>
<li><code>amplicon_type</code>: column containing the amplicon that was sequenced for the sample in question 
(e.g. 16S, 18S, RGN3)</li>
</ul>
<h6>Laboratory information system metadata (optional)</h6>
LIS metadata is defined in the <code>lis_columns</code> section of the input config. Column names that can be set are
<ul>
<li><code>sample_number</code>: field for the sample number</li>
<li><code>sample_type</code>: field for sample type identifier (as part of unique sample identifier)</li>
<li><code>isolate_number</code>: field for isolate number (as part of unique sample identifier)</li>
<li><code>patient_id</code>: field for patient identifier</li>
<li><code>date_received</code>: when the sample was received</li>
<li><code>material</code>: sample material</li>
<li><code>anatomy</code>: the anatomical site that was sampled</li>
<li><code>indication</code>: the indication for sampling</li>
</ul>
</details>


#### Laboratory information system settings
This pipeline may incorporate data from a laboratory information system (LIS) in the final report. LIS data is read from the file given under
`lis_report`.

If LIS data is given, the final report will contain additional data on each sample:
* sample date
* patient ID - note for a Danish context that this is *not* the patient's CPR number but a per-run patient ID meant to 
aid in identifying samples from the same patient in one run;
while the pipeline uses the patient's CPR number to unambiguously identify patients in this mapping
step if this is what is given in the LIS data, no mapping of patient ID to CPR number is stored by the pipeline
* sample material
* anatomical location
* indication for sampling

The pipeline was designed with OUH's MADS setup in mind. Report files therefore are 
assumed to
* be semicolon-separated
* be encoded using `latin-1` encoding
* contain fields for:
  * sample number (set under `sample_number` in the `lis_columns` section in the input settings) 
  * date received (set under `date_received` in the `lis_columns` section in the input settings)
  * patient identifier (set under `patient_id` in the `lis_columns` section in the input settings)
  * sample material/category (set under `material` in the `lis_columns` section in the input settings)
  * anatomical location (set under `anatomy` in the `lis_columns` section in the input settings)
  * indication for sampling (set under `indication` in the `lis_columns` section in the input settings)

An example of a Danish LIS report file can be found [here](docs/sample_mads.csv).  

If your laboratory information system differs, get in touch with this pipeline's maintainer(s) for help with
implementation - the idea is to make things as flexible as needed over time.

#### Defining sample number formats
In order to translate between different sample number formats in the runsheet, the output and 
optionally the LIS report, sample number formats must be described as regular expressions with named
 groups reflecting the components.

<details> <summary> Advanced sample number formats - translating between numeric and letter sample type </summary>

For sample type, if this can be marked with both a numeric code or a letter, 
it is possible to translate between these formats. This has the following requirements:
* the sample type must be defined as a named group named `sample_type` in the format specifications
where it is to be used
* the  format must be given as `"number"` for numeric codes or `"letter"` for letters
in `sample_numbers_in`, `sample_numbers_out` and `sample_numbers_output` 
in the `sample_number_settings` part of the config
* a mapping of numerical codes to letters must be given as `number_to_letter` 
in the `sample_number_settings` part of the config. For instance, if sample type F has the numeric
code 11 and sample type U has the numeric code 15, this mapping would be given as
```yaml
{"11": "F", "15": "U"}
```
</details>

<details> <summary> Advanced sample number formats - rearranging with named groups  </summary>
The use of named capturing groups allows for rearranging and dropping components as needed.
For instance:

```yaml
sample_number_settings:
  sample_number_format: '([FU]|1[15])([0-9]{8})'
  format_in_sheet: '(?P<sample_type>[FU]|1[15])(?P<sample_number>\d{6})(?P<sample_year>\d{2})'
  format_in_lis: '(?P<sample_type>[FU])(?P<sample_year>\d{2})(?P<sample_number>\d{6})'
  format_output: '(?P<sample_type>[FU]|1[15])(?P<sample_year>\d{2})(?P<sample_number>\d{6})'
```

In this example, the base `sample_number_format` stipulates that a sample number must start with
the letters F or U or the numeric code 11 or 15, followed by eight digits. \
The sample number format in the runsheet, `format_in_sheet`, breaks this up in named groups:
* `sample_type`: the letters F or U or the numeric code 11 or 15
* `sample_number`: the sample number itself, a six-digit identifier
* `sample_year`: the year component of the sample identifier (two digits)

The corresponding definition in the LIS is similar, but only allows the sample type to be given as 
letters and gives the year component *before* the sample number proper.

It is thus possible to match the sample number `F12345699` from the sample sheet to 
the sample number `F99123456` in the LIS.
</details>



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

### Input requirements
The pipeline requires the following inputs:
*  sequencing data in .fastq format (optionally compressed). One can either specify the `fastq_pass` directory containing all barcode directories, 
or a "base" directory with the default Nanopore output structure (`rawdata/*/fastq_pass/barcode*`).
* a "runsheet" with information about each sample in .xlsx format.
* and optionally:
  * a report from your laboratory information system, containing at least
    * sample number (set under `sample_number` in the `lis_columns` section in the input settings) 
    * date received (set under `date_received` in the `lis_columns` section in the input settings)
    * patient identifier (set under `patient_id` in the `lis_columns` section in the input settings)
    * sample material/category (set under `material` in the `lis_columns` section in the input settings)
    * anatomical location (set under `anatomy` in the `lis_columns` section in the input settings)
    * indication for sampling (set under `indication` in the `lis_columns` section in the input settings)

#### Minimal runsheet
As the pipeline was developed with the runsheets in use at KMA Odense in mind, the runsheet 
is expected to be an Excel sheet.

The runsheet is expected to consist of two sections:
* a three-row "header" section with the experiment name given in the third row; this column can be named by '
setting the `experiment_name` key in the input config
* a section of arbitrary length containing sample information:
  * sample number(s) (set under `sample_number` in the `run_sheet` section in the input settings)
  * the barcodes used for the samples (set under `barcode` in the `run_sheet` section in the input settings)
  Barcodes are expected to be in the format specified in the config
  (by default, RB \[rapid barcoding\] or NB \[native barcoding\] followed by the barcode number)
  * the amplicon that was sequenced for the sample in question (e.g. 16S, 18S, RGN3) 
  (set under `amplicon_type` in the `run_sheet` section in the input settings)
  The pipeline will only analyze samples with the amplicon type specified in the config under
  `amplicon_type` - if a runsheet contains samples for multiple amplicons, the pipeline will have to
  be started separately for each amplicon type with the respective config. 

An example for a runsheet can be found [here](docs/basic_runsheet.xlsx). Note that this runsheet
contains multiple amplicon types. Using a config that specifies `amplicon_type: 16S` will run the 
pipeline with samples F99123456-1 and F99123456-2 as well as a negative and positive control; with a
 config that specifies `amplicon_type: 18S`, the pipeline will  be run for samples F99123456-3 and 
a positive and negative control. \
**Note** that it is up to the user to supply the appropriate 
databases under `databases: emu_db: ...` for all amplicon types in use.

### Starting SnakeAmp before sequencing is done
SnakeAmp may be started before sequencing is complete; in this case, it will wait until 
sequencing is finished (as signaled by the presence of a `final_summary_*.txt` file) and then run 
all subsequent steps. 

## Output
The pipeline outputs a directory with filtered reads and read statistics for each sample, as well
as Emu results for all samples and a summary file. The structure is as follows:
```
+run_directory
+-sample_1
|   +-reads
+-sample_2
|   +-reads
+-...
+-emu
+-logs
+ RUN_NAME-AMPLICON_emu-combined.tsv
+ RUN_NAME-AMPLICON_emu-combined.xlsx
```
### Sample level results
* **reads**: filtered reads and QC results:
  * `SAMPLE_BARCODE.depleted.stats.tsv`: NanoStat output for reads after removal of human reads
  * `SAMPLE_BARCODE.filtered.fastq.gz`: reads remaining after filtering with filtlong and removal
  of human reads
  * `SAMPLE_BARCODE.kraken.tsv`: amount of reads mapping against the Kraken database used for 
  removing human reads
  * `SAMPLE_BARCODE.stats.tsv`: NanoStat output for the original, unfiltered reads

### Run level directories
* **emu**: Emu result files for all samples
* **logs**: sample- and pipeline-level logfiles:
  * **concat_fastq**: sample-level logs for initial concatenation of fastq files
  * **emu**: sample-level logs for taxonomic identification with emu
  * **filtlong**: sample-level logs for read filtering with filtlong
  * **kraken**: sample-level logs for human read removal with kraken
  * **snakemake.log**: Snakemake logfile for the entire run
  * **start_pipeline.log**: logfile for the start script, documenting input directory, wait time 
  until the pipeline is started, and the command used to start the pipeline

### Summary files
* `RUN_NAME-AMPLICON_emu-combined.xlsx`: a summary file combining Emu results for all samples,
as well as QC data, and, if given, information from the LIS (sample material, anatomical location, 
run-level patient identifier) for all samples. Abundance data is reported in several ways:
  * in the overview tab, both relative abundance and absolute counts are shown
  * in the abundance tab, only relative abundance is shown
  * in the count tab, only read counts are shown
* `RUN_NAME-AMPLICON_emu-combined.tsv`: a raw copy of the "overview" tab of the xlsx file

## License
The Snakeamp pipeline is released under version 3 of the GNU General Public License;
a copy of the license's text can be found [here](licenses/RSYD_BASIC_GPL).

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.