import pathlib
import re

import pandas as pd

import helpers
import pipeline_config
import snake_helpers

# set default config if none is provided - will be overridden in most setups
configfile: pipeline_config.default_config_file
# path to config file needs to be specified for other scripts
CONFIG_PATH = config["config_path"] if "config_path" in config \
    else pipeline_config.default_config_file


workdir: config["outdir"]

# set relevant dirs

RUNDIR = pathlib.Path(config["rundir"])
FASTQ_DIR = helpers.get_fastq_pass_dir(RUNDIR)
print(FASTQ_DIR)
sample_number_pattern = helpers.get_id_pattern(config["sample_number_settings"]["sample_number_format"],
        negative_control = config["sample_number_settings"]["negative_control"],
        positive_control = config["sample_number_settings"]["positive_control"])
sample_number_pattern = re.compile(f"^{sample_number_pattern.pattern}$")
input_format = re.compile(config["sample_number_settings"]["format_in_sheet"])
output_format = re.compile(config["sample_number_settings"]["format_output"])
prefix_translate = helpers.get_number_letter_combination(config["sample_number_settings"]["number_to_letter"],
                                                         config["sample_number_settings"]["sample_numbers_in"],
                                                         config["sample_number_settings"]["sample_numbers_output"])
(positive_control,
 negative_control) = helpers.get_control_patterns(config["sample_number_settings"]["negative_control"],
                                                  config["sample_number_settings"]["positive_control"])
wildcard_constraints:
    barcode_number = r"\d{2}",
    barcode = config["barcode_format"],
    #sample_number = sample_number_pattern
# TODO: we can absolutely solve this better - runsheets or such - use what's in place or have a new one?

sheet_data = pd.read_excel(config["runsheet"],usecols = "A:D",skiprows = 3,
                           dtype = {"Prøvenummer": str, "Eluat nr.": str})
sheet_data = sheet_data.dropna(subset=["Prøvenummer", "Barkode"])
sheet_data = sheet_data[sheet_data["Analyse"] == config["amplicon_type"]]
sheet_data["prøvenr"] = sheet_data["Prøvenummer"].apply(lambda sample_number:
                                                       helpers.translate_sample_number(sample_number,
                                                                                       input_format,
                                                                                       output_format,
                                                                                       prefix_translate,
                                                                                       positive_control,
                                                                                       negative_control))

BARCODE_PREFIX = config["barcode_prefix"]
ALL_IDS = list(sheet_data["prøvenr"])
print(ALL_IDS)
print(sheet_data["prøvenr"].str.match(sample_number_pattern, na=False))
ALL_BARCODES = list(sheet_data["Barkode"])

# we need to name some files after the experiment name
EXPERIMENT_NAME = helpers.extract_nanopore_run_name(pathlib.Path(config["runsheet"]))

rule all:
    input:
        all_results = f"{EXPERIMENT_NAME}_emu-combined.xlsx"  # TODO: experiment name!

rule concatenate_fastqs:
    params:
        barcode_dir = f"{FASTQ_DIR}/barcode{{barcode_number}}",
        file_format = lambda wildcards: "fastq.gz" if snake_helpers.is_gzipped(FASTQ_DIR,
                                                                               wildcards.barcode_number)
                                                    else "fastq",
        concatenate = lambda wildcards: "zcat" if snake_helpers.is_gzipped(FASTQ_DIR,
                                                                           wildcards.barcode_number)
                                               else "cat"
    output:
        concat_fasta = temp(f"{{sample_number}}_{BARCODE_PREFIX}{{barcode_number}}/reads/"
                            f"{{sample_number}}_{BARCODE_PREFIX}{{barcode_number}}.reads.fastq")
    message: f"# Concatenating fastq files for barcode {{wildcards.barcode_number}}...."
    log:
        f"logs/concat_fastq/{{sample_number}}_{BARCODE_PREFIX}{{barcode_number}}_log.txt"
    conda:
        "nanopore_qc_env"
    resources:
        mem_mb = 200
    shell:
        """
         find \
         {params.barcode_dir} \
            -type f \
            -name "*{params.file_format}" \
            -exec {params.concatenate} {{}} \; \
            > {output.concat_fasta} || touch {output.concat_fasta}

         if [ ! -s {output.concat_fasta} ]; then
            echo "No reads were found for {wildcards.barcode_number}. No consensus is generated."
            echo "No reads were found for {wildcards.barcode_number}. No consensus is generated." > {log}
         fi
        """

rule remove_human_reads:
    input:
        concat_fasta = "{sample_number}_{barcode}/reads/{sample_number}_{barcode}.reads.fastq"
    output:
        human_depleted = temp("{sample_number}_{barcode}/reads"
                              "/{sample_number}_{barcode}.depleted.fastq")
    params:
        kraken_db = pathlib.Path(config['databases']['human_reads']),
    conda: "kraken_env"
    resources:
        mem_mb = 5000  # database + a bit extra
    threads: workflow.cores
    shell:
        """
        kraken2 --db "{params.kraken_db}" --unclassified-out "{output.human_depleted}" \
        --output "-" --threads {threads} \
        {input.concat_fasta}
        """


rule clean_nanopore_reads:
    input:
        concat_fastq = "{sample_number}_{barcode}/reads/" \
                       "{sample_number}_{barcode}.depleted.fastq"
    output:
        filtered_fastq = temp("{sample_number}_{barcode}/reads/"
                              "{sample_number}_{barcode}.filtered.fastq")
    params:
        min_length = 100
    conda: "nanopore_qc_env" # TODO: set up env!
    shell:
        """
        filtlong --min_length {params.min_length} --keep_percent 95 {input.concat_fastq} >  {output.filtered_fastq}
        """


rule compress_nanopore_reads:
    input:
        filtered_fastq = "{sample_number}_{barcode}/reads/" \
                         "{sample_number}_{barcode}.filtered.fastq"
    output:
        compressed_fastq = "{sample_number}_{barcode}/reads/" \
                           "{sample_number}_{barcode}.filtered.fastq.gz"
    conda:
        "nanopore_qc_env"  # TODO: needs to have pigz
    threads: 2
    resources:
        mem_mb = 100
    shell:
        """
        pigz -p {threads} -c -n "{input.filtered_fastq}" > "{output.compressed_fastq}"
        """

rule fastq_to_fasta:
    input:
        compressed_fastq = "{sample_number}_{barcode}/reads/" \
                           "{sample_number}_{barcode}.filtered.fastq.gz"
    output:
        fasta_reads = "{sample_number}_{barcode}/reads/" \
                      "{sample_number}_{barcode}.filtered.fasta"
    conda:
        "nanopore_qc_env"
    shell:
        """
        seqtk seq -a "{input.compressed_fastq}" > "{output.fasta_reads}"
        """

rule run_emu:
    input:
        fasta_reads = "{sample_number}_{barcode}/reads/" \
                      "{sample_number}_{barcode}.filtered.fasta"
    output:
        relative_abundance = f"emu/{EXPERIMENT_NAME}_{{sample_number}}_{{barcode}}_rel-abundance.tsv"
    params:
        emu_db = config["databases"]["emu_db"],
        outdir = lambda wildcards, output: str(pathlib.Path(output.relative_abundance).parent),
        basename = f"{EXPERIMENT_NAME}_{{sample_number}}_{{barcode}}",
        # add very minimal results if emu fails
        fallback_header = r"tax_id\tabundance\testimated_counts\n"
    conda:
        "emu_env"
    threads: (workflow.cores / 4 ) if (workflow.cores / 4 ) <= 64 else 64
    log:
        "logs/emu/{sample_number}_{barcode}.log"
    shell:
        """
        emu abundance "{input.fasta_reads}" --db "{params.emu_db}" --keep-counts \
         --output-dir "{params.outdir}" --output-basename {params.basename} \
         --threads {threads} &> "{log}" || {{ printf "{params.fallback_header}" > "{output.relative_abundance}" ; \
          printf "unassigned\\t0.0\\t$(grep -P '(?<=Unassigned read count: )[0-9]+' {log:q} --only-matching)\\n" ; }}
        """

rule combine_emu:
    input:
        all_relative_abundance = expand(f"emu/{EXPERIMENT_NAME}_{{sample_number}}_{{barcode}}_rel-abundance.tsv",
                                        zip,
                                        sample_number=ALL_IDS, barcode=ALL_BARCODES)
    output:
        counts_combined = f"{EXPERIMENT_NAME}_emu-combined.xlsx",
        counts_raw = f"{EXPERIMENT_NAME}_emu-combined.tsv"
    params:
        emu_dir = "emu",
        basedir = workflow.current_basedir,
        configfile = CONFIG_PATH
    log:
        "logs/emu/combine_all.log"
    shell:
        """
        python3 {params.basedir}/summarize_emu.py "{params.emu_dir}" \
         --outfile "{output.counts_combined}" \
         --outfile_raw "{output.counts_raw}" \
         --workflow_config_file "{params.configfile}" &> "{log}"
        """

onsuccess:
    shell('mkdir -p logs; cat "{log}" >> "logs/snakemake.log"')

onerror:
    shell('mkdir -p logs; cat "{log}" >> "logs/snakemake.log"')
