import pathlib

import helpers
import snake_helpers

workdir: config["outdir"]

# set relevant dirs

RUNDIR = pathlib.Path(config["rundir"])
FASTQ_DIR = helpers.get_fastq_pass_dir(RUNDIR)

wildcard_constraints:
    barcode_number=r"\d{2}"
# TODO: we can absolutely solve this better - runsheets or such - use what's in place or have a new one?
BARCODES = glob_wildcards(f"{FASTQ_DIR}/barcode{{barcode_number}}/*").barcode_number
print(BARCODES)



rule all:
    input:
        all_results = "emu/emu-combined-tax_id.tsv"

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
        concat_fasta = temp("barcode{barcode_number}/reads/barcode{barcode_number}.reads.fastq")
    message: f"# Concatenating fastq files for barcode {{wildcards.barcode_number}}...."
    log:
        "logs/concat_fastq/barcode{barcode_number}_log.txt"
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

rule clean_nanopore_reads:
    input:
        concat_fastq = "barcode{barcode_number}/reads/barcode{barcode_number}.reads.fastq"
    output:
        filtered_fastq = temp("barcode{barcode_number}/reads/barcode{barcode_number}.filtered.fastq")
    params:
        min_length = 100
    conda: "nanopore_qc_env" # TODO: set up env!
    shell:
        """
        filtlong --min_length {params.min_length} --keep_percent 95 {input.concat_fastq} >  {output.filtered_fastq}
        """


rule compress_nanopore_reads:
    input:
        filtered_fastq = "barcode{barcode_number}/reads/barcode{barcode_number}.filtered.fastq"
    output:
        compressed_fastq = "barcode{barcode_number}/reads/barcode{barcode_number}.filtered.fastq.gz"
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
        compressed_fastq = "barcode{barcode_number}/reads/barcode{barcode_number}.filtered.fastq.gz"
    output:
        fasta_reads = "barcode{barcode_number}/reads/barcode{barcode_number}.filtered.fasta"
    conda:
        "nanopore_qc_env"
    shell:
        """
        seqtk seq -a "{input.compressed_fastq}" > "{output.fasta_reads}"
        """

# TODO: logging?
rule run_emu:
    input:
        fasta_reads = "barcode{barcode_number}/reads/barcode{barcode_number}.filtered.fasta"
    output:
        relative_abundance = "emu/barcode{barcode_number}_rel-abundance.tsv"
    params:
        emu_db = config["databases"]["emu_db"],
        outdir = lambda wildcards, output: str(pathlib.Path(output.relative_abundance).parent),
        basename = "barcode{barcode_number}"
    conda:
        "emu_env"
    threads: (workflow.cores / 4 ) if (workflow.cores / 4 ) <= 64 else 64
    shell:
        """
        emu abundance "{input.fasta_reads}" --db "{params.emu_db}" --keep-counts \
         --output-dir "{params.outdir}" --output-basename {params.basename} \
         --threads {threads}
        """

rule combine_emu:
    input:
        all_relative_abundance = expand("emu/barcode{barcode_number}_rel-abundance.tsv",
                                        barcode_number = BARCODES)
    output:
        counts_combined = "emu/emu-combined-tax_id.tsv"
    params:
        emu_dir = "emu",
        tax_rank = "tax_id"
    conda:
        "emu_env"  # minmap >= 2.22
    shell:
        """
        emu combine-outputs "{params.emu_dir}" {params.tax_rank}
        """
