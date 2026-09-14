#!/usr/bin/env bash

mkdir -p run_test/rawdata/testdir/fastq_pass/barcode01
touch run_test/rawdata/testdir/fastq_pass/barcode01/test.fastq.gz
mkdir -p run_test/rawdata/testdir/fastq_pass/barcode02
touch run_test/rawdata/testdir/fastq_pass/barcode02/test.fastq.gz
mkdir -p run_test/rawdata/testdir/fastq_pass/barcode03
touch run_test/rawdata/testdir/fastq_pass/barcode03/test.fastq.gz
mkdir -p run_test/rawdata/testdir/fastq_pass/barcode42
touch run_test/rawdata/testdir/fastq_pass/barcode42/test.fastq.gz
touch run_test/rawdata/testdir/final_summary_test.txt
snakemake --version # is snakemake even alive?
python3 /snake_qc/run_pipeline.py \
        --runsheet /snake_qc/test_pipeline/test_nanopore_runsheet_16s_only.xlsx \
        --rundir run_test/rawdata/testdir/fastq_pass \
        --workflow_config_file /snake_qc/test_pipeline/pipeline_testjob_no_lis_local.yml \
         --snake_flags "--dry-run " > >(tee -a local_log) 2> >(tee -a local_log) # don't save this at all for now so we can see what's happening
# do we need to wait here?
sleep 30
echo "printing local log"
cat local_log
grep -e "This was a dry-run (flag -n)" local_log