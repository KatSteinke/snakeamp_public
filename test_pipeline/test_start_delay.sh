#!/usr/bin/env bash

mkdir -p run_test/rawdata/testdir/fastq_pass
python3 /snake_qc/run_pipeline.py --dry_run \
        --runsheet /snake_qc/test_pipeline/test_nanopore_runsheet_16s_only.xlsx \
        --rundir run_test \
        --workflow_config_file /snake_qc/test_pipeline/pipeline_testjob_no_lis.yml &> delayed_start_log &
sleep 30
touch run_test/rawdata/testdir/final_summary_test.txt
sleep 30
cat delayed_start_log
grep -e "Found .+ after 30 seconds\." delayed_start_log
