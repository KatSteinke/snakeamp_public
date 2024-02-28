#!/usr/bin/env bash

mkdir -p run_test/testdir/rawdata/fastq_pass
python3 /snake_qc/run_pipeline.py --dry_run \
        --runsheet /snake_qc/test/data/sample_sheet_test/test_translate_runsheet_rearrange.xlsx \
        --rundir run_test \
        --workflow_config_file /snake_qc/config/pipeline_testjob.yml &> delayed_start_log &
sleep 30
touch run_test/testdir/rawdata/final_summary_test.txt
cat delayed_start_log
grep -e "Found pattern .+ after 30 seconds\." delayed_start_log
