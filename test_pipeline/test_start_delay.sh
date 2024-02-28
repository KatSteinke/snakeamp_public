#!/usr/bin/env bash

mkdir -p run_test/testdir/rawdata/fastq_pass
python3 /snake_qc/run_pipeline.py --runsheet /data/test_data/16S/test_data/20240124-runsheet-rearrange.xlsx \
        --rundir run_test \
        --workflow_config_file /snake_qc/config/pipeline_testjob.yml &> delayed_start_log &
sleep 30
touch run_test/testdir/rawdata/final_summary_test.txt
cat delayed_start_log
grep -e "Found pattern .+ after 30 seconds\." delayed_start_log
