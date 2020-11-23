#!/bin/bash

source /opt/ferramentas/scripts/setup.synopsys

export ISPD_CONTEST_ROOT=/home/gaflach/benchmarks/sizing/ispd2013
export ISPD_CONTEST_BENCHMARK=$1

pt_shell -f ./pt_load_scripts.tcl

RESULT_SLACK=$(grep "Constraint" -A 2 out.txt | grep "max_delay/setup" | awk '{print $2}')

echo "Benchmark: $1"
echo "Worst Slack: $RESULT_SLACK"
