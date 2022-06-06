#!/usr/bin/env bash
## This small script runs a module test with the sample data

###
## environment variable setting
export BEDDIR="test/data"
###

echo "[>..] test running this module with data in test/data"
## Remove old test results, if any; then create test/reults dir
rm -rf test/results
mkdir -p test/results
echo "[>>.] results will be created in test/results"
## Execute runmk.sh, it will find the basic example in test/data
## Move results from . to test/results
## results files are *genome1_vs_genome2_read_proportions*
./runmk.sh \
&& mv allvariants.bed variants_of_interest_by_type.csv boxplots.svg test/results \
&& echo "[>>>] Module Test Successful"
