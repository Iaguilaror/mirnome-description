#!/usr/bin/env bash
## This small script runs a module test with the sample data

###
## environment variable setting
export BEDFILE="test/reference/complete_mirna.bed"
###

echo "[>..] test running this module with data in test/data"
## Remove old test results, if any; then create test/reults dir
rm -rf test/results
mkdir -p test/results
echo "[>>.] results will be created in test/results"
## Execute runmk.sh, it will find the basic example in test/data
## Move results from . to test/results
## results files are *.filtered.vcf
./runmk.sh \
&& mv test/data/*.filtered.vcf test/results \
&& echo "[>>>] Module Test Successful"
