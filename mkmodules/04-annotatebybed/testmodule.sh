#!/usr/bin/env bash
## This small script runs a module test with the sample data

###
## environment variable setting
export BASEBED="test/reference/samplechr22_76g_PASS.filtered.bed"
###

echo "[>..] test running this module with data in test/data"
## Remove old test results, if any; then create test/reults dir
rm -rf test/results
mkdir -p test/results
echo "[>>.] results will be created in test/results"
## Execute runmk.sh, it will find the basic example in test/data
## Move results from . to test/results
## results files are *.variants.vcf
./runmk.sh \
&& mv test/data/*.variants.bed test/results \
&& echo "[>>>] Module Test Successful"
