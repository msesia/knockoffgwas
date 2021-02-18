#!/bin/bash
#
# Generate knockoff negative-controls
#
# Authors: Matteo Sesia
# Date:    07/21/2020

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

# Temporary storage of intermediate files
TMP_DIR="../tmp"
mkdir -p $TMP_DIR
mkdir -p $TMP_DIR"/knockoffs"

# Storage of output files
OUT_DIR="../results"
mkdir -p $OUT_DIR

# List of chromosomes
CHR_MIN=21
CHR_MAX=22

# Path to snpknock2 executable built as described above
SNPKNOCK2="../snpknock2/bin/snpknock2"

# Which operations should we perform?
FLAG_GENERATE_KNOCKOFFS=1

######################
# Generate knockoffs #
######################

if [[ $FLAG_GENERATE_KNOCKOFFS == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Generating knockoffs"
  echo "----------------------------------------------------------------------------------------------------"

  # Run snpknock2
  $SNPKNOCK2 \
    --bgen "../data/haplotypes/example_chr{"$CHR_MIN":"$CHR_MAX"}" \
    --keep "../data/qc/samples_qc.txt" \
    --extract "../data/qc/qc_chr{"$CHR_MIN":"$CHR_MAX"}.txt" \
    --map "../data/maps/genetic_map_chr{"$CHR_MIN":"$CHR_MAX"}.txt" \
    --part "../tmp/partitions/example_chr{"$CHR_MIN":"$CHR_MAX"}.txt" \
    --ibd "../data/ibd/example_chr{"$CHR_MIN":"$CHR_MAX"}.txt" \
    --K 10 \
    --cluster_size_min 1000 \
    --cluster_size_max 10000 \
    --hmm-rho 1 \
    --hmm-lambda 1e-3 \
    --windows 0 \
    --n_threads 4 \
    --seed 2020 \
    --compute-references \
    --generate-knockoffs \
    --out $TMP_DIR"/knockoffs/knockoffs_chr{"$CHR_MIN":"$CHR_MAX"}"

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping generation of knockoffs"
  echo "----------------------------------------------------------------------------------------------------"
fi
