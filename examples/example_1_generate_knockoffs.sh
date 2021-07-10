#!/bin/bash

# If you haven't already done so, first compile snpknock2 
# by entering the 'snpknock2' directory and running 'make'.

# Path to snpknock2 executable built as described above
SNPKNOCK2="../snpknock2/bin/snpknock2"

# Create directory for output
mkdir -p "../tmp/knockoffs"

# Run snpknock2
$SNPKNOCK2 \
  --bgen "../data/haplotypes/example_chr{21:22}" \
  --keep "../data/qc/samples_qc.txt" \
  --extract "../data/qc/qc_chr{21:22}.txt" \
  --map "../data/maps/genetic_map_chr{21:22}.txt" \
  --part "../data/partitions/example_chr{21:22}.txt" \
  --ibd "../data/ibd/example_chr{21:22}.txt" \
  --K 10 \
  --cluster_size_min 1000 \
  --cluster_size_max 10000 \
  --hmm-rho 1 \
  --hmm-lambda 1e-3 \
  --windows 0 \
  --n_threads 1 \
  --seed 2020 \
  --compute-references \
  --generate-knockoffs \
  --out "../tmp/knockoffs/knockoffs_chr{21:22}"

