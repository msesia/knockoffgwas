#!/bin/bash

# If you haven't already done so, compile snpknock2.
# To compile snpknock2, enter the 'snpknock2' directory and run 'make'.

# Path to snpknock2 executable built as described above
SNPKNOCK2="../snpknock2/bin/snpknock2"

./$SNPKNOCK2 \
  --bgen ../data/haplotypes/example_chr{21:22} \
