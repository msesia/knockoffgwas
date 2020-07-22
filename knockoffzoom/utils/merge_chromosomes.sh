#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Merge knockoffs for different chromosomes
#
# Authors: Matteo Sesia
# Date:    07/21/2020

OUT_BASENAME=$1
RESOLUTION=$2
CHR_MIN=$3
CHR_MAX=$4

# List of chromosomes
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

# Make list of chromosomes to be merged with the first one
MERGE_LIST=$OUT_BASENAME"_mergelist.txt"
rm -f $MERGE_LIST
touch $MERGE_LIST
for CHR in ${CHR_LIST[@]}; do
  if [ $CHR != $CHR_MIN ]; then
    # Basename for the augmented genotype file for this chromosome
    CHR_BASENAME="../tmp/knockoffs/knockoffs_chr"$CHR"_res"$RESOLUTION
    echo $CHR_BASENAME >> $MERGE_LIST
  fi
done

# Merge the augmented data from all chromosomes
CHR_BASENAME_FIRST="../tmp/knockoffs/knockoffs_chr"$CHR_MIN"_res"$RESOLUTION
plink \
  --bfile $CHR_BASENAME_FIRST \
  --merge-list $MERGE_LIST \
  --make-bed \
  --memory 5000 \
  --out $OUT_BASENAME

rm -f $OUT_BASENAME".log" $MERGE_LIST

# Make list of SNP groups
GRP_LIST=$OUT_BASENAME"_grp.txt"
rm -f $GRP_LIST
touch $GRP_LIST
for CHR in ${CHR_LIST[@]}; do
  # Basename for the augmented genotype file for this chromosome
  GRP_FILE="../tmp/knockoffs/knockoffs_chr"$CHR"_res"$RESOLUTION"_grp.txt"
  tail -n +2 $GRP_FILE >> $GRP_LIST
done
