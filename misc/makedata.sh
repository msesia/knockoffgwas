#!/bin/bash
#
#

#####################
# Make genetic data #
#####################
Rscript --vanilla make_genetic.R

CHR_LIST=$(seq 21 22)

######################
# Process haplotypes #
######################
INP_DIR="/scratch/groups/candes/ukbiobank_public/data/tmp"
HMM_DIR="/scratch/groups/candes/ukbiobank_public/data/hmm"
OUT_DIR="/scratch/groups/candes/ukbiobank_public/data/haplotypes"

mkdir -p $HMM_DIR
mkdir -p $OUT_DIR

for CHR in $CHR_LIST; do

  BASENAME=$INP_DIR"/example_chr"$CHR
  BASENAME_OUT=$OUT_DIR"/example_chr"$CHR

  plink2 \
    --haps $BASENAME".haps" \
    --legend $BASENAME".leg" $CHR \
    --sample $BASENAME".sample" \
    --export bgen-1.3 \
    --out $BASENAME_OUT

  rm $BASENAME_OUT".log"

done

###################
# Quality control #
###################
GEN_DIR="/scratch/groups/candes/ukbiobank_public/data/genotypes"
QC_DIR="/scratch/groups/candes/ukbiobank_public/data/qc"
mkdir -p $QC_DIR

# List of individuals that passed QC
awk '{ print $1,$2 }' $GEN_DIR"/example_chr21.fam" > $QC_DIR"/samples_qc.txt"

# List of variants that passed QC
rm -f $QC_DIR"/variants_qc.txt"
touch $QC_DIR"/variants_qc.txt"
for CHR in $CHR_LIST; do
  awk '{ print $2 }' $GEN_DIR"/example_chr"$CHR".bim" >> $QC_DIR"/variants_qc.txt"
done

########
# Copy #
########
cp -r "/scratch/groups/candes/ukbiobank_public/data/genotypes" "../data/"
cp -r "/scratch/groups/candes/ukbiobank_public/data/haplotypes" "../data/"
cp -r "/scratch/groups/candes/ukbiobank_public/data/qc" "../data/"
cp -r "/scratch/groups/candes/ukbiobank_public/data/hmm" "../data/"

###################
# Make phenotypes #
###################
Rscript --vanilla make_phenotypes.R
