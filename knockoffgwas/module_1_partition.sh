#!/bin/bash
#
# Partition the variants through adjacent clustering
#
# Authors: Matteo Sesia
# Date:    07/21/2020

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

# Temporary storage of intermediate files
TMP_DIR="../tmp"
mkdir -p $TMP_DIR
mkdir -p $TMP_DIR"/partitions"

# List of chromosomes
CHR_LIST=$(seq 21 22)

# List of resolutions in cM (note: this is also fixed inside the R script)
RESOLUTION_LIST=("0" "0.01" "0.05" "0.1" "0.2" "0.5" "1")

# Utility scripts
PARTITION_VARIANTS="Rscript --vanilla utils/partition.R"

# Which operations should we perform?
FLAG_PARTITION=1

################
# Partitioning #
################

if [[ $FLAG_PARTITION == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Partitioning variants"
  echo "----------------------------------------------------------------------------------------------------"

  for CHR in $CHR_LIST; do

    echo ""
    echo "Processing chromosome "$CHR" ..."
    echo ""

    # Input genotype files (PLINK format)
    GENO_BIM="../data/genotypes/example_chr"$CHR".bim"

    # List of variants that passed QC
    QC_VARIANTS="../data/qc/variants_qc.txt"

    # Genetic map file
    GEN_MAP="../data/maps/genetic_map_chr"$CHR".txt"

    # Basename for output dendrogram file
    OUT_FILE=$TMP_DIR"/partitions/example_chr"$CHR".txt"

    # Partition the variants at different resolutions
    $PARTITION_VARIANTS $GEN_MAP $GENO_BIM $QC_VARIANTS $OUT_FILE
  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping variant partitioning"
  echo "----------------------------------------------------------------------------------------------------"
fi
