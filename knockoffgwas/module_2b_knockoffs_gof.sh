#!/bin/bash
#
# Evaluate knockoff goodness of fit
#
# Authors: Matteo Sesia
# Date:    07/27/2021

set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

# Temporary storage of intermediate files
TMP_DIR="../tmp"
mkdir -p $TMP_DIR
mkdir -p $TMP_DIR"/knockoffs"
mkdir -p $TMP_DIR"/knockoffs_gof"
mkdir -p $TMP_DIR"/knockoffs_gof_plots"

# Storage of output files
OUT_DIR="../results"
mkdir -p $OUT_DIR

# List of chromosomes
CHR_MIN=21
CHR_MAX=22
CHR_LIST=($(seq $CHR_MIN $CHR_MAX))

# List of resolutions
RESOLUTION_LIST=("6" "5" "4" "3" "2" "1" "0")

# Utility scripts
PLOT_GOF="Rscript --vanilla utils/knockoffs_gof.R"

# Which operations should we perform?
FLAG_GOF_KNOCKOFFS=1

##########################
# Evaluate Knockoffs GOF #
##########################

if [[ $FLAG_GOF_KNOCKOFFS == 1 ]]; then
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Checking knockoffs goodness-of-fit"
  echo "----------------------------------------------------------------------------------------------------"

  for RESOLUTION in "${RESOLUTION_LIST[@]}"; do
    for CHR in "${CHR_LIST[@]}"; do
    echo ""
    echo "Computing goodness-of-fit diagnostics for chromosome "$CHR" at resolution "$RESOLUTION" ..."
    echo ""

    STATS_BASENAME=$TMP_DIR"/knockoffs_gof/knockoffs_chr"$CHR"_res"$RESOLUTION
    GROUPS_FILE=$TMP_DIR"/knockoffs/knockoffs_chr"$CHR"_res"$RESOLUTION"_grp.txt"
    OUT_BASENAME=$TMP_DIR"/knockoffs_gof_plots/knockoffs_chr"$CHR"_res"$RESOLUTION

    # Compute self-similarity diagnostics
    plink --bfile $TMP_DIR"/knockoffs/knockoffs_chr"$CHR"_res"$RESOLUTION \
          --keep-allele-order \
          --freq \
          --r in-phase --ld-window 2 --ld-window-kb 0 \
          --memory 9000 \
          --out $STATS_BASENAME"_self"

    # Compute exchangeability diagnostics
    plink --bfile $TMP_DIR"/knockoffs/knockoffs_chr"$CHR"_res"$RESOLUTION \
          --keep-allele-order \
          --r2 in-phase --ld-window 5000 --ld-window-kb 10000 --ld-window-r2 0.01 \
          --memory 9000 \
          --out $STATS_BASENAME

    # Make GOF plots
    $PLOT_GOF $CHR $RESOLUTION $STATS_BASENAME $GROUPS_FILE $OUT_BASENAME
    done
  done

else
  echo ""
  echo "----------------------------------------------------------------------------------------------------"
  echo "Skipping knockoffs goodness-of-fit"
  echo "----------------------------------------------------------------------------------------------------"
fi
