#!/bin/bash
#
#

# Range of chromosomes to include in the analysis
CHR_MIN=21
CHR_MAX=22
CHR_LIST=$(seq $CHR_MIN $CHR_MAX)

GENO_FILE="../data/genotypes/example_chr"
FAM_FILE="../data/genotypes/example_chr"$CHR_MIN".fam"

PHENO_FILE="../data/phenotypes/phenotypes.tab"
PHENO_NAME="y"

OUT_DIR="../data/lmm"
mkdir -p $OUT_DIR
STATS_FILE=$OUT_DIR"/example_lmm.txt"
CLUMP_BASENAME=$OUT_DIR"/example"

CLUMP_THRESHOLD=0.00000005 # 5e-8

####################
# Compute p-values #
####################

# Stuff for bolt
LD_TABLE="../aux/LDSCORE.1000G_EUR.tab.gz"
MAP_TABLE="../aux/genetic_map_hg17_withX.txt.gz"

bolt \
    --bed=$GENO_FILE"{$CHR_MIN:$CHR_MAX}.bed" \
    --bim=$GENO_FILE"{$CHR_MIN:$CHR_MAX}.bim" \
    --fam=$FAM_FILE \
    --maxMissingPerSnp=1 \
    --phenoFile=$PHENO_FILE \
    --phenoCol=$PHENO_NAME \
    --covarFile=$PHENO_FILE \
    --covarCol="sex" \
    --LDscoresFile=$LD_TABLE \
    --geneticMapFile=$MAP_TABLE \
    --statsFile=$STATS_FILE \
    --lmm \
    --numThreads=1

echo "Output file:"
echo $STATS_FILE

#####################
# Clump nearby loci #
#####################

# Create temporary stats file with modified header
STATS_FILE_TMP=$STATS_FILE".tmp"
head -n 1 $STATS_FILE > $STATS_FILE_TMP
sed -i 's/\S\+$/P/' $STATS_FILE_TMP
tail -n +2 $STATS_FILE >> $STATS_FILE_TMP

for CHR in $CHR_LIST; do
  # Clumping of LMM p-values
  CLUMP_FILE=$CLUMP_BASENAME"_lmm_clumped"
  echo "Clumping chromosome $CHR at threshold $CLUMP_THRESHOLD"
  plink --bfile $GENO_FILE$CHR \
        --clump $STATS_FILE_TMP \
        --clump-p1 $CLUMP_THRESHOLD \
        --clump-r2 0.01 \
        --clump-kb 5000 \
        --out $CLUMP_FILE"_chr"$CHR
  rm -f $CLUMP_FILE"_chr"$CHR".log"
  rm -f $CLUMP_FILE"_chr"$CHR".nosex"
done

# Erase temporary stats file with modified header
rm $STATS_FILE_TMP

# Combine results of clumping into a single file
echo "Combining clumping files at threshold $CLUMP_THRESHOLD"
CLUMP_FILE=$CLUMP_BASENAME"_lmm_clumped"
# Write header
HEADER=" CHR    F             SNP         BP        P    TOTAL   NSIG    S05    S01   S001  S0001    SP2"
echo $HEADER > $CLUMP_FILE".tab"
for CHR in $CHR_LIST; do
  if [ -s $CLUMP_FILE"_chr"$CHR".clumped" ]; then
    tail -n +2 $CLUMP_FILE"_chr"$CHR".clumped" >> $CLUMP_FILE".tab"
    rm $CLUMP_FILE"_chr"$CHR".clumped"
  fi
done
# Remove empty lines
sed -i '/^$/d' $CLUMP_FILE".tab"
echo "Results written on "$CLUMP_FILE".tab"

# Parse clumped p-values and summarise discoveries
OUT_FILE=$CLUMP_BASENAME"_lmm_regions.txt"
Rscript --vanilla summarise_lmm.R
