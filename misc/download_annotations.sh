#!/bin/bash
#
#

URL_REMOTE="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/"
FUNC_FILE="wgEncodeBroadHmmGm12878HMM.txt"
GENE_FILE="ncbiRefSeq.txt"

mkdir -p "../data"
mkdir -p "../data/annotations"

cd "../data/annotations"

# Download functional annotation file
if [ ! -f $FUNC_FILE ]; then
  rm -f $FUNC_FILE".gz"
  wget $URL_REMOTE$FUNC_FILE".gz"
  gunzip $FUNC_FILE".gz"
fi

# Download gene list
if [ ! -f $GENE_FILE ]; then
  rm -f $GENE_FILE".gz"
  wget $URL_REMOTE$GENE_FILE".gz"
  gunzip $GENE_FILE".gz"
fi
