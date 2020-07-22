#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Visualize KnockoffZoom discoveries on a toy dataset
#
# Authors: Matteo Sesia
# Date:    04/24/2019

#########
# Setup #
#########

# Print header
printf "KnockoffZoomView v0.1 (24 Apr 2019) \n"
printf "https://bitbucket.org/msesia/knockoffzoom \n"
printf "(C) 2019 Matteo Sesia, Eugene Katsevich (Stanford University)   GNU General Public License v3 \n\n"

# Setup spinner for long jobs
source "misc/spinner.sh"

# Log file
LOG_FILE="knockoffzoomview.log"
rm -f $LOG_FILE
touch $LOG_FILE
echo "Log file: "$LOG_FILE

######################
# Check dependencies #
######################

printf "\nSetup\n"
# System dependencies
check_dependency () {
  CMD=$1
  if [ ! -x "$(command -v $CMD)" ]; then
    echo "Error: command $CMD not available"
    exit
  fi
}
DEPENDENCY_LIST=("R")
start_spinner " - Checking system dependencies..."
for DEPENDENCY in "${DEPENDENCY_LIST[@]}"; do    
  check_dependency $DEPENDENCY &>> "../"$LOG_FILE
done
stop_spinner $?

# R libraries
start_spinner " - Checking R library dependencies..."
Rscript --vanilla "visualization/check_packages.R" &>> "../"$LOG_FILE
stop_spinner $?

# Variant annotations
start_spinner " - Checking variant annotations..."
cd misc
./download_annotations.sh
cd ..
stop_spinner $?

#########################
# Visualize discoveries #
#########################

# Enter the directory where the scripts are stored
cd visualization

printf "\nStarting Shiny app...\n"
R -e "shiny::runApp('app.R', launch.browser=T)" &>> "../"$LOG_FILE
printf "Exiting.\n"
