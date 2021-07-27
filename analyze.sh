#!/bin/bash
# UK Biobank GWAS
#
# Class: script
#
# Run KnockoffGWAS on a toy dataset
#
# Authors: Matteo Sesia
# Date:    04/24/2019

#########
# Setup #
#########

# Print header
printf "KnockoffGWAS (21 July 2020) \n"
printf "https://bitbucket.org/msesia/knockoffgwas \n"
printf "(C) 2020 Matteo Sesia   GNU General Public License v3 \n\n"

# Setup spinner for long jobs
source "misc/spinner.sh"

# Log file
LOG_FILE="knockoffgwas.log"
rm -f $LOG_FILE
touch $LOG_FILE
echo "Log file: "$LOG_FILE

# Enter the directory where the scripts are stored
cd knockoffgwas

######################
# Check dependencies #
######################

printf "\nSetup\n"
# System dependencies
ERROR=0
check_dependency () {
  CMD=$1
  if [ ! -x "$(command -v $CMD)" ]; then
    echo -e "Error: command $CMD not available"
    ERROR=1
  fi
}
DEPENDENCY_LIST=("plink" "R" "../snpknock2/bin/snpknock2")
start_spinner " - Checking system dependencies..."
for DEPENDENCY in "${DEPENDENCY_LIST[@]}"; do  
  check_dependency $DEPENDENCY &>> "../"$LOG_FILE
done
stop_spinner $ERROR

# R libraries
start_spinner " - Checking R library dependencies..."
Rscript --vanilla "utils/check_packages.R" &>> "../"$LOG_FILE
stop_spinner $?

####################
# Run KnockoffGWAS #
####################

printf "\nData analysis\n"

# Module 1: partition the genome into LD blocks
start_spinner ' - Running module 1 (partitioning the genome)...'
./module_1_partition.sh &>> "../"$LOG_FILE
stop_spinner $?

# Module 2: generate the knockoffs
start_spinner ' - Running module 2 (generating knockoffs)...'
./module_2_knockoffs.sh &>> "../"$LOG_FILE
stop_spinner $?

# Module 3: compute the test statistics
start_spinner ' - Running module 3 (computing test statistics)...'
./module_3_statistics.sh &>> "../"$LOG_FILE
stop_spinner $?

# Module 4: report significant findings
start_spinner ' - Running module 4 (applying the knockoff filter)...'
./module_4_discover.sh &>> "../"$LOG_FILE
stop_spinner $?

#####################
# Summarize results #
#####################

printf "\nResults written in 'results/'\n"
