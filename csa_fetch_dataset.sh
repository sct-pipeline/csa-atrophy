#!/bin/bash
#
#  fetching data from openneuro if data does not already exist
#
# Usage:
#   ./CSA-dataset
#
# Author: Paul Bautin


# Uncomment for full verbose
#set -v

# Immediately exit if error
set -e

# Exit if user presses CTRL+C (Linux) or CMD+C (OSX)
trap "echo Caught Keyboard Interrupt within script. Exiting now.; exit" INT


# FUNCTIONS
# ==============================================================================
# input data to be processed
SUBJECTS_DIR=(
  "sub-amu01"
  "sub-amu02"
)



for SUBJECT in ${SUBJECTS_DIR[@]}; do
  if [ ! -e $SUBJECT ]; then
    echo "check if $SUBJECT exists if not downloading from openneuro.org"
    aws s3 sync --no-sign-request s3://openneuro.org/ds001919/$SUBJECT data/$SUBJECT
  fi
done
