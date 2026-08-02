#!/bin/csh

#########################################################################
# File Name: GO_vConTACT3.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Gene-sharing network analysis using vConTACT3. (Bolduc et al., 2025)"
        echo ""
        echo "Requires: vConTACT3 v3.1.6, vConTACT3 database and a nucleotide FASTA file"
        echo ""
        echo "Usage: GO_vConTACT3.sh <fasta_file> <path_output_files> <path_database> "
        echo ""
    exit 0
endif

set fasta_file = "$1"
set path_output_files = "$2"
set path_database = "$3"

#########################################################################


echo "Running vConTACT3"

vcontact3 run \
    --nucleotide $fasta_file \
    --output $path_output_files \
    --db-path $path_database \
    --db-version 230 \

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
