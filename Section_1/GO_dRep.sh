#!/bin/csh

#########################################################################
# File Name: GO_dRep.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Dereplication of viral genomes using dRep. (Olm et al., 2017)"
        echo ""
        echo "Requires: dRep v3.6.2 and *.fna files"
        echo ""
        echo "Usage: GO_dRep.sh <path_fasta_files> <path_output_files> <minimum_length>"
        echo ""
    exit 0
endif

set path_fasta_files = "$1"
set path_output_files = "$2"
set minimum_length = "$3"

#########################################################################


echo "Running dRep"

dRep dereplicate $path_output_files \
    -g "$path_fasta_files"/*.fna \
    -pa 1 \
    -l $minimum_length \
    --skip_plots \
    --ignoreGenomeQuality

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
