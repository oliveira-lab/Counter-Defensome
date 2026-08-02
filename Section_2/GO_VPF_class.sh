#!/bin/csh

#########################################################################
# File Name: GO_VPF_Class.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Taxonomic classification of viral sequences using VPF-class. (Pons et al., 2021)"
        echo ""
        echo "Requires: VPF-class, VPF-class data index and a nucleotide FASTA file"
        echo ""
        echo "Usage: GO_VPF_Class.sh <fasta_file> <path_output_files> <data_index>"
        echo ""
    exit 0
endif

set fasta_file = "$1"
set path_output_files = "$2"
set data_index = "$3"

#########################################################################


echo "Running VPF-class"

stack exec -- vpf-class \
    --data-index $data_index \
    -i $fasta_file \
    -o $path_output_files

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
