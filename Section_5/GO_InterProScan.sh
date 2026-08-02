#!/bin/csh

#########################################################################
# File Name: GO_InterProScan.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Functional annotation of protein sequences using InterProScan. (Jones et al., 2014)"
        echo ""
        echo "Requires: InterProScan v5.75-106.0 and *.faa files"
        echo ""
        echo "Usage: GO_InterProScan.sh <path_fasta_files> <path_output_files> "
        echo ""
    exit 0
endif

set path_fasta_files = "$1"
set path_output_files = "$2"

#########################################################################


echo "Running InterProScan"

foreach f ("$path_fasta_files"/*.faa)
    set fasta_filename = `basename $f`
    set output_file = "$path_output_files/$fasta_filename:r.tsv"

    interproscan -i $f -f tsv -o $output_file -dp
end

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
