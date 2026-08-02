#!/bin/csh

#########################################################################
# File Name: GO_CheckV.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Assessment of viral genome quality using CheckV. (Nayfach et al., 2021)"
        echo ""
        echo "Requires: CheckV v1.0.3, CheckV database v1.5 and *.fna files"
        echo ""
        echo "Usage: GO_CheckV.sh <path_fasta_files> <path_output_files> <path_database>"
        echo ""
    exit 0
endif

set path_fasta_files = "$1"
set path_output_files = "$2"
set path_database = "$3"

#########################################################################


echo "Running CheckV"

setenv CHECKVDB "$path_database"

foreach f ("$path_fasta_files"/*.fna)
    set fasta_filename = `basename $f`
    set output_dir = "$path_output_files/$fasta_filename:r"

    checkv end_to_end $f $output_dir
end

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
