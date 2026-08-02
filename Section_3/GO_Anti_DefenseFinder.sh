#!/bin/csh

#########################################################################
# File Name: GO_AntiDefenseFinder.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Detection of anti-defense genes / systems using DefenseFinder anti-defense mode. (Tesson et al., 2025)"
        echo ""
        echo "Requires: DefenseFinder v2.0.0 and *.faa files"
        echo ""
        echo "Usage: GO_AntiDefenseFinder.sh <path_fasta_files> <path_output_files>"
        echo ""
    exit 0
endif

set path_fasta_files = "$1"
set path_output_files = "$2"

#########################################################################


echo "Running DefenseFinder anti-defense mode"

foreach f ("$path_fasta_files"/*.faa)
    set fasta_filename = `basename $f`
    set output_dir = "$path_output_files/$fasta_filename:r"

    defense-finder run $f --antidefensefinder-only --preserve-raw -o $output_dir
end

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
