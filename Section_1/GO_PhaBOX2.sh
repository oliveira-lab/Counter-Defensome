#!/bin/csh

#########################################################################
# File Name: GO_PhaBOX2.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Identification and characterization of viral contigs using PhaBOX2. (Shang et al., 2023)"
        echo ""
        echo "Requires: PhaBOX2 v2.1.12, PhaBOX2 database and *.fna files"
        echo ""
        echo "Usage: GO_PhaBOX2.sh <path_fasta_files> <path_output_files> <path_database>"
        echo ""
    exit 0
endif

set path_fasta_files = "$1"
set path_output_files = "$2"
set path_database = "$3"

#########################################################################


echo "Running PhaBOX2"

foreach f ("$path_fasta_files"/*.fna)
    set fasta_filename = `basename $f`
    set output_dir = "$path_output_files/$fasta_filename:r"

    phabox2 --dbdir $path_database --contigs $f --outpth $output_dir 
end

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
