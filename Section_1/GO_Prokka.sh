#!/bin/csh

#########################################################################
# File Name: GO_Prokka.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Annotation of viral genomes using Prokka. (Seemann, 2014)"
        echo ""
        echo "Requires: Prokka v1.14.5 and *.fna files"
        echo ""
        echo "Usage: GO_Prokka.sh <path_fasta_files> <path_output_files> "
        echo ""
    exit 0
endif

set path_fasta_files = "$1"
set path_output_files = "$2"

#########################################################################


echo "Running Prokka"

foreach f ("$path_fasta_files"/*.fna)
    set fasta_filename = `basename $f`
    set prefix = "$fasta_filename:r"
    set output_dir = "$path_output_files/$prefix"

    prokka $f --outdir $output_dir --prefix $prefix
end

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
