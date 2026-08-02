#!/bin/csh

#########################################################################
# File Name: GO_iPHoP.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Prediction of phage hosts using iPHoP. (Roux et al., 2023)"
        echo ""
        echo "Requires: iPHoP v1.4.1, iPHoP database and *.fna files"
        echo ""
        echo "Usage: GO_iPHoP.sh <path_fasta_files> <path_output_files> <path_database>"
        echo ""
    exit 0
endif

set path_fasta_files = "$1"
set path_output_files = "$2"
set path_database = "$3"

#########################################################################


echo "Running iPHoP"

setenv IPHOP_DB_DIR "$path_database"

foreach f ("$path_fasta_files"/*.fna)
    set fasta_filename = `basename $f`
    set output_dir = "$path_output_files/$fasta_filename:r"

    iphop predict --fa_file $f -d $IPHOP_DB_DIR --out_dir $output_dir
end

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
