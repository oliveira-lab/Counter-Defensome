#!/bin/csh

#########################################################################
# File Name: GO_DIAMOND.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Protein similarity search using DIAMOND blastp. (Buchfink et al., 2021)"
        echo ""
        echo "Requires: DIAMOND v2.0.15, a DIAMOND database and *.faa files"
        echo ""
        echo "Usage: GO_DIAMOND.sh <path_fasta_files> <path_output_files> <diamond_database>"
        echo ""
    exit 0
endif

set path_fasta_files = "$1"
set path_output_files = "$2"
set diamond_database = "$3"

#########################################################################


echo "Running DIAMOND"

module load diamond/2.0.15

foreach f ("$path_fasta_files"/*.faa)
    set fasta_filename = `basename $f`
    set output_file = "$path_output_files/$fasta_filename:r.diamond.tsv"

    diamond blastp \
        --db $diamond_database \
        -q $f \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
        -o $output_file
end

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
