#!/bin/csh

#########################################################################
# File Name: GO_MMseqs2.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Protein similarity search against RefSeq viral proteins using MMseqs2. (Steinegger and Söding, 2017)"
        echo ""
        echo "Requires: MMseqs2 and protein FASTA files"
        echo ""
        echo "Usage: GO_MMseqs2.sh <query_faa> <reference_faa> <output_prefix> <temporary_directory> "
        echo ""
    exit 0
endif

set query_faa = "$1"
set reference_faa = "$2"
set output_prefix = "$3"
set temporary_directory = "$4"

#########################################################################


echo "Running MMseqs2"

module load mmseqs2/17.gitb804f

mmseqs createdb $reference_faa "$output_prefix"_reference_db
mmseqs createdb $query_faa "$output_prefix"_query_db

mmseqs search \
    "$output_prefix"_query_db \
    "$output_prefix"_reference_db \
    "$output_prefix"_search \
    $temporary_directory \
    -e 1e-5 \

mmseqs convertalis \
    "$output_prefix"_query_db \
    "$output_prefix"_reference_db \
    "$output_prefix"_search \
    "$output_prefix".m8 \
    --format-output "query,target,evalue,bits,qcov,tcov,pident,alnlen"

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
