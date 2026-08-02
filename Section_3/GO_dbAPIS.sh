#!/bin/csh

#########################################################################
# File Name: GO_dbAPIS.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Detection of anti-defense proteins using the dbAPIS workflow. (Yan et al., 2024)"
        echo ""
        echo "Requires: HMMER v3.4, DIAMOND v2.0.15, dbAPIS.hmm, APIS_db and parse_annotation_result.sh"
        echo ""
        echo "dbAPIS database files can be downloaded from: https://github.com/azureycy/dbAPIS"
        echo ""
        echo "Usage: GO_dbAPIS.sh <path_fasta_files> <path_output_files> <dbAPIS_hmm> <APIS_database> <parse_script>"
        echo ""
        exit 0
endif

set path_fasta_files = "$1"
set path_output_files = "$2"
set dbAPIS_hmm = "$3"
set APIS_database = "$4"
set parse_script = "$5"

#########################################################################


echo "Running dbAPIS"

module load hmmer/3.4
module load diamond/2.0.15

foreach f ("$path_fasta_files"/*.faa)
    set fasta_filename = `basename $f`
    set prefix = "$fasta_filename:r"
    set hmm_output = "$path_output_files/${prefix}_hmm.out"
    set diamond_output = "$path_output_files/${prefix}_diamond.out"

    hmmscan --domtblout $hmm_output --noali $dbAPIS_hmm $f

    diamond blastp \
        --db $APIS_database \
        -q $f \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
        -o $diamond_output

    bash $parse_script $hmm_output $diamond_output
end

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
