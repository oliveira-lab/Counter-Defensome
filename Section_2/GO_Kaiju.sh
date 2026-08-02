#!/bin/csh

#########################################################################
# File Name: GO_Kaiju.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 31/07/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Taxonomic classification of viral contigs using Kaiju. (Menzel et al., 2016)"
        echo ""
        echo "Requires:"
        echo "  - Kaiju v1.7.3"
        echo "  - Nucleotide FASTA files (*.fna)"
        echo "  - A Kaiju protein database index (*.fmi)"
        echo "  - nodes.dmp: NCBI taxonomic tree (taxon IDs, parents and ranks) (from kaiju)"
        echo "  - names.dmp: NCBI taxon names associated with taxon IDs (from kaiju)"
        echo ""
        echo "Usage:"
        echo "  GO_Kaiju.sh <path_fasta_files> <path_output_files> <nodes.dmp> <names.dmp> <kaiju_database.fmi>"
        echo ""
        exit 0
endif

set path_fasta_files = "$1"
set path_output_files = "$2"
set nodes_file = "$3"
set names_file = "$4"
set database_file = "$5"

#########################################################################


echo "Running Kaiju"

module load kaiju/1.7.3

foreach f ("$path_fasta_files"/*.fna)
    set fasta_filename = `basename $f`
    set prefix = "$fasta_filename:r"
    set kaiju_output = "$path_output_files/$prefix.out"
    set table_output = "$path_output_files/$prefix.out.tsv"

    kaiju -t $nodes_file -f $database_file -i $f -o $kaiju_output
    kaiju2table -t $nodes_file -n $names_file -r species -p -e -o $table_output $kaiju_output
end

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################
