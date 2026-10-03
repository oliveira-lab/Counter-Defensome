#!/bin/csh

#########################################################################
# File Name: GO_SegB_raw_files.sh
# Author(s): Angelina Beavogui
# Institution: Genoscope, Evry, France
# Mail: beavogui67@gmail.com
# Date: 25/09/2026
#########################################################################


if ("$1" == "-h" || "$1" == "-help") then
        echo ""
        echo "Search for SegB homologs and tRNAs in one viral genome."
        echo ""
        echo "Requires: HMMER, tRNAscan-SE, a SegB HMM profile, a protein FASTA file and a nucleotide FASTA file"
        echo ""
        echo "Usage: GO_SegB_raw_files.sh <SegB_profile.hmm> <genome.faa> <genome.fna> <path_output_files>"
        echo ""
    exit 0
endif

set hmm_profile = "$1"
set protein_file = "$2"
set genome_file = "$3"
set path_output_files = "$4"
set genome_name = `basename "$protein_file"`
set genome_name = "$genome_name:r"

#########################################################################


mkdir -p "$path_output_files"

echo "Running hmmsearch"

hmmsearch -T 20 --domT 20 --domtblout "$path_output_files/$genome_name.domtblout" "$hmm_profile" "$protein_file" > "$path_output_files/$genome_name.hmmsearch.txt"

echo "Running tRNAscan-SE"

tRNAscan-SE -o "$path_output_files/$genome_name.tRNA.out" -j "$path_output_files/$genome_name.tRNA.json" "$genome_file" >& "$path_output_files/$genome_name.tRNAscan.log"

tput setaf 2; echo "Done!"; tput sgr0

#########################################################################