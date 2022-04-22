#!/bin/bash

FASTA=$1
#OUTPUT=$2

# NOTE: " " (space), ";" (semicolon), ":" (colon), "," (comma), "()" (parentheses), "'" (quote). 

echo "Cleaning: "$FASTA

sed -i 's/ /_/g' $FASTA 
sed -i 's/;/_/g' $FASTA 
sed -i 's/:/_/g' $FASTA 
sed -i 's/,/_/g' $FASTA 
#sed -i 's/(/_/g' $FASTA 
#sed -i 's/)/_/g' $FASTA 
#sed -i 's/'/_/g' $FASTA 

exit 0
