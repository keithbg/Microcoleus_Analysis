#!/usr/bin bash

## Remove the gene locus from the fasta headers from Roary output
## This is necessary for input into PopGenome, so that genes within the same genome can be identified

## Roary header: >PH2015_01U_Oscillatoriales_45_14_00178
## Output header: >PH2015_01U_Oscillatoriales_45_14


for file in $(ls -1 *.aln)
  do
  sed -i .bck 's/_[^_]*$//' "$file"
done

rm *.bck
