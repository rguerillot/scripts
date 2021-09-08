#!/usr/bin/env bash

#Author: R. Guerillot
#Description: Use trimal -sgc option to calculate %gap, %ambiguous and %snp at every position of core.full.aln from snippy 

#$aln = core.full.aln

aln=$1

echo "processing $aln file"

echo "get %gap at each position"
trimal -in $aln -sgc > ${aln}.pct_gap

echo "get %ambiguous or gap at each position"
(sed 's|N|-|g' | sed 's|n|-|g') < ${aln} > ${aln}_ambiguous_and_gap_to-
trimal -in ${aln}_ambiguous_and_gap_to- -sgc > ${aln}.pct_gap_ambiguous

echo "get %ambiguous at each position"
sed 's|-|T|g' < ${aln} > ${aln}_gap_to_T
(sed 's|N|-|g' | sed 's|n|-|g') < ${aln}_gap_to_T > ${aln}_ambiguous_to-
trimal -in ${aln}_ambiguous_to- -sgc > ${aln}.pct_ambiguous

echo "get %snp at each position"
sed 's|n|T|g' < ${aln}_gap_to_T > ${aln}_gap_and_n_to_T
tr '[:lower:]' '-' < ${aln}_gap_and_n_to_T > ${aln}_snp_to-
trimal -in ${aln}_snp_to- -sgc > ${aln}.pct_snp

echo "generate artemis userplots"
cut -f3 ${aln}.pct_gap | grep -v "|" | grep -v "+" | awk '{printf "%i\n", $1}' > ${aln}.pct_gap.userplot
cut -f3 ${aln}.pct_gap_ambiguous | grep -v "|" | grep -v "+" | awk '{printf "%i\n", $1}' > ${aln}.pct_gap_ambiguous.userplot
cut -f3 ${aln}.pct_ambiguous | grep -v "|" | grep -v "+" | awk '{printf "%i\n", $1}' > ${aln}.pct_ambiguous.userplot
cut -f3 ${aln}.pct_snp | grep -v "|" | grep -v "+" | awk '{printf "%i\n", $1}' > ${aln}.pct_snp.userplot

echo "clean intermediate files"
rm ${aln}_ambiguous_and_gap_to- ${aln}_gap_to_T ${aln}_ambiguous_to- ${aln}_gap_and_n_to_T ${aln}_snp_to-
