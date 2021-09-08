#!/bin/bash

FA_TO_SAMPLE="../../core.full.aln"
N_SAMPLE=100
SIZE_SAMPLE=100
MIN_SNPS=10

for i in $(seq 6 "$N_SAMPLE")
do
	# Randomly sample fasta sequence in core.full.aln
	cat $FA_TO_SAMPLE | awk '/^>/ { if(i>0) printf("\n"); i++; printf("%s\t",$0); next;} {printf("%s",$0);} END { printf("\n");}' | shuf | head -n $SIZE_SAMPLE | awk '{printf("%s\n%s\n",$1,$2)}' > core.full.aln.samp${SIZE_SAMPLE}_${i}
	# Run gubbins on subsampled fasta
	run_gubbins.py --min_snps $MIN_SNPS --prefix gubbins.samp${SIZE_SAMPLE}_${i} --threads 70 core.full.aln.samp${SIZE_SAMPLE}_${i} 
done


