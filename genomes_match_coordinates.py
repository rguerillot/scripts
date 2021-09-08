#!/usr/bin/env python3

# author: R. Guerillot

# Parse mauve whole genome alignment and output a table with all matching coordinates and SNP/INS/INDEL annotations
## run mauve Aligner first:
### progressiveMauve --output=NRS384_N315.xmfa NRS384_w_annot_FPR3757.fasta N315.fasta

from Bio import AlignIO

alignments= AlignIO.parse(open("./mauveAligner_results/NRS384_N315.xmfa"), "mauve")

for alignment in alignments:
    if len(alignment) == 2:
        recordA = alignment[0]
        recordA_name = recordA.id.split("/")[0]
        recordA_start = recordA.annotations["start"] + 1
        recordA_end = recordA.annotations["end"] + 1
        recordA_seq = recordA.seq
        recordB = alignment[1]
        recordB_name = recordB.id.split("/")[0]
        recordB_start = recordB.annotations["start"] + 1
        recordB_end = recordB.annotations["end"] + 1
        recordB_seq = recordB.seq
        #print(recordB_name + " " + str(recordB_start) + " " + str(recordB_end)+ "\n" + recordB_seq[0:10])
        


print("\t".join(["NRS384_coor", "N315_coor","NRS384_base","N315_base","mutation"]))

for i in range(len(recordA_seq)):
    if recordA_seq[i] == recordB_seq[i]:
        mutation = ""
    if recordA_seq[i] != recordB_seq[i]:
        mutation = "SNP"
    if recordB_seq[i] == "-":
        mutation = "DEL"
    if recordA_seq[i] == "-":
        mutation = "INS"
    if recordA_seq[i] != "-":
        print("\t".join([str(recordA_start), str(recordB_start), recordA_seq[i],recordB_seq[i], mutation]))
        recordA_start +=1
    if recordB_seq[i] != "-":
        recordB_start += 1
    
            
        
    
