#!/usr/bin/env python3
"""gbktofasta.py: Convert one or multiple genbank file to fata format"""

__author__      = "Romain Guerillot"

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

print('genbank conversion to fasta\nusage gbktofasta <gbk file(s)>')

for gb in sys.argv[1:]:
    print("processing: " + gb)
    my_gbk_name = gb.replace(".gbk", "").replace(".genbank", "").replace(".gbff", "").replace(".gb", "")
    SeqIO.convert(gb, "genbank", my_gbk_name + ".fa", "fasta")

