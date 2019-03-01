#!/usr/bin/env python
from plumbum import local
from plumbum.cmd import bwa, samtools
import sys
import os
import argparse
import subprocess
import shlex

#argument parser
parser = argparse.ArgumentParser("Subsample reads from paired end reads")
parser.add_argument("R1", type=str, help="read 1 fastq")
parser.add_argument("R2", type=str, help="read 2 fastq")
parser.add_argument("start", type=str, help="smallest subsample size")
parser.add_argument("stop", help="biggest subsample size")
parser.add_argument("increment", help="increment size")
parser.add_argument("outdir", help="outdir")


args = parser.parse_args()
            
def deleteContent(pfile):
    pfile.seek(0)
    pfile.truncate()
        
def main():
    # create outputdir if it doesn't exit
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # run seqtk
    for nb_reads in range(int(args.start), int(args.stop) + int(args.increment),  int(args.increment)):
        nb_reads = str(nb_reads)
        print "---sampling " + nb_reads + " reads---"
        os.system("seqtk sample -s 42 " + args.R1 + " " + nb_reads + " > " + args.outdir + "/" + nb_reads + "_R1_subsample.fq")
        os.system("seqtk sample -s 42 " + args.R2 + " " + nb_reads + " > " + args.outdir + "/" + nb_reads + "_R2_subsample.fq")

if __name__ == '__main__':
    main()
