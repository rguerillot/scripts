#!/usr/bin/env python
from plumbum import local
from plumbum.cmd import bwa, samtools
import sys
import os
import argparse
import subprocess
import shlex

#argument parser
parser = argparse.ArgumentParser("Align reads to to a reference and output coverage and non reference allele count for each position")
parser.add_argument("-1", "--R1", type=str, help="read 1 fastq")
parser.add_argument("-2", "--R2", default = "", type=str, help="read 2 fastq (only use -1 if not paired end)")
parser.add_argument("reffa", type=str, help="fasta reference file")
parser.add_argument("refgff", type=str, help="gff3 reference file")
parser.add_argument("outdir", help="output directory")
parser.add_argument("-t", "--thread", default = "20", help="number of threads (default: 20")
parser.add_argument("-n", "--name", default = "sample", help="sample name to add in table (default: sample)")
parser.add_argument("-c", "--mincov", default = "0", help="minimum coverage to call snp (default: 0)")
parser.add_argument("-m", "--mapqual", default = "0", help="minimum read mapping quality (default: 0)")
parser.add_argument("-b", "--basequal", default = "0", help="minimum base quality to use for count (default: 0)")
parser.add_argument("-f", "--feature", default = "CDS", help="feature to use for counting non ref snp (default: CDS)")

args = parser.parse_args()
            
def deleteContent(pfile):
    pfile.seek(0)
    pfile.truncate()

def parse_acgt():
    with open( args.outdir + "/" + args.name + ".acgt", "r") as fin:
        lines = fin.readlines()
        sample = args.outdir.split("/")[-1]
        tab_file = "sample" + "\t" + str(lines.pop(0)).replace("\n","") + "\t" + "cov" + "\t" + "ref_call" + "\t" + "non_ref_call" + "\t" + "maj_call" + "\t" + "non_maj_call" + "\n"
        bed_file = ""
        for line in lines:
            values = line.split("\t")
            chrom = values[0]
            pos = values[1]
            refb = values[2]
            A = int(values[4])
            C = int(values[5])
            G = int(values[6])
            T = int(values[7])
            cov = int(A) + int(C) + int(G) + int(T)
            maj = int(max([A,C,G,T]))
            non_maj = cov - maj
            if refb == "A":
                ref = A
                non_ref = C+G+T
            if refb == "C":
                ref = C
                non_ref = A+G+T
            if refb == "G":
                ref = G
                non_ref = A+C+T
            if refb == "T":
                ref = T
                non_ref = A+C+G
            tab_file += sample + "\t" + "\t".join(values).replace("\n", "") + "\t" + str(cov) + "\t" + str(ref) + "\t" + str(non_ref) + "\t" + str(maj) + "\t" + str(non_maj) + "\n"
            bed_file += "\t".join([chrom,pos,pos,sample,str(non_ref),str(non_maj),str(cov)]) + "\n"
    with open(args.outdir + "/" + args.name + ".nref", 'w') as fout:
        fout.write(tab_file)
    with open(args.outdir + "/" + args.name + ".bed", 'w') as fout:
        fout.write(bed_file)
    
        
def main():
    # create outputdir if it doesn't exit
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # indexing reference
    print "---Indexing reference---"
    subprocess.call(["bwa", "index", args.reffa])

    # run bwa-mem
    print "---Aligning reads to reference---"
    os.system("bwa mem -t " + args.thread + " " + args.reffa + " " +  args.R1 + " " + args.R2 + " | " + "samtools sort -@" + args.thread + " -O BAM -o " + args.outdir + "/" + args.name + ".bam" + " - ")
    
    ## sam/bam conversion and sorting
    #print "---Converting bam to sam and sorting---"
    #os.system("samtools view -@ " + args.thread + " -q " + args.mapqual + " -b " + args.outdir + "/" + args.name + ".sam" + " | samtools sort -@ " + args.threads + " > " + args.outdir + "/" + args.name + ".bam")

    # run mpileup and sequenza-utils pilup2acgt
    print "---Counting all base call---"
    os.system("samtools mpileup -d 10000000 -Q 0 -f " + args.reffa + " " + args.outdir + "/" + args.name + ".bam" + "|" + "sequenza-utils.py pileup2acgt -p " + args.thread + " -q " + args.basequal + " -n " + args.mincov + " -" + "> " +  args.outdir + "/" + args.name + ".acgt") 
    # count non_ref call
    parse_acgt()

    # get CDS - snp intersection with gff
    #print "---Counting non-reference, non-majoritary and total base call for each " + args.feature + "---" 
    #os.system( "grep \#\# "+ args.refgff + " >" + args.outdir + "/reference.gff")
    #os.system( "grep "+ args.feature + " " + args.refgff + ">>"+ args.outdir + "/reference.gff")
    #os.system( "bedtools intersect -a " + args.outdir + "/" + args.name + ".bed " + "-b " + args.outdir + "/reference.gff " + "-wb >" + args.outdir + "/" + args.name + ".intersect")
    #os.system("/home/rguerillot/bin/groupBy -i " + args.outdir + "/" + args.name + ".intersect -g 4,16 -c 5 | sort | /home/rguerillot/bin/groupBy -g 1,2 -c 3 >" + args.outdir + "/" + args.name + "_" + args.feature + "_nref.count")
    #os.system("/home/rguerillot/bin/groupBy -i " + args.outdir + "/" + args.name + ".intersect -g 4,16 -c 6 | sort | /home/rguerillot/bin/groupBy -g 1,2 -c 3 >" + args.outdir + "/" + args.name + "_" + args.feature + "_nmaj.count")
    #os.system("/home/rguerillot/bin/groupBy -i " + args.outdir + "/" + args.name + ".intersect -g 4,16 -c 7 | sort | /home/rguerillot/bin/groupBy -g 1,2 -c 3 >" + args.outdir + "/" + args.name + "_" + args.feature + "_cov.count")
    #os.system("paste " + args.outdir + "/" + args.name + "_" + args.feature + "_nref.count " + args.outdir + "/" + args.name + "_" + args.feature + "_nmaj.count "+ args.outdir + "/" + args.name + "_" + args.feature + "_cov.count | cut -f 1,2,3,6,9 > " + args.outdir + "/" + args.name + ".count")

if __name__ == '__main__':
    main()
