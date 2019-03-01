#!/usr/bin/env python3
__author__ = "Romain Guerillot"
__version__ = "0.1"
__email__ = "guerillot.romain@gmail.com"

import sys
import os
from subprocess import Popen, PIPE
import shutil
import argparse
import re

def copy_reference_make_index(ref, outpath):
    if not os.path.exists(outpath + '/reference'):
        os.mkdir(outpath + '/reference')
#    os.chdir(outpath + '/reference')
    shutil.copy(ref, outpath + '/reference/ref.fa')
    cmd = 'bwa index ' + outpath + '/reference/ref.fa'
    print(cmd)
    os.system(cmd)
#    os.chdir('..')    

def run_bwa(outpath, R1, R2, threads, refpath):
    if R2:
        cmd = 'bwa mem -t ' + threads + ' ' + refpath + ' ' + R1 + ' ' + R2 + ' | samtools sort > ' + outpath + '/all_reads.bam'
    else:
        cmd = 'bwa mem -t ' + threads + ' ' + refpath + ' ' + R1 + ' | samtools sort > ' + outpath + '/all_reads.bam'
    print(cmd)
    os.system(cmd)

def run_minimap2(outpath, R1, threads, ref):
    cmd = 'minimap2 -ax map-pb -t ' + threads + ' ' + ref + ' ' + R1 + ' | samtools sort > ' + outpath + '/all_reads.bam' 
    print(cmd)
    os.system(cmd)

def get_coverage(refpath, outpath, in_filename, out_filename):
    # get coverage at each position from bam alignment and generate userplot for artemis
    cmd = 'bedtools genomecov -d -ibam ' + outpath + '/' + in_filename + ' > ' + outpath + '/' + out_filename
    print(cmd)
    os.system(cmd)
    cmd = 'cut -f3 ' + outpath + '/' + out_filename + '>' + outpath + '/' + out_filename + '.userplot'
    print(cmd)
    os.system(cmd)

def get_split(outpath):
    # convert bam to sam, extract split reads, convert to bam and sort
    cmd = 'samtools view -h ' + outpath + '/all_reads.bam | ' + 'extractSplitReads_BwaMem -i stdin | samtools view -bS | samtools sort - > ' + outpath + '/split_reads.bam'
    print(cmd)
    os.system(cmd)
    # index sorted bam 
    cmd = 'samtools index ' + outpath + '/split_reads.bam ' + outpath + '/split_reads.bai'
    print(cmd)
    os.system(cmd)

def get_discordant(outpath):
    cmd = 'samtools view -b -F 1294 ' + outpath + '/all_reads.bam | samtools sort - > ' + outpath + '/discordant_reads.bam'
    print(cmd)
    os.system(cmd)
    # index sorted bam 
    cmd = 'samtools index ' + outpath + '/discordant_reads.bam ' + outpath + '/discordant_reads.bai'
    print(cmd)
    os.system(cmd)

def main():
    #parse arguments
    parser = argparse.ArgumentParser(description='Align reads to a reference genome and output genome coverage of all reads, split reads and discordant reads (userplot files can be imported into artemis')
    parser.add_argument('-r', metavar = 'REFERENCE', required = True, 
                help = 'Reference in fasta format')
    parser.add_argument('-o', metavar = 'OUTDIR', required = True,
                        help = 'Output directory')
    parser.add_argument('-R1', metavar = 'FASTQ', required = False,
                        help = 'Read 1')
    parser.add_argument('-R2', metavar = 'FASTQ', required = False,
                        help = 'Read 2 (not required for single end reads)', default = None)
    parser.add_argument('-pacbio', metavar = 'FASTQ', required = False,                                                                                                                                                              
                        help = 'To use for pacbio reads (align with minimap2)', default = None)
    parser.add_argument('-t', '--threads',
                        help = 'Number of CPUS (default 8)',
                        default = '8')
    parser.add_argument('-s', '--skip',
                        help = 'skip bwa alignment',
                        action = 'store_true')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

    
    # get absolute paths:
    refpath = os.path.abspath(args.o + "/reference/ref.fa")
    outpath = os.path.abspath(args.o)
    if args.R1:
        R1path = os.path.abspath(args.R1)
    else:
        R1path = ""
    if args.R2:
        R2path = os.path.abspath(args.R2)
    else:
        R2path = ""

    # create out directory if absent
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    # don't skip bwa alignment if skip is false 
    if not args.skip:
        if args.pacbio:
            run_minimap2(outpath, args.pacbio, args.threads, args.r)
        else:
            copy_reference_make_index(args.r, outpath)
            run_bwa(outpath, R1path, R2path, args.threads, refpath)
    
    # get splitted reads
    get_split(outpath)

    # get dicordant reads
    get_discordant(outpath)

    # get coverage and artemis userplot
    get_coverage(refpath, outpath, "all_reads.bam", "all_reads.cov")
    get_coverage(refpath, outpath, "split_reads.bam", "split_reads.cov")
    get_coverage(refpath, outpath, "discordant_reads.bam", "discordant_reads.cov")

if __name__ == '__main__':
    sys.exit(main())
