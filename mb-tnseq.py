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
    shutil.copy(ref, outpath + '/reference/ref.fa')
    cmd = 'bwa index ' + outpath + '/reference/ref.fa'
    print(cmd)
    os.system(cmd)

def extract_bc(R1, barcode, outpath, R2 = None):
    if R2 != None:
        cmd = 'umi_tools extract -p ' + barcode + ' -I ' + R1 + ' -S ' + outpath + '/R1_bc.fq.gz --read2-in=' + R2 + ' --read2-out=' + outpath + '/R2_bc.fq.gz'
    else:
        cmd = 'umi_tools extract -p ' + barcode + ' -I ' + R1 + ' -S ' + outpath + '/R1_bc.fq.gz'
    print(cmd)
    os.system(cmd)

def run_bwa(outpath, R1, threads, refpath, R2 = None
):
    if R2 != None:
        cmd = 'bwa mem -t ' + threads + ' ' + refpath + ' ' + R1 + ' ' + R2 + ' | samtools sort > ' + outpath + '/all_reads.bam'
    else:
        cmd = 'bwa mem -t ' + threads + ' ' + refpath + ' ' + R1 + ' | samtools sort > ' + outpath + '/all_reads.bam'
    print(cmd)
    os.system(cmd)
    #index bam
    cmd = 'samtools index ' + outpath + '/all_reads.bam'
    print(cmd)
    os.system(cmd)

def dedup_bam(outpath):
    cmd = 'umi_tools dedup --output-stats=deduplicated --stdin=' + outpath + '/all_reads.bam ' + '--stdout=' + outpath + '/all_reads_dedup.bam'
    print(cmd)
    os.system(cmd)
    
def bam_to_bed(outpath):
    cmd = 'bedtools bamtobed -i ' + outpath + '/all_reads.bam > ' + outpath + '/all_reads.bed'
    print(cmd)
    os.system(cmd)

    cmd = 'bedtools bamtobed -i ' + outpath + '/all_reads_dedup.bam > ' + outpath + '/all_reads_dedup.bed'
    print(cmd)
    os.system(cmd)
    
def bedtools_closest(outpath, bedfile):
    cmd = 'bedtools closest -iu -b ' + bedfile + ' -a ' + outpath + '/all_reads.bed -D a | awk -F \'\\t\' \'{ if($12!=-1)print $0}\'| awk -F \'\\t\' \'{ if($12<600)print $0}\' > ' + outpath + '/NEB_read_count_raw.tab'
    print(cmd)
    os.system(cmd)

    cmd = 'bedtools closest -iu -b ' + bedfile + ' -a ' + outpath + '/all_reads_dedup.bed -D a | awk -F \'\\t\' \'{ if($12!=-1)print $0}\'| awk -F \'\\t\' \'{ if($12<600)print $0}\' > ' + outpath + '/NEB_read_dedup_count_raw.tab'
    print(cmd)                                                          
    os.system(cmd)

def count_insertion(outpath):
     cmd = 'cut -f10 ' + outpath + '/NEB_read_count_raw.tab | sort | uniq -c | sort -n | awk \'{$1=$1;print}\' > ' + outpath + '/NEB_read_count.tab'
     print(cmd)
     os.system(cmd)

     cmd = 'cut -f10 ' + outpath + '/NEB_read_dedup_count_raw.tab | sort | uniq -c | sort -n | awk \'{$1=$1;print}\' > ' + outpath + '/NEB_read_dedup_count.tab'
     print(cmd)
     os.system(cmd)
    
def get_coverage(refpath, outpath, in_filename, out_filename, stranded = False):
    # get coverage at each position from bam alignment and generate userplot for artemis
    cmd = 'bedtools genomecov -d -ibam ' + outpath + '/' + in_filename + ' > ' + outpath + '/' + out_filename + ".cov"
    print(cmd)
    os.system(cmd)
    cmd = 'cut -f3 ' + outpath + '/' + out_filename + ".cov" + '>' + outpath + '/' + out_filename + '.userplot'
    print(cmd)
    os.system(cmd)
    if stranded == True:
        cmd = 'bedtools genomecov -d -strand + -ibam ' + outpath + '/' + in_filename + ' > ' + outpath + '/' + out_filename + "_plus.cov"
        print(cmd)
        os.system(cmd)
        cmd = 'cut -f3 ' + outpath + '/' + out_filename + "_plus.cov" + '>' + outpath + '/' + out_filename + '_plus.userplot'
        print(cmd)
        os.system(cmd)
        cmd = 'bedtools genomecov -d -strand - -ibam ' + outpath + '/' + in_filename + ' > ' + outpath + '/' + out_filename + "_minus.cov"
        print(cmd)
        os.system(cmd)
        cmd = 'cut -f3 ' + outpath + '/' + out_filename + "_minus.cov" + '>' + outpath + '/' + out_filename + '_minus.userplot'
        print(cmd)
        os.system(cmd)

def main():
    #parse arguments
    parser = argparse.ArgumentParser(description='MB-Tnseq analysis pipeline: 1. Extract barcode from fastq sequence and add to header (umi_tools extract); 2. Map reads to genome (bwa); 3. Deduplicate reads with same barrcode at same coordinate (umitools dedup); 4. Associate reads to Tn insertion provided in a bed file (bedtools closest); 5. Count number of reads (<=> unique barcoded bacteria) for every Tn insertion')

    parser.add_argument('-r', metavar = 'REFERENCE', required = True,
                        help = 'Reference genome in fasta format (for NEBRASKA Tn library use USA300_FPR3757.fa)')
    parser.add_argument('-a', metavar = 'ANNOTATION', required = True,
                        help = 'Annotation file with the coordinate of the Tn insertion of the mutant library in bed format (for NEBRASKA Tn library use NEB_pos.bed)')
    parser.add_argument('-b', metavar = 'BARCODE', required = False,
                        help = 'barcode patern for umitools extract (default for NEBRASKA Tn library NNNNNNNN',
                        default="NNNNNNNN")
    parser.add_argument('-o', metavar = 'OUTDIR', required = True,
                        help = 'Output directory')
    parser.add_argument('-R1', metavar = 'FASTQ', required = False,
                        help = 'Read 1')
    #parser.add_argument('-R2', metavar = 'FASTQ', required = False,
    #                    help = 'Read 2 (not required for single end reads)', default = None)
    parser.add_argument('-t', '--threads',
                        help = 'Number of CPUS (default 8)',
                        default = '8')
    #parser.add_argument('-sa', '--skip_aln',
    #                    help = 'skip bwa alignment',
    #                    action = 'store_true')
    #parser.add_argument('-st', '--strand',
    #                    help = 'report coverage for each strand',
    #                    action = 'store_true')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

    
    # get absolute paths:
    refpath = os.path.abspath(args.o + "/reference/ref.fa")
    outpath = os.path.abspath(args.o)
    annotpath = os.path.abspath(args.a)
    if args.R1:
        R1path = os.path.abspath(args.R1)
    else:
        R1path = ""
    #if args.R2:
    #    R2path = os.path.abspath(args.R2)
    #else:
    #    R2path = ""

    # create out directory if absent
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    # extract barcode
    extract_bc(R1path, args.b, outpath)

    # copy reference genome, run bwa alignments and index bam file 
    copy_reference_make_index(args.r, outpath)
    #run_bwa(outpath, outpath+"/R1_bc.fq.gz", outpath+"/R2_bc.fq.gz", args.threads, refpath)
    run_bwa(outpath, outpath+"/R1_bc.fq.gz", args.threads, refpath)

    # run deduplication or mapped reads using barcodes
    dedup_bam(outpath)

    # get strand coverage and artemis userplot
    get_coverage(refpath, outpath, "all_reads.bam", "all_reads", True)
    get_coverage(refpath, outpath, "all_reads.bam", "all_reads", False)
    get_coverage(refpath, outpath, "all_reads_dedup.bam", "all_reads_dedup", True)
    get_coverage(refpath, outpath, "all_reads_dedup.bam", "all_reads_dedup", False)

    # convert bam to bed
    bam_to_bed(outpath)

    # asign every read to insertion
    bedtools_closest(outpath, annotpath)

    # count reads for every insertion
    count_insertion(outpath)
    
if __name__ == '__main__':
    sys.exit(main())
