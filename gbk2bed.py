#!/usr/bin/env python3
# Author: Romain Guerillot
# Date: September 2020


from Bio import SeqIO
import os
import sys
import pdb
import argparse


parser = argparse.ArgumentParser(description='Convert gbk to bed format')
parser.add_argument('gbk', type=str,
                    help='genbank file')
parser.add_argument('--feature_type', type=str,
                    default = "gene",
                    help='eg. gene, CDS, tRNA...')  
parser.add_argument('--feature_qualifier', type=str,
                    default = "locus_tag",
                    help='eg. locus_tag, old_locus_tag, name...')  

args = parser.parse_args()

def main():
    outf = open(os.path.basename(args.gbk).split(".")[0] + "__" + args.feature_type + "__" + args.feature_qualifier + '.bed', 'w')
    for record in SeqIO.parse(open(args.gbk, "r"), "genbank") :
        for feature in record.features:
            if feature.type == args.feature_type:
                start = feature.location.start.position
                stop = feature.location.end.position
                try:
                    name = feature.qualifiers[args.feature_qualifier][0]
                except:
                    sys.stderr.write("Warning: feature at position " + str(start) + " doesn't have the feature qualifier " + args.feature_qualifier + " and has been skipped\n")
                    pass
                if feature.strand < 0:
                    strand = "-"
                else:
                    strand = "+"
                bed_line = "chrom1\t{0}\t{1}\t{2}\t100\t{3}\n".format(start, stop, name, strand)
                print(bed_line, end = "")
                outf.write(bed_line)
    outf.close()


if __name__ == '__main__':
    main()
