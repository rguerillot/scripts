#!/home/linuxbrew/.linuxbrew/opt/python/bin/python3.7

'''
    Uses python3.
    email: romain.guerillot@hotmail.fr
    Authors: Romain Guerillot
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import sys
import subprocess
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Alphabet import SingleLetterAlphabet
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import operator
from io import StringIO
import os, errno
import argparse
import shlex
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature

#argument parser
parser = argparse.ArgumentParser(description = "Script to identify and annotate inverted repeats (potentially promoting chromosomal inversion) from bacterial complete genome assemblies")
parser.add_argument("genbank", type=str, help="Complete chromosome assembly in genbank format")
parser.add_argument("outdir", help="Output directory")
parser.add_argument("-i", "--identity", type = int,  default = 90, help="Minimum %% identity of the inverted repeats (default: 90)")
parser.add_argument("-s", "--size", type = int, default = 1000, help="Minimum size in bp of the inverted repeats (default: 1000)")
parser.add_argument("-f", "--filter", type = str, default = "yes", help="Set to no if you want inverted repeats detected in the same arm of the chromosome (default: \"yes\")")
args = parser.parse_args()

#reverse complement DNA sequence
def revcomp(seq):
    for base in seq:
        if base not in 'ATCGatcg':
            print("Error: NOT a DNA sequence")
            return None
    seq1 = 'ATCGTAGCatcgtagc'
    seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
    return "".join([seq_dict[base] for base in reversed(seq)])

#convert genbank file to fasta file
def gbk_to_fasta(gbk_file, outputfolder):
    input_handle = open(gbk_file, "r")
    output_handle = open(outputfolder + "/" + "genome.fa", "w")
    first_record = SeqIO.parse(input_handle, "genbank").__next__()
    SeqIO.write(first_record, output_handle, "fasta")
    input_handle.close()
    output_handle.close()

#make blast database
def make_blast_db(fasta_file, outputfolder):
    print("Making blast database...")
    cmd = "makeblastdb -in " + str(fasta_file) + " -input_type fasta -dbtype nucl -title blast.bdb -parse_seqids -out " + outputfolder + "/blast_db"
    print("Running " + cmd)
    os.system(cmd)

#do blastn search
def local_blastn(my_query, my_sbjct, outputfolder):
    make_blast_db(my_sbjct, outputfolder)
    print("Running Blast ...")
    cmd = "blastn -db " + outputfolder + "/blast_db " + "-query " + my_query + " -outfmt 5 " + "-out " + outputfolder + "/blast_result.xml"
    print("Running " + cmd)    
    os.system(cmd)

#identify inverted repeats from blast results in xml format and create new features
def get_inverted_repeats(xml_file, size, identity, framing):
    count = 0
    new_features=[]
    fasta_seq=""
    print("Parsing blast result: " + xml_file)
    result_handle=open(str(xml_file),"r")
    blast_record = NCBIXML.parse(result_handle)
    for record in blast_record:
        for alignment in record.alignments:
            rep_terminus_coor = alignment.length/2
            for hsp in alignment.hsps:
                ident = float(hsp.identities)/float(hsp.align_length)*100
                if int(hsp.align_length) >= size and ident >= identity: # filter alignment size and identity
                    if (hsp.query_start - hsp.query_end) < 0 and (hsp.sbjct_start - hsp.sbjct_end) > 0:  # filter inverted matches
                        if framing == "no":
                            count +=1
                            fasta1 = ">" + alignment.title + "_" + str(count) + "For_length:" + str(hsp.align_length) + "_id:" + str(int(ident*100)) + "_coor:" + str(hsp.query_start) + ":" + str(hsp.query_end) + "\n" + hsp.query.replace("-", "") + "\n"
                            fasta2 = ">" + alignment.title + "_" + str(count) + "Rev_length:" + str(hsp.align_length) + "_id:" + str(int(ident*100)) + "_coor:" + str(hsp.sbjct_end) + ":" + str(hsp.sbjct_start) + "\n" + hsp.sbjct.replace("-", "") + "\n"
                            fasta_seq += fasta1
                            fasta_seq += fasta2
                            feat1 = create_feature(hsp.query_start, hsp.query_end, get_strand(hsp.query_start, hsp.query_end), ftype = "repeat_region", fqual = {"locus_tag":"inverted repeat " + str(count)})  
                            feat2 = create_feature(hsp.sbjct_end, hsp.sbjct_start, get_strand(hsp.sbjct_start, hsp.sbjct_end), ftype = "repeat_region", fqual = {"locus_tag":"inverted repeat " + str(count)})
                            new_features.append(feat1)
                            new_features.append(feat2)
                        elif framing == "yes" and hsp.query_end < rep_terminus_coor and hsp.sbjct_end > rep_terminus_coor: # filter matches that frame replication terminus
                            count +=1
                            fasta1 = ">" + alignment.title + "_" + str(count) + "For_length:" + str(hsp.align_length) + "_id:" + str(int(ident*100)) + "_coor:" + str(hsp.query_start) + ":" + str(hsp.query_end) + "\n" + hsp.query.replace("-", "") + "\n"
                            fasta2 = ">" + alignment.title + "_" + str(count) + "Rev_length:" + str(hsp.align_length) + "_id:" + str(int(ident*100)) + "_coor:" + str(hsp.sbjct_end) + ":" + str(hsp.sbjct_start) + "\n" + hsp.sbjct.replace("-", "") + "\n"
                            fasta_seq += fasta1
                            fasta_seq += fasta2
                            feat1 = create_feature(hsp.query_start, hsp.query_end, get_strand(hsp.query_start, hsp.query_end), ftype = "repeat_region", fqual = {"locus_tag":"inverted repeat " + str(count)})  
                            feat2 = create_feature(hsp.sbjct_end, hsp.sbjct_start, get_strand(hsp.sbjct_start, hsp.sbjct_end), ftype = "repeat_region", fqual = {"locus_tag":"inverted repeat " + str(count)})
                            new_features.append(feat1)
                            new_features.append(feat2)
    print("Found " + str(int(count/2)) + " inverted repeats") 
    return (new_features, fasta_seq)

#define DNA strand from start and end coordinate
def get_strand(start, end):
    if start - end > 0:
        strand = -1
    if start - end < 0:
        strand = 1
    return strand

#create biopython SeqFeature objects                                                
def create_feature(fstart, fend, fstrand, ftype, fqual):
    my_feature_location = FeatureLocation(fstart, fend, strand = fstrand)
    my_feature = SeqFeature(my_feature_location,type=ftype, qualifiers=fqual)
    return my_feature

#add feature to original gbk file
def add_features_to_gbk(gbk_file_path, file_name, new_features, outputfolder):
    input_handle = open(gbk_file_path, "r")
    output_handle = open(outputfolder + "/" + file_name + "_inverted_repeat.gbk", "w")
    first_record = SeqIO.parse(input_handle, "genbank").__next__()
    first_record.features.extend(new_features)
    first_record.name = file_name[0:15]
    SeqIO.write(first_record, output_handle, "genbank")
    input_handle.close()
    output_handle.close()

def main():
    # argument and variables
    args = parser.parse_args()
    genbank = os.path.abspath(args.genbank)
    outdir = os.path.abspath(args.outdir)

    # create outputdir if it doesn't exit
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # parse filename
    file_name = os.path.basename(genbank)
        
    # convert gbk to fasta
    gbk_to_fasta(genbank, outdir) # convert first genbank record and write it as fasta file
    
    # run recipocal blastn 
    local_blastn(outdir + "/genome.fa", outdir + "/genome.fa", outdir) # blast query file against itself
    
    # get inverted repeats features
    ir_features = get_inverted_repeats(outdir + "/blast_result.xml", args.size, args.identity, args.filter)
    
    # write new genbank files with inverted repeats annotated
    add_features_to_gbk(genbank, file_name, ir_features[0], outdir)
    
    # write mutifasta with nverted repeats sequences
    with open(outdir + "/" + file_name + "_inverted_repeat.fa", "w") as text_file:
        text_file.write(ir_features[1])


if __name__ == '__main__':
    main() 
