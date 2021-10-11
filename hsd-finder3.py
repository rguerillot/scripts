from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from operator import attrgetter
import copy

#cd-hit -i REBASE_all_prot.fa -o REBASE_all_prot_r95.fa -c 0.95
#nohup nice tblastn -db ./all_Sepi_blastdb -query REBASE_all_prot_r95.fa -out REBASE_r95_vs_all_Sepi.xml -outfmt 5 -num_threads 12

class hsd_hit:
    def __init__(self, query_name, sbjct_name, hsp_prot_seq, strain, contig, query_start_end, sbjct_start_end, bit_score,  id_frac, posit_frac, coverage, evalue, cluster = []):
        self.query_name=query_name
        self.sbjct_name=sbjct_name
        self.hsp_prot_seq=hsp_prot_seq
        self.strain=strain
        self.contig=contig
        self.query_start_end=query_start_end
        self.sbjct_start_end=sbjct_start_end
        self.bit_score = bit_score
        self.id_frac=id_frac
        self.posit_frac=posit_frac
        self.coverage=coverage
        self.evalue=evalue
        self.cluster=cluster
        

        if "restriction" in self.query_name:
            self.RMsu = "Restriction"
        elif "methyltransferase" in self.query_name:
            self.RMsu = "Methyltransferase"
        elif "specificity" in self.query_name:
            self.RMsu = "Specificity"
        else:
            self.RMsu = "Unknown"
        if "TypeI" in self.query_name and "TypeII" not in self.query_name and "TypeIII" not in self.query_name and "TypeIV" not in self.query_name:
            self.RMtype="TypeI"
        elif "TypeII" in self.query_name and "TypeIII" not in self.query_name:
            self.RMtype="TypeII"
        elif "TypeIII" in self.query_name :
            self.RMtype="TypeIII"
        elif "TypeIV" in self.query_name :
            self.RMtype="TypeIV"
        else:
            self.RMtype="Unknown"

    def __str__(self):
        string = ""
        dico = self.__dict__
        for key in dico:
            string += str(dico[key]) + "\t" 
        return string
    
    def print_attribut_name(self):
        dico = self.__dict__
        string=""
        for key in dico:
            string += str(key) + "\t" 
        return string

    def __eq__(self, other):
        return self.__dict__ == other.__dict__
    
    def __ne__(self, other):
        """Define a non-equality test"""
        return not self.__eq__(other)

    def get_prot_seq(self, gbkfile):
        records = SeqIO.parse("gbkfile", "genbank")
        for feature in records.features:
            if my_snp in feature:
                print("%s %s" % (feature.type, feature.qualifiers.get('db_xref')))
        pass

class Interval(object):
    """
    Represents an interval. 
    Defined as half-open interval [start,end), which includes the start position but not the end.
    Start and end do not have to be numeric types. 
    """
    
    
    def __init__(self, start, end):
        "Construct, start must be <= end."
        if start > end:
            raise ValueError('Start (%s) must not be greater than end (%s)' % (start, end))
        self._start = start
        self._end = end
        
         
    start = property(fget=lambda self: self._start, doc="The interval's start")
    end = property(fget=lambda self: self._end, doc="The interval's end")
     

    def __str__(self):
        "As string."
        return '[%s,%s)' % (self.start, self.end)
    
    
    def __repr__(self):
        "String representation."
        return '[%s,%s)' % (self.start, self.end)
    
    
    def __cmp__(self, other):
        "Compare."
        if None == other:
            return 1
        start_cmp = cmp(self.start, other.start)
        if 0 != start_cmp:
            return start_cmp
        else:
            return cmp(self.end, other.end)


    def __hash__(self):
        "Hash."
        return hash(self.start) ^ hash(self.end)
    
    
    def intersection(self, other):
        "Intersection. @return: An empty intersection if there is none."
        if self > other:
            other, self = self, other
        if self.end <= other.start:
            return Interval(self.start, self.start)
        return Interval(other.start, self.end)


    def hull(self, other):
        "@return: Interval containing both self and other."
        if self > other:
            other, self = self, other
        return Interval(self.start, other.end)
    

    def overlap(self, other):
        "@return: True iff self intersects other."
        if self > other:
            other, self = self, other
        return self.end > other.start
         

    def __contains__(self, item):
        "@return: True iff item in self."
        return self.start <= item and item < self.end
         

    def zero_in(self):
        "@return: True iff 0 in self."
        return self.start <= 0 and 0 < self.end
         

    def subset(self, other):
        "@return: True iff self is subset of other."
        return self.start >= other.start and self.end <= other.end
         

    def proper_subset(self, other):
        "@return: True iff self is proper subset of other."
        return self.start > other.start and self.end < other.end
         

    def empty(self):
        "@return: True iff self is empty."
        return self.start == self.end
         

    def singleton(self):
        "@return: True iff self.end - self.start == 1."
        return self.end - self.start == 1
    
    
    def separation(self, other):
        "@return: The distance between self and other."
        if self > other:
            other, self = self, other
        if self.end > other.start:
            return 0
        else:
            return other.start - self.end    


def parse_hsd_blast_xml(xml_file, min_match_length=0, min_coverage=0, min_identity=0, min_positive_fraction=0, max_evalue=0.005):
    
    list_hsd_hit=[]
    result_handle=open(str(xml_file),"r")
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
   # i = 0
   # while i < 1000:
   #     i += 1
   #     blast_record = blast_records.next()
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                query_length = float(blast_record.query_length)
                match_length = len(hsp.match)
                coverage = len(hsp.match)/query_length
                positives = float(hsp.positives)
                identical =  float(hsp.identities)
                positive_fraction = positives/len(hsp.query)
                identical_fraction = positives/len(hsp.query)
                if match_length > min_match_length and identical_fraction > min_identity and coverage > min_coverage and positive_fraction > min_positive_fraction and hsp.expect < max_evalue:
                    prot = hsp.sbjct
                    query_name = blast_record.query
                    sbjct_name = str(alignment.title)
                    hit_seq = hsp.sbjct
                    bitscore= hsp.bits
                    contig = alignment.title.split("_")[2].split(" ")[0]
                    strain = alignment.title.split("_")[1]
                    queryse = Interval(hsp.query_start, hsp.query_end)
                    sbjctse = Interval(hsp.sbjct_start, hsp.sbjct_end)
                    filtered_hit = hsd_hit(query_name, sbjct_name, prot, strain, contig, queryse, sbjctse, bitscore, identical_fraction, positive_fraction, coverage, hsp.expect)
                    list_hsd_hit.append(filtered_hit)
    return list_hsd_hit


class strain_hsd:
    def __init__(self, strain, hsd1="", hsd2="", hsd3=""):
        self.strain = strain
        self.hsd1 = hsd1
        self.hsd2 = hsd2
        self.hsd3 = hsd3

def get_strain_dic_of_hsd(list_of_hsd_hit, min_id_frac, min_coverage):
    dict_strain={}
    for hsd in list_of_hsd_hit :
        if hsd.id_frac > min_id_frac  and hsd.coverage > min_coverage  :
            strain = hsd.strain
            hsdgene= hsd.hsd_type + hsd.hsd_allele #+ "(cov" + str(int(hsd.coverage*100)) + "%)" + "(id" + str(int(hsd.id_frac*100)) + "%)"# + str(hsd.query_start_end)
            if dict_strain.has_key(strain):
                dict_strain[strain]=dict_strain[strain] + "\t" + hsdgene
            else:
                dict_strain[strain]=hsdgene
    return dict_strain

def get_strain_hsd_tab():
    list_strain_hsd=[]
    for key in dict_strain:
        strain = key
        hsd1=""
        hsd2=""
        hsd3=""
        if "hsdM1" in dict_strain[key] :
            if hsd1 != "":
                hsd1+="M"
            else:
                hsd1="M"
        if "hsdR1" in dict_strain[key] :
            if hsd1 != "":
                hsd1+="R"
            else:
                hsd1="R"
        if "hsdS1" in dict_strain[key] :
            if hsd1 != "":
                hsd1+="S"
            else:
                hsd1="S"
        if "hsdM2" in dict_strain[key] :
            if hsd2 != "":
                hsd2+="M"
            else:
                hsd2="M"
        if "hsdR2" in dict_strain[key] :
            if hsd2 != "":
                hsd2+="R"
            else:
                hsd2="R"
        if "hsdS2" in dict_strain[key] :
            if hsd2 != "":
                hsd2+="S"
            else:
                hsd2="S"
        if "hsdS3" in dict_strain[key] :
            if hsd3 != "":
                hsd3+="S"
            else:
                hsd3="S"
    return list_strain_hsd.append(strain_hsd(strain, hsd1, hsd2, hsd3))


def get_fasta_in_multifasta(multifasta, header_name):

    fasta_handle=(open(multifasta, "r"))
    frecord_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))
    for key in frecord_dict :
        if key == queryname:
            fasta =  ">" + key + "\n" + frecord_dict[key].seq + "\n"
            return fasta
    fast_handle.close()

def get_unique_hsd_hit(list_hsd):
    unique_hsd=[]
    for hsd in list_hsd:
        if hsd in unique_hsd:
            continue
        else:
            unique_hsd.append(hsd)
    return unique_hsd    

def get_best_hit_at_pos(list_hsd):
    best_hsd_list = [hsd for hsd in list_hsd if iambesthsd(hsd, list_hsd)]
    return best_hsd_list

def iambesthsd(my_hsd, hsd_list):
    for hsd in hsd_list:
        if my_hsd.strain == hsd.strain and my_hsd.contig == hsd.contig and my_hsd.RMtype == hsd.RMtype and my_hsd.RMsu == hsd.RMsu and my_hsd.sbjct_start_end.overlap(hsd.sbjct_start_end) and my_hsd.bit_score < hsd.bit_score:
            return False
            #if len(my_hsd.hsp_prot_seq) < len(hsd.hsp_prot_seq): #exclude smaller overlapping blast hit
            #    return False
            #if len(my_hsd.hsp_prot_seq) == len(hsd.hsp_prot_seq) and my_hsd.evalue > hsd.evalue:
            #    return False #exclude hit with higher evalue if hit of same size
            #if my_hsd.evalue == hsd.evalue and len(my_hsd.hsp_prot_seq) == len(hsd.hsp_prot_seq) and my_hsd.id_frac < hsd.id_frac:
            #    return  False
            #if my_hsd.evalue == hsd.evalue and len(my_hsd.hsp_prot_seq) == len(hsd.hsp_prot_seq) and my_hsd.id_frac == hsd.id_frac and my_hsd.posit_frac < hsd.posit_frac:
            #    return  False
    return True
                    
def get_clustered_hsd(list_hsd, max_distance, cluster_size):
    clustered_hsd_list = [hsd for hsd in list_hsd if iamclustered(hsd, list_hsd, max_distance, cluster_size)]
    return clustered_hsd_list

def iamclustered(my_hsd, hsd_list, max_distance, cluster_size):
    cluster_size=1
    for hsd in hsd_list:
        if my_hsd.strain == hsd.strain and my_hsd.contig == hsd.contig and my_hsd.sbjct_start_end.separation(hsd.sbjct_start_end) < max_distance and my_hsd != hsd:
            cluster_size += 1
            if cluster_size >= cluster_size:
                return True
    return False

class RM:
    def __init__(self, hsd1, hsd2=None, hsd3=None):
        
        if hsd2 and not hsd3:
            sorthsd=sorted([hsd1, hsd2], key = attrgetter("RMsu"))
            self.hsd1 = sorthsd[0]
            self.hsd2 = sorthsd[1]
            self.hsd3 = None
            self.RMsus = [sorthsd[0].RMsu, sorthsd[1].RMsu]
            self.name = [sorthsd[0].query_name, sorthsd[1].query_name]
            self.size = 2
            if hsd1.RMtype == "TypeII" or hsd1.RMtype == "TypeIII":
                self.complete = True
            else:
                self.complete = None
        elif hsd2 and hsd3:
            sorthsd=sorted([hsd1, hsd2, hsd3], key = attrgetter("RMsu"))
            self.hsd1 = sorthsd[0]
            self.hsd2 = sorthsd[1]
            self.hsd3 = sorthsd[2]
            self.RMsus = [sorthsd[0].RMsu, sorthsd[1].RMsu, sorthsd[2].RMsu] 
            self.name = [sorthsd[0].query_name, sorthsd[1].query_name, sorthsd[2].query_name]
            self.size = 3
            self.complete = True
        elif hsd1.RMtype == "TypeIV":
            self.hsd1 = hsd1
            self.hsd2 = None
            self.hsd3 = None
            self.RMsus = [hsd1.RMsu]
            self.name = [hsd1.query_name]
            self.size = 1
            self.complete = True
        else :
            self.hsd1 = hsd1
            self.hsd2 = None
            self.hsd3 = None
            self.RMsus = [hsd1.RMsu]
            self.name = [hsd1.query_name]
            self.size = 1
            self.complete = False
        self.strain = hsd1.strain
        self.contig = hsd1.contig
        self.type = hsd1.RMtype
        self.recompleted=False
        self.count=None
        if self.size == 1:
            self.coor = self.hsd1.sbjct_start_end
        if self.size == 2:
            self.coor = self.hsd1.sbjct_start_end.hull(self.hsd2.sbjct_start_end)
        if self.size == 3:
            self.coor = self.hsd1.sbjct_start_end.hull(self.hsd2.sbjct_start_end.hull(self.hsd3.sbjct_start_end))
    
    def __repr__(self):
        return self.strain + "\t" + self.type + "\t" + str(self.RMsus) + "\t" + str(self.name) + "\t" + str(self.complete) + "\t" + str(self.size) + "\t" + str(self.coor) + "\n"

    def __eq__(self, other):
        return self.name == other.name and self.coor == other.coor and self.strain == other.strain and self.contig == other.contig and self.type == other.type and self.RMsus == other.RMsus
    def __ne__(self, other):
        return (not self.__eq__(other))

    def __hash__(self):
        return hash(self.__repr__())

    def recomplete(self, list_hsd, list_prototype_RM):
        compRM = False
        if self.complete == False:
            if self.type == "TypeI":
                for RM_ in list_prototype_RM:
                    if self.type == RM_.type:
                        hsdtofindA=None
                        hsdtofindB=None
                        matchinghsd=None
                        for numhsd, hsdname in enumerate(RM_.name):
                            if self.name == hsdname :
                                if numhsd == 0:
                                    matchinghsd=RM_.hsd1
                                    hsdtofindA=RM_.hsd2
                                    hsdtofindB=RM_.hsd3
                                if numhsd == 1:
                                    matchinghsd=RM_.hsd2
                                    hsdtofindA=RM_.hsd1
                                    hsdtofindB=RM_.hsd3
                                if numhsd == 2:
                                    matchinghsd=RM_.hsd3
                                    hsdtofindA=RM_.hsd1
                                    hsdtofindB=RM_.hsd2
                                for hsd in list_hsd:
                                    findhsdA=None
                                    findhsdB=None
                                    if self.hsd1.strain == hsd.strain and hsdtofindA.query_name == hsd.query_name:
                                        findhsdA=hsd
                                    if self.hsd1.strain == hsd.strain and hsdtofindB.query_name == hsd.query_name:
                                        findhsdB=hsd
                                if findhsdA and findhsdB:
                                    compRM = RM(self.hsd1, findhsdA, findhsdB)
                                    compRM.recompleted = True
                                    compRM.complete = True
                                    return compRM
                                if findhsdA and not findhsdB:
                                    compRM = RM(self.hsd1, findhsdA)
                                    compRM.recompleted = True
                                    compRM.complete = False
                                if findhsdB and not findhsdA:
                                    compRM = RM(self.hsd1, findhsdB)
                                    compRM.recompleted = True
                                    compRM.complete = False
            
            if self.type == "TypeII" or self.type == "TypeIII":
                for RM_ in list_prototype_RM:
                    hsdtofindA=None
                    matchinghsd=None
                    if RM_.type == self.type:
                        for numhsd, hsdname in enumerate(RM_.name):
                            if self.name == hsdname :
                                if numhsd == 0:
                                    matchinghsd=RM_.hsd1
                                    hsdtofindA=RM_.hsd2
                                if numhsd == 1:
                                    matchinghsd=RM_.hsd2
                                    hsdtofindA=RM_.hsd1
                                for hsd in list_hsd:
                                    findhsdA=None
                                    if self.hsd1.strain == hsd.strain and hsdtofindA.query_name == hsd.query_name:
                                        findhsdA=hsd
                                        compRM = RM(self.hsd1, findhsdA, findhsdB)
                                        compRM.recompleted = True
                                        return compRM
        if compRM:
            return compRM
        else :
            return self
                                    


    def my_occurence_in_listRM(self,list_RM):
        count=0
        for RM in list_RM:
            if self.name == RM.name:
                count+=1
        self.count = count
        

def get_complete_RM(list_hsd, max_distance, sorting_attribute=None):

    list_hsd_tI = [hsd for hsd in list_hsd if hsd.RMtype == "TypeI"]

    list_hsd_tII_III = [hsd for hsd in list_hsd if hsd.RMtype == "TypeII" or hsd.RMtype == "TypeIII"]

    RM_IV = [RM(hsd1) for hsd1 in list_hsd if hsd1.RMtype == "TypeIV"]
    
    RM_II_III = [RM(hsd1,hsd2) for hsd1 in list_hsd_tII_III for hsd2 in list_hsd_tII_III if hsd1.strain == hsd2.strain and hsd1.contig == hsd2.contig and hsd1.RMsu != hsd2.RMsu and hsd1.RMtype == hsd2.RMtype and hsd1.sbjct_start_end.separation(hsd2.sbjct_start_end) < max_distance]
    
    RM_I = [RM(hsd1,hsd2,hsd3) for hsd1 in list_hsd_tI for hsd2 in list_hsd_tI for hsd3 in list_hsd_tI if hsd1.strain == hsd2.strain == hsd3.strain and hsd1.contig == hsd2.contig == hsd3.contig and hsd1.RMsu != hsd2.RMsu != hsd3.RMsu != hsd1.RMsu and hsd1.sbjct_start_end.separation(hsd2.sbjct_start_end) < max_distance and hsd2.sbjct_start_end.separation(hsd3.sbjct_start_end) < max_distance]
                      
    return sorted(list(set(RM_I + RM_II_III + RM_IV)) , key=attrgetter("strain", "type", "name"))
   
def get_hsd_not_in_RM(listhsd, list_RM):
    list_hsd = copy.deepcopy(listhsd)

    hsd_in_all_RM = []
    for RM in list_RM:
        if RM.hsd1:
            hsd_in_all_RM.append(RM.hsd1)
        if RM.hsd2:
            hsd_in_all_RM.append(RM.hsd2)
        if RM.hsd3:
            hsd_in_all_RM.append(RM.hsd3)
            
    for hsd in list_hsd: 
        if hsd in hsd_in_all_RM:
            list_hsd.remove(hsd)

    return list_hsd

def get_single_RM(list_hsd):
    single_RM_list=[]
    for hsd in list_hsd:
        single_RM_list.append(RM(hsd))
    return single_RM_list
  
def recomplete_RM(single_RM_list, list_hsd, prototype_RM_list):
    recompleted_RM=[]
    for RM in single_RM_list:
        recompleted_RM.append(RM.recomplete(list_hsd, prototype_RM_list))
    return recompleted_RM

def get_RM_prototype(list_RM):
    name_proto = []
    RM_proto = []
    for RM in list_RM:
        if RM.name not in name_proto:
            RM_proto.append(RM)
            name_proto.append(RM.name)
    return RM_proto
                
def get_unique_cluster(list_hsd):
    list_unique_hsd_name = []
    for hsd in list_hsd:
        if hsd.query_name not in list_unique_hsd_name :
            list_unique_hsd_name.append(hsd.query_name)
    return list_unique_hsd_name

def tri_hsd_by_strain_contig_midpos_evalue(hsd_list):
    return sorted(hsd_list, key=attrgetter("strain", "contig", "sbjct_start_end", "evalue"))

def extract_REBASE_hsd (REBASE_all_prot, unique_cluster):
    multifasta = ""
    REBASE_handle = open(REBASE_all_prot,"r")
    fasta_records = SeqIO.parse(REBASE_handle, "fasta")
    unique_cluster_handle = open(unique_cluster, "r")
    hsds = unique_cluster_handle.read().splitlines() 
    for fasta in fasta_records:
        if fasta.id in hsds:
            multifasta += fasta.format("fasta")
    REBASE_handle.close() 
    unique_cluster_handle.close()
    return multifasta

def parse_RM_tab_file(tab_file):
    handle = open(tab_file, "r")
    RMs = handle.readlines()
    listRM = []
    for RM in RMs :
        listRM.append(RM.split("\t"))
    return listRM
        #strain = line[0]
        #type = line[1]
        #Rmsus = line[2]
        #name = line[3]
        #complete = line[4]
        #size = line[5]
        #coor = line[6]
        #count = line[7]

def RM_to_itol(list_best_hsd_hit, list_RM_proto):
    list_RM_proto = sorted(list_RM_proto, key=attrgetter("type", "name", "count"))
    dict_strain = {}
    count = 0
    print len(list_best_hsd_hit)
    for hsd in list_best_hsd_hit:
        
        if not dict_strain.has_key(hsd.strain) :
            dict_strain[hsd.strain] = [hsd.query_name]
        else: 
            foo = dict_strain[hsd.strain]
            foo.append(hsd.query_name)
            dict_strain[hsd.strain] = foo
            
    itol_tab = "Strain" + "\t"

    for RM in list_RM_proto:
        itol_tab += str(RM.name).replace("[u'","").replace("']","") + "--COUNT = " + str(RM.count) + "\t"
        
    itol_tab += "\n"
    
    for strain in dict_strain:
        itol_tab += strain + "\t"
        for RM in list_RM_proto:
            if set(RM.name).issubset(set(dict_strain[strain])) or set(dict_strain[strain]).issubset(set(RM.name)):
                itol_tab += "1" + "\t"
                count += 1
            else :
                itol_tab += "0" + "\t"
        itol_tab += "\n"

    return (itol_tab, count)

def main():

    print "Parsing blast file ..."
    
    listhsd=parse_hsd_blast_xml("REBASE_r95_eval10-4_vs_all_Sepi.xml", min_match_length=200, min_coverage=0.5, min_identity=0, min_positive_fraction=0.5, max_evalue=0.0001)

    print "Sorting all hit ..."

    listhsd = tri_hsd_by_strain_contig_midpos_evalue(listhsd)
    all_hsd_hit = open("hsd_all.result", "w")
    all_hsd_hit.write(str(listhsd[0].print_attribut_name()) + "\n")
    for hsd in listhsd:
        all_hsd_hit.write(str(hsd.__str__()) + "\n")
    print str(len(listhsd)) + " hsd found"
    all_hsd_hit.close()


    print "Getting best hit ..."
    
    listbesthsd = get_best_hit_at_pos(listhsd)
    savebesthsd = listbesthsd
    result_best_hit_at_pos = open("hsd_best_hit.result","w")
    result_best_hit_at_pos.write(str(listbesthsd[0].print_attribut_name()) + "\n")
    for hsd in listbesthsd:
        result_best_hit_at_pos.write(str(hsd.__str__()) + "\n")
    print str(len(listbesthsd)) + " best hsd found"
    result_best_hit_at_pos.close()

    print "Searching complete RM system ..."

    completeRM = get_complete_RM(listbesthsd, 10000)
    result_RM = open("RM_complete.result","w")
    for RM in completeRM:
        RM.my_occurence_in_listRM(completeRM)
        result_RM.write(RM.strain + "\t" + RM.type + "\t" + str(RM.RMsus) + "\t" + str(RM.name) + "\t" + str(RM.complete) + "\t" + str(RM.size) + "\t" + str(RM.coor) + "\t" + str(RM.count) +"\n")
    print str(len(completeRM)) + " RM found"
    result_RM.close()

    print "Searching RM prototype ..."

    protoRM = get_RM_prototype(completeRM)
    RM_proto = open("RM_complete_prototype.result","w")
    for RM in protoRM: 
        RM_proto.write(RM.strain + "\t" + RM.type + "\t" + str(RM.RMsus) + "\t" + str(RM.name) + "\t" + str(RM.complete) + "\t" + str(RM.size) + "\t" + str(RM.coor) + "\t" + str(RM.count) +"\n")
    print str(len(protoRM)) + " RM found"
    RM_proto.close()


    print "Searching RM split in different contigs ..."

    remaining_hsd=get_hsd_not_in_RM(listbesthsd, completeRM)
    print str(len(listbesthsd)) + " best hsd found"
    single_RM = get_single_RM(remaining_hsd)
    recompleted_single_RM = [RM.recomplete(remaining_hsd, protoRM) for RM in single_RM]
    recomp_RM = [RM for RM in recompleted_single_RM if RM.recompleted]
    recomplete_and_complete_RM = completeRM + recomp_RM
    print str(len(recomp_RM)) + " splitted RM found"
    
    RM_complete_recomplete = open("RM_complete_recomplete.result","w")
    for RM in recomplete_and_complete_RM:
        RM.my_occurence_in_listRM(recomplete_and_complete_RM)
        RM_complete_recomplete.write(RM.strain + "\t" + RM.type + "\t" + str(RM.RMsus) + "\t" + str(RM.name) + "\t" + str(RM.complete) + "\t" + str(RM.size) + "\t" + str(RM.coor) + "\t" + str(RM.count) +"\n")
    RM_complete_recomplete.close()
    
    print "Updating uncomplete RM ..."

    remaining_hsd=get_hsd_not_in_RM(listbesthsd, recomplete_and_complete_RM)# Updating remaining hsd and RM
    single_RM = get_single_RM(remaining_hsd)
    print str(len(remaining_hsd)) + " isolated hsd remaining (not identified in RM)"

    RM_uncomplete = open("RM_uncomplete.result","w")
    for RM in single_RM :
        RM.my_occurence_in_listRM(single_RM)
        RM_uncomplete.write(RM.strain + "\t" + RM.type + "\t" + str(RM.RMsus) + "\t" + str(RM.name) + "\t" + str(RM.complete)+ "\t" + str(RM.size) + "\t" + str(RM.coor) + "\t" + str(RM.count) +"\n")
    RM_uncomplete.close()
    
    print "Counting occurence of RM prototype ..."

    RM_proto = open("RM_complete_prototype.result","w")
    for RM in protoRM:
        RM.my_occurence_in_listRM(completeRM)
        RM_proto.write(RM.strain + "\t" + RM.type + "\t" + str(RM.RMsus) + "\t" + str(RM.name) + "\t" + str(RM.complete) + "\t" + str(RM.size) + "\t" + str(RM.coor) + "\t" + str(RM.count) +"\n")
    RM_proto.close()

    print "Writing itol tabular file ..."
    
    itol_f = open("itol_tab.result", "w")
    itolout = RM_to_itol(savebesthsd, protoRM)
    itol_f.write(itolout[0])
    itol_f.close()
    print str(len(listbesthsd) - itolout[1]) + " unidentified hsd remaining"
 
main()



