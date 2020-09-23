'''
Created on 8 Sep 2020

@author: joshsinger
'''
import logging
from abc import ABC, abstractmethod
from gluepy.nt_utils import all_ambig_nts, ambig_nt_to_concrete_nts, concrete_nts_to_ambig_nt


logging.basicConfig()
logger = logging.getLogger('translation_utils')
# logger.setLevel(logging.INFO)


# standard text format for each possible codon table
# see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes
standard_codon_table = """
  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M------**--*----M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""
    
concrete_nt_triplet_to_aa = {}
aa_to_concrete_nt_triplets = {}

def init_from_codon_table(codon_table):
    # find the 4 relevant lines from the table
    lines = codon_table.splitlines()
    base1_line = next(line for line in lines if "Base1" in line)
    base1_line = base1_line[base1_line.rfind(" ")+1:]
    base2_line = next(line for line in lines if "Base2" in line)
    base2_line = base2_line[base2_line.rfind(" ")+1:]
    base3_line = next(line for line in lines if "Base3" in line)
    base3_line = base3_line[base3_line.rfind(" ")+1:]
    aas_line = next(line for line in lines if ("AAs" in line))
    aas_line = aas_line[aas_line.rfind(" ")+1:]
    # populate maps in both directions
    for i in range(0,64):
        nt_triplet = "".join([base1_line[i], base2_line[i], base3_line[i]])
        aa = aas_line[i]
        concrete_nt_triplet_to_aa[nt_triplet] = aa
        aa_to_concrete_nt_triplets.setdefault(aa, []).append(nt_triplet)

init_from_codon_table(standard_codon_table)

class TranslationException(Exception):
    pass

class AminoAcid(ABC):
    def __init__(self, codon_label):
        self.codon_label = codon_label
    @abstractmethod
    def to_aa_char(self):
        pass
    @staticmethod
    def list_to_string(amino_acid_list):
        return ''.join([amino_acid.to_aa_char() for amino_acid in amino_acid_list])

class DeletedAminoAcid(AminoAcid):
    def __init__(self, codon_label, am_nts):
        super().__init__(codon_label)
        self.am_nts = am_nts
    def to_aa_char(self):
        return "-"


class LabeledAminoAcid(AminoAcid):
    def __init__(self, codon_label, am_tr_result):
        super().__init__(codon_label)
        self.am_tr_result = am_tr_result
    def to_aa_char(self):
        return self.am_tr_result.aa



# Encapsulates possible / definite translations, based on a 
# nucleotide triplet
# which may contain ambiguity codes.

class AmbigTripletResult:
    def __init__(self, ambig_triplet, aa, c_triplets, definite_aas, possible_aas):
        self.ambig_triplet = ambig_triplet
        self.concrete_triplets = c_triplets
        self.aa = aa
        self.definite_aas = definite_aas
        self.possible_aas = possible_aas

# find the set of all possible concrete NT triplets for the ambiguous NT triplet
# mapping each of these to AA, this gives the set of possible AAs for the ambiguous NT triplet
# for each AA in turn
#    - delete the concrete NT triplets which code for that AA from the concrete triplets set.
#    - if the set is now empty then the AA is "definitely present".
#      otherwise
#    - recompute the ambiguous triplet for the remaining concrete triplets.
#    - if this is different from the original ambiguous triplet, 
#            have a look at the position(s) where the recomputed triplet differs from the original triplet
#            If there is an N in the original triplet at any of these positions, but a non-N in the recomputed triplet
#            then the AA is not "definitely present".
#            otherwise it is "definitely present", 
#      otherwise it is merely "possibly present".
        
def find_concrete_triplets(concretes1, concretes2, concretes3):
    c_triplets = []
    for concrete1 in concretes1:
        for concrete2 in concretes2:
            for concrete3 in concretes3:
                c_triplets.append(''.join([concrete1, concrete2, concrete3]))      
    return c_triplets

def unique_sorted(a_list):
    result = list(set(a_list))
    result.sort()
    return result

def result_for_ambig_triplet(ambig_triplet):
    #logger.info("ambig_triplet:"+ambig_triplet)
    concretes1 = ambig_nt_to_concrete_nts[ambig_triplet[0]]
    concretes2 = ambig_nt_to_concrete_nts[ambig_triplet[1]]
    concretes3 = ambig_nt_to_concrete_nts[ambig_triplet[2]]
    c_triplets = unique_sorted(find_concrete_triplets(concretes1, concretes2, concretes3))
    #logger.info("c_triplets:"+str(c_triplets))
    possible_aas = unique_sorted([concrete_nt_triplet_to_aa[triplet] for triplet in c_triplets])
    #logger.info("possible_aas:"+str(possible_aas))
    definite_aas = []
    for possible_aa in possible_aas:
        is_definite = False
        coding_triplets = aa_to_concrete_nt_triplets[possible_aa]
        #logger.info("coding_triplets for "+possible_aa+":"+str(coding_triplets))
        remaining_c_triplets = [triplet for triplet in c_triplets if triplet not in coding_triplets]
        #logger.info("remaining_c_triplets for "+possible_aa+":"+str(remaining_c_triplets))
        if len(remaining_c_triplets) == 0:
            is_definite = True
        else:
            remaining_ctrip_ambigs = []
            for i in range(3):
                concretes_at_i = [r_c_triplet[i] for r_c_triplet in remaining_c_triplets]
                ambig_at_i = concrete_nts_to_ambig_nt[''.join(unique_sorted(concretes_at_i))]
                remaining_ctrip_ambigs.append(ambig_at_i)
            remaining_ambig_triplet = ''.join(remaining_ctrip_ambigs)
            #logger.info("remaining_ambig_triplet for "+possible_aa+":"+remaining_ambig_triplet)
            if remaining_ambig_triplet != ambig_triplet:
                is_definite = True
                for i in range(3):
                    if ambig_triplet[i] == 'N' and remaining_ambig_triplet[i] != 'N':
                        is_definite = False
            #logger.info(possible_aa+" is_definite:"+str(is_definite))
        if is_definite:
            definite_aas.append(possible_aa)                    
    if len(possible_aas) == 1 and len(definite_aas) == 1:
        aa = possible_aas[0]
    else:
        aa= 'X'
    return AmbigTripletResult(ambig_triplet, aa, c_triplets, ''.join(definite_aas), ''.join(possible_aas))

amb_tripl_to_result = {}

for ambig1 in all_ambig_nts:
    for ambig2 in all_ambig_nts:
        for ambig3 in all_ambig_nts:
            ambig_triplet = ''.join([ambig1, ambig2, ambig3])
            amb_tripl_to_result[ambig_triplet] = result_for_ambig_triplet(ambig_triplet)





def translate_region(almt_row_region):
    region_length = len(almt_row_region)
    if region_length % 3 != 0:
        raise TranslationException("Attempt to translate region with length not multiple of 3")
    result = []
    nt_coord = 0
    codon_label_int = 1
    while nt_coord < len(almt_row_region):
        codon_label = str(codon_label_int)
        am_nts = almt_row_region[nt_coord:nt_coord+3]
        if am_nts.find("-") != -1:
            result.append(DeletedAminoAcid(codon_label, am_nts))
        else:
            result.append(LabeledAminoAcid(codon_label, amb_tripl_to_result[am_nts]))
        nt_coord += 3
    return result
    