'''
Created on 10 Sep 2020

@author: joshsinger
'''

from gluepy.nt_utils import all_ambig_nts, ambig_nt_to_concrete_nts
from gluepy.aa_utils import aa_to_concrete_nt_triplets, concrete_nt_triplet_to_aa

# Encapsulates possible / definite translations, based on a 
# nucleotide triplet
# which may contain ambiguity codes.

class AmbigTripletResult:
    def __init__(self, ambig_triplet, aa, concrete_triplets, definite_aas, possible_aas):
        self.ambig_triplet = ambig_triplet
        self.concrete_triplets = concrete_triplets
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
    concrete_triplets = []
    for concrete1 in concretes1:
        for concrete2 in concretes2:
            for concrete3 in concretes3:
                concrete_triplets.append(''.join([concrete1, concrete2, concrete3]))      
    return concrete_triplets

def unique_sorted(a_list):
    result = list(set(a_list))
    result.sort()
    return result


ambig_triplet_to_result = {}

for ambig1 in all_ambig_nts:
    concretes1 = ambig_nt_to_concrete_nts[ambig1]
    for ambig2 in all_ambig_nts:
        concretes2 = ambig_nt_to_concrete_nts[ambig2]
        for ambig3 in all_ambig_nts:
            concretes3 = ambig_nt_to_concrete_nts[ambig3]
            
            ambig_triplet = ''.join([ambig1, ambig2, ambig3])
            concrete_triplets = find_concrete_triplets(concretes1, concretes2, concretes3)
            possible_aas = ''.join(unique_sorted([concrete_nt_triplet_to_aa[triplet] for triplet in concrete_triplets]))
            if len(possible_aas) == 1:
                aa = possible_aas[0]
            else:
                aa= 'X'

            ambig_triplet_to_result[ambig_triplet] = AmbigTripletResult(ambig_triplet, aa, concrete_triplets, '', possible_aas)
