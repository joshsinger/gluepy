'''
Created on 23 Sep 2020

@author: joshsinger
'''


'''
Given a genome model and a coding region name, return a generator of tuples representing definite AA replacements
e.g. (S, D, 614, G, True)
which would be S:D614G plus a boolean indicating whether the AA replacement is the unique definite AA
at that location. AA replacements are relative to the reference sequence
'''

from builtins import map
from itertools import chain

from gluepy.aa_utils import amb_tripl_to_result


def aa_replacements_generator(genome_model, region_name, fasta_sequence):
    qry_nt_chars = fasta_sequence.nt_chars
    qry_region_nts = genome_model.extract_region(qry_nt_chars, region_name)
    ref_region_nts = genome_model.ref_regions[region_name]
    ref_region_aas = genome_model.ref_translations[region_name]

    # codon index
    codon_idx_range = range(0, len(ref_region_nts)//3)
    # generator of generators
    # transform each codon index to a generator of pairs (possibly 0, 1 or more)
    gtors_gtor = map( lambda i: aa_reps_at_coord(region_name, ref_region_nts, ref_region_aas[i], qry_region_nts, i*3, str(i+1)), codon_idx_range)
    return chain.from_iterable(gtors_gtor) # chain the generators together into one big generator
    
# generator yielding 0 or more pairs for a specific codon    
# uses a simple optimisation where if the nucleotides are the same, no replacement is present
def aa_reps_at_coord(region_name, ref_region_nts, ref_region_aa, qry_region_nts, nt_start_coord, codon_label):
    ref_triplet = ref_region_nts[nt_start_coord:nt_start_coord+3]
    qry_triplet = qry_region_nts[nt_start_coord:nt_start_coord+3]
    if qry_triplet == ref_triplet:
        return
    if qry_triplet.find("-") != -1:
        return
    qry_ambig_triplet_result = amb_tripl_to_result[qry_triplet]
    is_unique = ( len(qry_ambig_triplet_result.definite_aas) == 1 )
    for definite_aa in qry_ambig_triplet_result.definite_aas:
        if definite_aa != ref_region_aa:
            yield (region_name, ref_region_aa, codon_label, definite_aa, is_unique)