'''
Created on 21 Sep 2020

@author: joshsinger
'''

from builtins import map
from itertools import chain

from gluepy.fasta_utils import ob_nt_char
from gluepy.nt_utils import ambig_nt_to_concrete_nts

'''
Given a genome model and an alignment row, return a generator of pairs,
each pair being (<SNP_id>, <is_ambig>)
where <SNP_id> is e.g. C14805T and <is_ambig> is a boolean indicating whether the SNP was generated as a result of a 
non-N ambiguity code
The SNPs are relative to the reference sequence and region_name
'''
def snps_generator(genome_model, region_name, fasta_sequence):
    qry_nt_chars = fasta_sequence.nt_chars
    ref_nt_chars = genome_model.reference_nts
    # ob_range for the named region 
    spanning_range = genome_model.get_region(region_name).spanning_range()
    # generator of generators
    # transform each ref coord to a generator of csv rows (possibly 0, 1 or more)
    gtors_gtor = map( lambda i: snps_at_coord(ref_nt_chars, qry_nt_chars, i), range(spanning_range.start, spanning_range.end+1))
    #gtors_gtor = map( lambda i: (yield i), range(spanning_range.start, spanning_range.end+1))
    return chain.from_iterable(gtors_gtor) # chain the generators together into one big generator
    
    
def snps_at_coord(ref_nt_chars, qry_nt_chars, ref_nt_coord):
    qry_nt_char = ob_nt_char(qry_nt_chars, ref_nt_coord)  
    ref_nt_char = ob_nt_char(ref_nt_chars, ref_nt_coord)
    if qry_nt_char != '-' and qry_nt_char != 'N' and qry_nt_char != ref_nt_char:
        qry_concrete_nts = ambig_nt_to_concrete_nts[qry_nt_char]
        for qry_concrete_nt in qry_concrete_nts:
            if qry_concrete_nt != ref_nt_char:
                snp_id = ref_nt_char + str(ref_nt_coord) + qry_concrete_nt
                is_ambig = len(qry_concrete_nts) > 1
                yield (snp_id,is_ambig)