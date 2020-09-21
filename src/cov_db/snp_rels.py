'''
Created on 21 Sep 2020

@author: joshsinger
'''
'''
Given a genome model and an alignment row, return a generator of Neo4J import CSV row strings for 
the relationship between the sequence and its SNPs relative to the reference sequence within region_name
'''

from builtins import map
from itertools import chain

from gluepy import genome_model
from gluepy.fasta_utils import ob_nt_char
from gluepy.nt_utils import ambig_nt_to_concrete_nts


def snp_csv_row_generator(genome_model, region_name, fasta_sequence):
    seq_id = fasta_sequence.seq_id
    qry_nt_chars = fasta_sequence.nt_chars
    ref_nt_chars = genome_model.reference_nts
    # ob_range for the named region 
    spanning_range = genome_model.get_region(region_name).spanning_range()
    # generator of generators
    # transform each ref coord to a generator of csv rows (possibly 0, 1 or more)
    gtors_gtor = map( lambda i: snp_csv_rows(seq_id, ref_nt_chars, qry_nt_chars, i), range(spanning_range.start, spanning_range.end+1))
    #gtors_gtor = map( lambda i: (yield i), range(spanning_range.start, spanning_range.end+1))
    return chain.from_iterable(gtors_gtor) # chain the generators together into one big generator
    
    
def snp_csv_rows(seq_id, ref_nt_chars, qry_nt_chars, ref_nt_coord):
    qry_nt_char = ob_nt_char(qry_nt_chars, ref_nt_coord)  
    ref_nt_char = ob_nt_char(ref_nt_chars, ref_nt_coord)
    if qry_nt_char != '-' and qry_nt_char != 'N' and qry_nt_char != ref_nt_char:
        qry_concrete_nts = ambig_nt_to_concrete_nts[qry_nt_char]
        for qry_concrete_nt in qry_concrete_nts:
            if qry_concrete_nt != ref_nt_char:
                snp_id = ref_nt_char + str(ref_nt_coord) + qry_concrete_nt
                is_ambig = len(qry_concrete_nts) > 1
                yield seq_id + "," + snp_id + "," + str(is_ambig)