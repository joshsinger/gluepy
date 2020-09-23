'''
Created on 11 Sep 2020

@author: joshsinger
'''

import unittest

from cov.cov_genome_model import build_cov_genome_model
from gluepy.aa_utils import AminoAcid, translate_region
from gluepy.fasta_utils import fasta_file_to_dict


class Test(unittest.TestCase):
    def test1(self):
        genome_model = build_cov_genome_model()
        fasta_dict = fasta_file_to_dict("cog_alignment_small.fasta")
        for fasta_seq in fasta_dict.values():
            print(fasta_seq.seq_id+":"+AminoAcid.list_to_string(translate_region(genome_model.extract_region(fasta_seq.nt_chars, "nsp12"))))
