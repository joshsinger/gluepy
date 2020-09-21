'''
Created on 11 Sep 2020

@author: joshsinger
'''

import unittest

from gluepy.fasta_utils import fasta_file_to_dict
from cov.cov_genome_model import cov_genome_model
from gluepy.aa_utils import AminoAcid


class Test(unittest.TestCase):
    def test1(self):
        fasta_dict = fasta_file_to_dict("cog_alignment_small.fasta")
        for fasta_seq in fasta_dict.values():
            print(fasta_seq.seq_id+":"+AminoAcid.list_to_string(cov_genome_model.translate_region(fasta_seq.nt_chars, "nsp12")))