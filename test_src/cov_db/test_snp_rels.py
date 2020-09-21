'''
Created on 21 Sep 2020

@author: joshsinger
'''
import unittest

from gluepy.genome_model import GenomeModel, GenomeRegion
from gluepy.ob_range import OneBasedRange
from gluepy.fasta_utils import FastaSequence
from cov_db.snp_rels import snp_csv_row_generator


class Test(unittest.TestCase):


    def testSnpRels(self):
        
        genome_model = GenomeModel("ATCGGCTGAGATGAGCCTAA", [
            GenomeRegion("csr", False, [OneBasedRange(5,20)])
        ])
        
        seq1 = FastaSequence("seq1", None, "ATCGGCTTAGATGAGCCTAA")
        self.assertEqual(["seq1,G8T,False"], list(snp_csv_row_generator(genome_model, "csr", seq1)))

        seq2 = FastaSequence("seq2", None, "ATCGGCTGAGATGAGCCTAA")
        self.assertEqual([], list(snp_csv_row_generator(genome_model, "csr", seq2)))


        seq3 = FastaSequence("seq2", None, "ATCGGCTGAGAYGAGCCTAA")
        self.assertEqual(['seq2,T12C,True'], list(snp_csv_row_generator(genome_model, "csr", seq3)))

        seq4 = FastaSequence("seq2", None, "ATCGGCNNNNNNGAGCCTAA")
        self.assertEqual([], list(snp_csv_row_generator(genome_model, "csr", seq4)))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSnpRels']
    unittest.main()