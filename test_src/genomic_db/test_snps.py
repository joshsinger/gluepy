'''
Created on 21 Sep 2020

@author: joshsinger
'''
import unittest

from gluepy.genome_model import GenomeModel, GenomeRegion
from gluepy.ob_range import OneBasedRange
from gluepy.fasta_utils import FastaSequence
from genomic_db.snps import snps_generator


class Test(unittest.TestCase):


    def testSnpRels(self):
        
        genome_model = GenomeModel("ATCGGCTGAGATGAGCCTAA", [
            GenomeRegion("csr", False, [OneBasedRange(5,20)])
        ])
        
        seq1 = FastaSequence("seq1", None, "ATCGGCTTAGATGAGCCTAA")
        self.assertEqual([("G8T",False)], list(snps_generator(genome_model, "csr", seq1)))

        seq2 = FastaSequence("seq2", None, "ATCGGCTGAGATGAGCCTAA")
        self.assertEqual([], list(snps_generator(genome_model, "csr", seq2)))


        seq3 = FastaSequence("seq2", None, "ATCGGCTGAGAYGAGCCTAA")
        self.assertEqual([("T12C",True)], list(snps_generator(genome_model, "csr", seq3)))

        seq4 = FastaSequence("seq2", None, "ATCGGCNNNNNNGAGCCTAA")
        self.assertEqual([], list(snps_generator(genome_model, "csr", seq4)))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSnpRels']
    unittest.main()