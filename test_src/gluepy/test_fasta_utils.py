'''
Created on 8 Sep 2020

@author: joshsinger
'''
import unittest

from gluepy.fasta_utils import fasta_file_to_dict


class Test(unittest.TestCase):


    def testFastaFileToDict1(self):
        fasta_dict = fasta_file_to_dict("testFastaFileToDict1.fasta")
        self.assertEquals(len(fasta_dict), 2)
        self.assertEquals(fasta_dict['Seq1'].description, 'desc1')
        self.assertEquals(fasta_dict['Seq1'].nt_chars, 'ACTGCTA----TCTGACC---NNTCTG')
        self.assertIsNone(fasta_dict['Seq2'].description)
        self.assertEquals(fasta_dict['Seq2'].nt_chars, 'ACTGCTATTTTTCRGACC---NNTCTG')


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()