'''
Created on 14 Sep 2020

@author: joshsinger
'''
import time
import unittest
from gluepy import fasta_utils


class Test(unittest.TestCase):


    def testSync(self):
        t0 = time.time()
        path = "/Users/joshsinger/coguk/cog_2020-09-07_all_alignment.fasta"
        with open(path) as file_object:
            lines = file_object.readlines()
        t1 = time.time() - t0
        print("Time elapsed: ", t1) 
        print("Num lines: ", len(lines)) 

    def testLoadFasta(self):
        t0 = time.time()
        path = "/Users/joshsinger/coguk/cog_2020-09-07_all_alignment.fasta"
        fasta_dict = fasta_utils.fasta_file_to_dict(path)
        t1 = time.time() - t0
        print("Time elapsed: ", t1) 
        print("Num sequences: ", len(fasta_dict)) 


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()