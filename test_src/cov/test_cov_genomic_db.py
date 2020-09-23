'''
Created on 23 Sep 2020

@author: joshsinger
'''
import pickle
import unittest

from cov.cov_genome_model import build_cov_genome_model
from genomic_db.genomic_db import build_db


class Test(unittest.TestCase):


    def testBuildCovGenomicDb(self):
        almt_file_path = "/Users/joshsinger/coguk/cog_2020-09-07_all_alignment.fasta"
        output_dir = "/Users/joshsinger/coguk/graphdb_csvs"
        genome_model = build_cov_genome_model()
        
        pickle.dumps(genome_model)
        build_db(almt_file_path, genome_model, output_dir)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testBuildCovGenomicDb']
    unittest.main()