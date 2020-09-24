'''
Created on 23 Sep 2020

@author: joshsinger
'''
import unittest

from cov.cov_genome_model import build_cov_genome_model
from genomic_db.genomic_db import build_db


class Test(unittest.TestCase):


    def testBuildCovGenomicDb(self):
        almt_file_path = "/Users/joshsinger/coguk/cog_2020-09-07_all_alignment.fasta"
        output_dir = "/Users/joshsinger/coguk/graphdb_csvs"
        genome_model = build_cov_genome_model()
        snps_region = "CSR"
        num_workers = 16
        aa_regions = [c_region.name for c_region in genome_model.get_regions() if c_region.is_coding]
        seq_id_mapper = lambda old_id: old_id.split("/")[1]
        build_db(almt_file_path, seq_id_mapper, 
                 genome_model, snps_region, aa_regions, num_workers, output_dir)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testBuildCovGenomicDb']
    unittest.main()