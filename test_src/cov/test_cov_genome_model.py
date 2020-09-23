'''
Created on 11 Sep 2020

@author: joshsinger
'''

import unittest

from cov.cov_genome_model import build_cov_genome_model


class Test(unittest.TestCase):
    def test1(self):
        genome_model = build_cov_genome_model()
        self.assertEqual("SADAQSFLNGFAV", genome_model.ref_translations.get("nsp11"))