'''
Created on 7 Sep 2020

@author: joshsinger
'''
import logging

from gluepy.aa_utils import AminoAcid, translate_region
from gluepy.ob_range import OneBasedRange


logging.basicConfig()
logger = logging.getLogger('genome_model')
# logger.setLevel(logging.INFO)

class GenomeModelException(Exception):
    pass

class GenomeModel:
    def check_almt_row(self, almt_row):
        if len(almt_row) != len(self.reference_nts):
            raise GenomeModelException("Alignment row length differs from genome model reference length")

    def get_region(self, region_name):
        return next(g_region for g_region in self.genome_regions if g_region.name == region_name)

    def get_regions(self):
        return self.genome_regions

    def extract_region(self, almt_row, region_name):
        self.check_almt_row(almt_row)
        g_region = self.get_region(region_name)
        return ''.join([ob_range.pull_from_str(almt_row) for ob_range in g_region.one_based_ranges])
    
    def __init__(self, reference_nts, genome_regions):
        self.reference_nts = reference_nts
        self.genome_regions = genome_regions
        # list of coding region names
        c_reg_names = [c_region.name for c_region in genome_regions if c_region.is_coding]
        
        # build dictionary mapping coding region name to reference NTs of that region
        self.ref_regions = {r_name: extracted_region for (r_name, extracted_region) in 
            map(lambda c_reg_name: (c_reg_name, self.extract_region(self.reference_nts, c_reg_name)),
                c_reg_names)}
        
        # build dictionary mapping coding region name to AA translation of that region
        self.ref_translations = {r_name: translation for (r_name, translation) in 
             map(lambda c_reg_name : (c_reg_name, AminoAcid.list_to_string(
                 translate_region(self.ref_regions[c_reg_name]))), 
                 c_reg_names)}


class GenomeRegion:
    def __init__(self, name, is_coding, one_based_ranges):
        self.name = name
        self.is_coding = is_coding
        self.one_based_ranges = one_based_ranges
        if is_coding:
            ob_length_sum = sum([ob_range.length() for ob_range in one_based_ranges])
            if ob_length_sum % 3 != 0:
                raise GenomeModelException("Total length of one-based ranges for region '"+name+"' is not a multiple of 3")
        
    def spanning_range(self):
        return OneBasedRange(min([ob_range.start for ob_range in self.one_based_ranges]),
                max([ob_range.end for ob_range in self.one_based_ranges]))
    

