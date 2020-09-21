'''
Created on 7 Sep 2020

@author: joshsinger
'''
import logging
from gluepy.aa_utils import DeletedAminoAcid, LabeledAminoAcid
from gluepy.translation import result_for_ambig_triplet
from gluepy.ob_range import OneBasedRange

logging.basicConfig()
logger = logging.getLogger('genome_model')
# logger.setLevel(logging.INFO)

class GenomeModelException(Exception):
    pass

class GenomeModel:
    def __init__(self, reference_nts, genome_regions):
        self.reference_nts = reference_nts
        self.genome_regions = genome_regions

    def check_almt_row(self, almt_row):
        if len(almt_row) != len(self.reference_nts):
            raise GenomeModelException("Alignment row length differs from genome model reference length")

    def get_region(self, region_name):
        return next(g_region for g_region in self.genome_regions if g_region.name == region_name)

    def extract_region(self, almt_row, region_name):
        self.check_almt_row(almt_row)
        g_region = self.get_region(region_name)
        return ''.join([ob_range.pull_from_str(almt_row) for ob_range in g_region.one_based_ranges])
        
    def translate_region(self, almt_row, region_name):
        g_region = self.get_region(region_name)
        if not g_region.is_coding:
            raise GenomeModelException("Cannot translate non-coding region "+region_name)
        extracted_region = self.extract_region(almt_row, region_name)
        result = []
        i = 0
        while i < len(extracted_region):
            codon_label = str(i)
            am_nts = extracted_region[i:i+3]
            if am_nts.find("-") != -1:
                result.append(DeletedAminoAcid(region_name, codon_label, am_nts))
            else:
                result.append(LabeledAminoAcid(region_name, codon_label, result_for_ambig_triplet(am_nts)))
            i += 3
        return result

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
    

