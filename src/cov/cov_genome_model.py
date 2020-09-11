'''
Created on 11 Sep 2020

@author: joshsinger
'''
import pkgutil

from gluepy.fasta_utils import fasta_bytes_to_dict
from gluepy.genome_model import GenomeModel, GenomeRegion
from gluepy.ob_range import OneBasedRange


ref_fasta_seq = fasta_bytes_to_dict(pkgutil.get_data("cov", "EPI_ISL_402125.fasta"))["EPI_ISL_402125"]

cov_genome_model = GenomeModel(ref_fasta_seq.nt_chars, [
    GenomeRegion("nsp1", True, [OneBasedRange(266,805)]),
    GenomeRegion("nsp2", True, [OneBasedRange(806,2719)]),
    GenomeRegion("nsp3", True, [OneBasedRange(2720,8554)]),
    GenomeRegion("nsp4", True, [OneBasedRange(8555,10054)]),
    GenomeRegion("nsp5", True, [OneBasedRange(10055,10972)]),
    GenomeRegion("nsp6", True, [OneBasedRange(10973,11842)]),
    GenomeRegion("nsp7", True, [OneBasedRange(11843,12091)]),
    GenomeRegion("nsp8", True, [OneBasedRange(12092,12685)]),
    GenomeRegion("nsp9", True, [OneBasedRange(12686,13024)]),
    GenomeRegion("nsp10", True, [OneBasedRange(13025,13441)]),
    GenomeRegion("nsp11", True, [OneBasedRange(13442,13480)]),
    GenomeRegion("nsp12", True, [OneBasedRange(13442,13468),OneBasedRange(13468,16236)]),
    GenomeRegion("nsp13", True, [OneBasedRange(16237,18039)]),
    GenomeRegion("nsp14", True, [OneBasedRange(18040,19620)]),
    GenomeRegion("nsp15", True, [OneBasedRange(19621,20658)]),
    GenomeRegion("nsp16", True, [OneBasedRange(20659,21552)]),
    GenomeRegion("S", True, [OneBasedRange(21563,25384)]),
    GenomeRegion("ORF 3a", True, [OneBasedRange(25393,26220)]),
    GenomeRegion("E", True, [OneBasedRange(26245,26472)]),
    GenomeRegion("M", True, [OneBasedRange(26523,27191)]),
    GenomeRegion("ORF 6", True, [OneBasedRange(27202,27387)]),
    GenomeRegion("ORF 7a", True, [OneBasedRange(27394,27759)]),
    GenomeRegion("ORF 7b", True, [OneBasedRange(27756,27887)]),
    GenomeRegion("ORF 8", True, [OneBasedRange(27894,28259)]),
    GenomeRegion("N", True, [OneBasedRange(28274,29533)]),
    GenomeRegion("ORF 10", True, [OneBasedRange(29558,29674)]),
    GenomeRegion("CSR", False, [OneBasedRange(266,29674)])
])