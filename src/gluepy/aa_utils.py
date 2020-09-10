'''
Created on 8 Sep 2020

@author: joshsinger
'''
# standard text format for each possible codon table
# see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes
standard_codon_table = """
  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M------**--*----M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""
	
concrete_nt_triplet_to_aa = {}
aa_to_concrete_nt_triplets = {}

def init_from_codon_table(codon_table):
	# find the 4 relevant lines from the table
	lines = codon_table.splitlines()
	base1_line = next(line for line in lines if "Base1" in line)
	base1_line = base1_line[base1_line.rfind(" ")+1:]
	base2_line = next(line for line in lines if "Base2" in line)
	base2_line = base2_line[base2_line.rfind(" ")+1:]
	base3_line = next(line for line in lines if "Base3" in line)
	base3_line = base3_line[base3_line.rfind(" ")+1:]
	aas_line = next(line for line in lines if ("AAs" in line))
	aas_line = aas_line[aas_line.rfind(" ")+1:]
	# populate maps in both directions
	for i in range(0,64):
		nt_triplet = "".join([base1_line[i], base2_line[i], base3_line[i]])
		aa = aas_line[i]
		concrete_nt_triplet_to_aa[nt_triplet] = aa
		aa_to_concrete_nt_triplets.setdefault(aa, []).append(nt_triplet)

init_from_codon_table(standard_codon_table)
