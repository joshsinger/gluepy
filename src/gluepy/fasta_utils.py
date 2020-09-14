'''
Created on 7 Sep 2020

@author: joshsinger
'''
import logging
from io import StringIO

logging.basicConfig()
logger = logging.getLogger('fasta_utils')
# logger.setLevel(logging.INFO)

class FastaException(Exception):
    pass

class FastaSequence:
    def __init__(self, seq_id, description, nt_chars):
        self.seq_id = seq_id
        self.nt_chars = nt_chars
        self.description = description



def fasta_file_to_dict(path):
    with open(path) as file_object:
        lines = file_object.readlines()
    return fasta_lines_to_dict(lines)

def fasta_bytes_to_dict(fasta_bytes):
    return fasta_lines_to_dict(fasta_bytes.decode("utf-8").splitlines())
        
def fasta_lines_to_dict(lines):
    fasta_dict = {} # dictionary mapping seq_id to FastaSequence
    # nt_char_sections contains the whitespace-free nucleotide lines of the current sequence
    seq_id, description, file_str = None, None, StringIO()
    line_number = 1
    seqs_added = 0
    
    for line in lines:
        if line.startswith('>'):
            if(not seq_id is None):
                add_fasta_sequence(fasta_dict, seq_id, description, file_str)
                seqs_added += 1
                if(seqs_added % 100 == 0):
                    logger.info("Sequences added "+str(seqs_added))
                seq_id, description, file_str = None, None, StringIO()
            first_spc_idx = line.find(' ')
            if(first_spc_idx == -1):
                seq_id = line[1:].rstrip()
            else:
                seq_id = line[1:first_spc_idx]
                description = line[first_spc_idx+1:].rstrip()
        else:
            add_nt_line(file_str, line, line_number)
        line_number += 1
    if(not seq_id is None):
        add_fasta_sequence(fasta_dict, seq_id, description, file_str)
        seqs_added += 1
    logger.info("Total sequences added "+str(seqs_added))
    return fasta_dict

nt_normalisation_switcher = {}

for nt_char in "ACGTRYKMSWBDHVN-":
    nt_normalisation_switcher[nt_char] = nt_char
nt_normalisation_switcher["?"] = "N"
for nt_char in "acgtrykmswbdhvn-":
    nt_normalisation_switcher[nt_char] = nt_char.upper()
for nt_char in "Uu":
    nt_normalisation_switcher[nt_char] = "T"
    
    
# normalise NT characters in line, capitalising, removing whitespace
# raise FastaException on illegal char
def add_nt_line(file_str, line, line_number):
    for char in line:
        if not char.isspace():
            norm_char = nt_normalisation_switcher[char]
            if(norm_char is None):
                raise FastaException("Illegal character '"+char+"' in FASTA line "+str(line_number));
            else:
                file_str.write(norm_char)
    
def add_fasta_sequence(fasta_dict, seq_id, description, file_str):
    if(seq_id in fasta_dict):
        raise FastaException("Duplicate seq_id "+seq_id)
    nt_chars = file_str.getvalue()
    fasta_dict[seq_id] = FastaSequence(seq_id, description, nt_chars)
    
