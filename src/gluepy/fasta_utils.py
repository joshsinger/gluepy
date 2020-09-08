'''
Created on 7 Sep 2020

@author: joshsinger
'''
import re

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
        
def fasta_lines_to_dict(lines):
    fasta_dict = {} # dictionary mapping seq_id to FastaSequence
    # nt_char_sections contains the whitespace-free nucleotide lines of the current sequence
    seq_id, description, nt_char_list = None, None, []
    line_number = 1;
    
    for line in lines:
        if line.startswith('>'):
            if(not seq_id is None):
                add_fasta_sequence(fasta_dict, seq_id, description, nt_char_list)
                seq_id, description, nt_char_list = None, None, []
            first_spc_idx = line.find(' ')
            if(first_spc_idx == -1):
                seq_id = line[1:].rstrip()
            else:
                seq_id = line[1:first_spc_idx]
                description = line[first_spc_idx+1:].rstrip()
        else:
            add_nt_line(nt_char_list, line, line_number)
        line_number += 1
    if(not seq_id is None):
        add_fasta_sequence(fasta_dict, seq_id, description, nt_char_list)
    return fasta_dict

# normalise NT characters in line, capitalising, removing whitespace
# raise FastaException on illegal char
def add_nt_line(nt_char_list, line, line_number):
    for char in line:
        if "ACGTRYKMSWBDHVN-".find(char) != -1:
            nt_char_list.append(char)
        elif "acgtrykmswbdhvn".find(char) != -1:
            nt_char_list.append(char.upper())
        elif "Uu".find(char) != -1:
            nt_char_list.append("T")    
        elif char.isspace():
            pass
        else:
            raise FastaException("Illegal character '"+char+"' in FASTA line "+line_number);
    
def add_fasta_sequence(fasta_dict, seq_id, description, nt_char_list):
    if(seq_id in fasta_dict):
        raise FastaException("Duplicate seq_id "+seq_id)
    nt_chars = ''.join(nt_char_list)
    fasta_dict[seq_id] = FastaSequence(seq_id, description, nt_chars)
    
