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

# generator which takes a list of lines constituting a multi-FASTA file 
# the generator produces pairs, each pair contains the header string 
# and a list of nucleotide lines for a single FASTA sequence.
def fasta_to_single_seq_gtor(lines):
    header, nt_strings_list = None, []
    for line in lines:
        line = line.rstrip();
        if len(line) > 0:
            if line.startswith('>'):
                if(header is None):
                    header = line
                else:
                    yield (header, nt_strings_list)
                    header, nt_strings_list = line, []
            else:
                nt_strings_list.append(line)
    if(not header is None):
        yield (header, nt_strings_list)

def single_seq_pair_to_fasta_seq(header, nt_strings_list):
    first_spc_idx = header.find(' ')
    if(first_spc_idx == -1):
        seq_id = header[1:].rstrip()
        description = None
    else:
        seq_id = header[1:first_spc_idx]
        description = header[first_spc_idx+1:].rstrip()
    file_str_io = StringIO()
    for nt_string in nt_strings_list:
        add_nt_line(file_str_io, nt_string)
    return FastaSequence(seq_id, description, file_str_io.getvalue())


def fasta_bytes_to_dict(fasta_bytes):
    return fasta_lines_to_dict(fasta_bytes.decode("utf-8").splitlines())
        
def fasta_lines_to_dict(lines):
    fasta_dict = {} # dictionary mapping seq_id to FastaSequence
    sequences_processed = 0
    for (header, nt_strings_list) in fasta_to_single_seq_gtor(lines):
        fasta_seq = single_seq_pair_to_fasta_seq(header, nt_strings_list)
        seq_id = fasta_seq.seq_id
        if(seq_id in fasta_dict):
            raise FastaException("Duplicate seq_id "+seq_id)
        fasta_dict[seq_id] = fasta_seq
        sequences_processed += 1
        if sequences_processed % 100 == 0:
            logger.info("Processed "+str(sequences_processed)+" sequences")
    logger.info("Processed "+str(sequences_processed)+" sequences")
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
def add_nt_line(file_str, line):
    for char in line:
        if not char.isspace():
            norm_char = nt_normalisation_switcher[char]
            if(norm_char is None):
                raise FastaException("Illegal character '"+char+"' in FASTA");
            else:
                file_str.write(norm_char)
    

    
# pull a single nucleotide character from a string, based on a one-based coordinate
def ob_nt_char(nt_chars, ob_coord):
    return nt_chars[ob_coord-1]

