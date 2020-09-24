'''
Created on 22 Sep 2020

@author: joshsinger
'''
from itertools import chain
import logging
from multiprocessing import Pool
import os
import shutil

from genomic_db.aa_reps import aa_replacements_generator
from genomic_db.snps import snps_generator
from gluepy.fasta_utils import single_seq_pair_to_fasta_seq, \
    fasta_to_single_seq_gtor


logging.basicConfig()
logger = logging.getLogger('genomic_db')
logger.setLevel(logging.INFO)


def process_single_seq_pair(single_seq_pair):
    (header, nt_strings_list) = single_seq_pair
    fasta_seq = single_seq_pair_to_fasta_seq(header, nt_strings_list)
    return {
        "seq_id": fasta_seq.seq_id,
        "seq_snps": list(snps_generator(my_genome_model, my_snps_region, fasta_seq)),
        "seq_aa_reps": list(chain.from_iterable(
            map(lambda aa_reg_name: aa_replacements_generator(my_genome_model, aa_reg_name, fasta_seq), my_aa_regions)))
    }
    
def consume_seq_result(context_mgr, seq_id_mapper, seq_result):
    seq_id = seq_id_mapper(seq_result["seq_id"])
    context_mgr.sequence_file.write("{},Sequence\n".format(seq_id))
    for (ref_nt_char, ref_nt_coord, qry_concrete_nt, is_ambig) in seq_result["seq_snps"]:
        snp_name = "{}{}{}".format(ref_nt_char, ref_nt_coord, qry_concrete_nt)
        if snp_name not in context_mgr.snp_names:
            context_mgr.snp_names.add(snp_name)
            context_mgr.snp_file.write("{},{},{},{},Snp\n".format(snp_name, ref_nt_char, ref_nt_coord, qry_concrete_nt))
        context_mgr.contains_snp_file.write("{},{},{},CONTAINS_SNP\n".format(seq_id, snp_name, is_ambig))
    for (region_name, ref_region_aa, codon_label, definite_aa, is_unique) in seq_result["seq_aa_reps"]:
        aa_rep_name = "{}:{}{}{}".format(region_name, ref_region_aa, codon_label, definite_aa)
        if aa_rep_name not in context_mgr.aa_rep_names:
            context_mgr.aa_rep_names.add(aa_rep_name)
            context_mgr.aa_replacement_file.write("{},{},{},{},AaReplacement\n".format(aa_rep_name, ref_region_aa, codon_label, definite_aa))
        context_mgr.contains_aa_rep_file.write("{},{},{},CONTAINS_AA_REP\n".format(seq_id,aa_rep_name,is_unique))
        
def init_worker(genome_model, snps_region, aa_regions):
    global my_genome_model
    my_genome_model = genome_model
    global my_snps_region
    my_snps_region = snps_region
    global my_aa_regions
    my_aa_regions = aa_regions

class ContextMgr():
    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.snp_names = set()
        self.aa_rep_names = set()
    def __enter__(self):
        self.sequence_file = open(self.output_dir+"/Sequence.csv", "w")
        self.sequence_file.write("sampleID:ID,:LABEL\n")
        self.snp_file = open(self.output_dir+"/Snp.csv", "w")
        self.snp_file.write("name:ID,refNucleotide,refCoordinate:int,snpNucleotide,:LABEL\n")
        self.aa_replacement_file = open(self.output_dir+"/AaReplacement.csv", "w")
        self.aa_replacement_file.write("name:ID,refAminoAcid,codonLabel:int,repAminoAcid,:LABEL\n")
        self.contains_snp_file = open(self.output_dir+"/CONTAINS_SNP.csv", "w")
        self.contains_snp_file.write(":START_ID,:END_ID,isAmbiguous:boolean,:TYPE\n")
        self.contains_aa_rep_file = open(self.output_dir+"/CONTAINS_AA_REP.csv", "w")
        self.contains_aa_rep_file.write(":START_ID,:END_ID,isUnique:boolean,:TYPE\n")
        return self
    def __exit__(self, _type, _value, _traceback):
        if self.sequence_file is not None:
            self.sequence_file.close()
        if self.snp_file is not None:
            self.snp_file.close()
        if self.aa_replacement_file is not None:
            self.aa_replacement_file.close()
        if self.contains_snp_file is not None:
            self.contains_snp_file.close()
        if self.contains_aa_rep_file is not None:
            self.contains_aa_rep_file.close()
    def flush_files(self):
        self.sequence_file.flush()
        self.snp_file.flush()
        self.aa_replacement_file.flush()
        self.contains_snp_file.flush()
        self.contains_aa_rep_file.flush()

def build_db(almt_file_path, seq_id_mapper, genome_model, snps_region, aa_regions, num_workers, output_dir):    
    with open(almt_file_path) as file_object:
        lines = file_object.readlines()
    
    shutil.rmtree(output_dir, ignore_errors=True)
    os.makedirs(output_dir)
    
    processed = 0
    with ContextMgr(output_dir) as context_mgr:
        with Pool(processes=num_workers,initializer=init_worker,initargs=(genome_model,snps_region,aa_regions)) as pool:
            for seq_result in pool.imap_unordered(process_single_seq_pair, fasta_to_single_seq_gtor(lines), 100):
                consume_seq_result(context_mgr, seq_id_mapper, seq_result)
                processed += 1
                if processed % 100 == 0:
                    context_mgr.flush_files()
                    logger.info("Processed "+str(processed)+" sequences")
            logger.info("Processed "+str(processed)+" sequences")
            pool.close()
            pool.join()
