'''
Created on 22 Sep 2020

@author: joshsinger
'''
from multiprocessing import Pool
import os
import shutil
from genomic_db.snps import snps_generator
from gluepy.fasta_utils import single_seq_pair_to_fasta_seq, \
    fasta_to_single_seq_gtor
import logging

logging.basicConfig()
logger = logging.getLogger('genomic_db')
logger.setLevel(logging.INFO)


num_workers = 16

def process_single_seq_pair(single_seq_pair):
    (header, nt_strings_list) = single_seq_pair
    fasta_seq = single_seq_pair_to_fasta_seq(header, nt_strings_list)
    return {
        "seq_id": fasta_seq.seq_id,
        "seq_snps": list(snps_generator(my_genome_model, "CSR", fasta_seq))
    }
    
def consume_seq_result(seq_result, snps_file):
    seq_id = seq_result["seq_id"]
    for (snp_id, is_ambig) in seq_result["seq_snps"]:
        snps_file.write("{},{},{}\n".format(seq_id,snp_id,is_ambig))
        
def init_worker(genome_model):
    global my_genome_model
    my_genome_model = genome_model



def build_db(almt_file_path, genome_model, output_dir):    
    with open(almt_file_path) as file_object:
        lines = file_object.readlines()
    
    shutil.rmtree(output_dir, ignore_errors=True)
    os.makedirs(output_dir)
    
    processed = 0
    with open(output_dir+"/snps.csv", "w") as snps_file:
        snps_file.write("seq_id,snp_id,is_ambig\n")
        with Pool(processes=num_workers,initializer=init_worker,initargs=(genome_model,)) as pool:
            for seq_result in pool.imap_unordered(process_single_seq_pair, fasta_to_single_seq_gtor(lines), 100):
                consume_seq_result(seq_result, snps_file)
                processed += 1
                if processed % 100 == 0:
                    snps_file.flush()
                    logger.info("Processed "+str(processed)+" sequences")
            logger.info("Processed "+str(processed)+" sequences")
            pool.close()
            pool.join()
