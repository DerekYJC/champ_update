from champ import link
import logging
import os
import pickle
from champ.link import FastqFiles

log = logging.getLogger(__name__)


#in controller/
            
def main(clargs):
    min_len = clargs.min_len
    max_len = clargs.max_len
    fastq_filenames = [os.path.join(clargs.fastq_directory, directory) for directory in os.listdir(clargs.fastq_directory)]
    layers = ['floor', 'ceiling']
    for i in range(2):
        log.debug('Constructing the all seq name files for {} layer'.format(layers[i]))   
        usable_read = lambda record_id: link.determine_side(record_id) == str(i+1)
        output_directory = clargs.output_directory
        
        fastq_files = FastqFiles(fastq_filenames)
        read_names_given_seq = {}
        
        with open(clargs.log_p_file_path) as f:
            log_p_struct = pickle.load(f)
        
        if not os.path.exists(os.path.join(output_directory, 'read_names_of_all_seq_{}_{}.txt'.format(clargs.chip_id, layers[i]))):
            read_names_given_seq = link.determine_sequences_of_read_names(log_p_struct, fastq_files, usable_read, min_len, max_len)
            link.write_read_names_by_sequence(read_names_given_seq, os.path.join(output_directory, 'read_names_of_all_seq_{}_{}.txt'.format(clargs.chip_id, layers[i])))
            log.debug('Writing all seq name of {} layer for {} chip!'.format(layers[i], clargs.chip_id))