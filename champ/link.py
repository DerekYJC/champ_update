from Bio import SeqIO
from champ.adapters_cython import simple_hamming_distance
from collections import defaultdict
import editdistance
import gzip
import itertools
import logging
import numpy as np
import os
import pickle
import pysam
import random
import subprocess
import yaml
import re

log = logging.getLogger(__name__)

#in champ/

class FastqFiles(object):
    """ Sorts compressed FastQ files provided to us from the Illumina sequencer. """
    def __init__(self, filenames):
        self._filenames = list(self._filter_names(filenames))

    def __iter__(self):
        for f in self._filenames:
            yield f

    def __len__(self):
        return len(self._filenames)

    @property
    def alignment_length(self):
        paired_length = len([(f1, f2) for f1, f2 in self.paired])
        single_length = len([f for f in self.single])
        return paired_length + single_length

    @property
    def paired(self):
        for f1, f2 in self._sort_filenames(paired=True):
            yield f1, f2

    @property
    def single(self):
        for f in self._sort_filenames(paired=False):
            yield f

    def _filter_names(self, data):
        # eliminate filenames that can't possibly be fastq files of interest
        for filename in reversed(data):
            if not filename.endswith('fastq.gz'):
                continue
            if '_I1_' in filename or '_I2_' in filename or '_I1.' in filename or '_I2.' in filename:
                continue
            yield filename

    def _sort_filenames(self, paired=True):
        # yield filenames that are the given type (single or paired)
        for filename in self._filenames:
            if '_R1_' in filename or '_R1.' in filename:
                pair = filename.replace('_R1_', '_R2_').replace('_R1.', '_R2.')
                if paired and pair in self._filenames:
                    yield filename, pair
                elif not paired and pair not in self._filenames:
                    yield filename


class FastqReadClassifier(object):
    def __init__(self, bowtie_path):
        clean_path = bowtie_path.rstrip(os.path.sep)
        self.name = os.path.basename(clean_path)
        self._common_command = ('bowtie2', '--local', '-p 15', '--no-unal', '-x %s' % clean_path)

    def paired_call(self, fastq_file_1, fastq_file_2):
        command = self._common_command + ('-1 ' + fastq_file_1,
                                          '-2 ' + fastq_file_2,
                                          '-S chimp.sam',
                                          '2>&1 | tee error.txt')
        return self._run(command)

    def single_call(self, fastq_file):
        command = self._common_command + ('-U ' + fastq_file,)
        return self._run(command)

    def _run(self, command):
        with open('/dev/null', 'w+') as devnull:
            shell_options = dict(shell=True, stderr=devnull, stdout=devnull)
            subprocess.call(' '.join(command), **shell_options)
            sam_command = 'samtools view -bS chimp.sam | samtools sort - final'
            subprocess.call(sam_command, **shell_options)
            subprocess.call('samtools index final.bam', **shell_options)
            for r in pysam.Samfile('final.bam'):
                yield r.qname
        for temp_file in ('chimp.sam', 'final.bam', 'error.txt', 'final.bam.bai'):
            try:
                os.unlink(temp_file)
            except (OSError, IOError):
                log.warn("Unable to delete temp file: %s. "
                         "Was it not created? You may be missing FASTQ reads." % temp_file)

def get_max_ham_dists(min_len, max_len):
    dists = defaultdict(list)
    for _ in xrange(50000):
        ref_seq = rand_seq(max_len)
        new_seq = rand_seq(max_len)
        for i in range(min_len, max_len+1):
            dists[i].append(simple_hamming_distance(ref_seq[:i], new_seq[:i]))
    max_ham_dists = [min(np.percentile(dists[i], 0.1), int(i/4)) for i in range(min_len, max_len+1)]
    return max_ham_dists

def determine_sequences_of_read_names(log_p_struct, fastq_files, usable_read, min_len, max_len):
    # --------------------------------------------------------------------------------
    # Pair fpaths and classify seqs
    # --------------------------------------------------------------------------------
    max_ham_dists = get_max_ham_dists(min_len, max_len)
#    log.debug("Max ham dists: %s" % str(max_ham_dists))
    read_names_given_seq = defaultdict(list)
    for fpath1, fpath2 in fastq_files.paired:
        log.debug('{}, {}'.format(*map(os.path.basename, (fpath1, fpath2))))
        discarded = 0
        total = 0
        for i, (rec1, rec2) in enumerate(
                itertools.izip(parse_fastq_lines(fpath1),
                               parse_fastq_lines(fpath2))
        ):
            if not usable_read(rec1.id):
                continue
            total += 1
            seq = classify_seq(rec1, rec2, min_len, max_len, max_ham_dists, log_p_struct)
            if seq:
                read_names_given_seq[seq].append(str(rec1.id))
            else:
                discarded += 1
        found = total - discarded
        print('Found {} of {} ({:.1f}%)'.format(found, total, 100 * found / float(total)))
    return read_names_given_seq
    
def determine_side(record_id):
    """ 
    DNA is sequenced on both sides of the chip, however the TIRF microscope can only see one side, so we want to 
    be able to ignore reads that we can't see just to save time/memory. 
    
    """
    return record_id.split(":")[4][0]
    
def classify_seq(rec1, rec2, min_len, max_len, max_ham_dists, log_p_struct):
    bases = set('ACGT')
    ML_bases = []
    # Store as strings
    seq1 = str(rec1.seq)
    seq2_rc = str(rec2.seq.reverse_complement())
    loc_max_len = min(max_len, len(seq1), len(seq2_rc))

    # Find aligning sequence, indels are not allowed, starts of reads included
    
    sig_lens = [i for i, max_ham in zip(range(min_len, loc_max_len + 1), max_ham_dists)
                if simple_hamming_distance(seq1[:i], seq2_rc[-i:]) < max_ham]
    if len(sig_lens) != 1:
        return ''.join(seq1[min_len:max_len])

    seq2_len = sig_lens[0]
    seq2_match = seq2_rc[-seq2_len:]
    seq1_match = seq1[:seq2_len]

    # Get corresponding quality scores
    quals1 = rec1.letter_annotations['phred_quality'][:seq2_len]
    quals2 = rec2.letter_annotations['phred_quality'][::-1][-seq2_len:]
 #   print quals1, quals2

    # Build consensus sequence
    #ML_bases = []
    for r1, q1, r2, q2 in zip(seq1_match, quals1, seq2_match, quals2):
        if r1 in bases and r1 == r2:
            ML_bases.append(r1)
        elif set([r1, r2]) <= bases and q1 > 2 and q2 > 2:
            r1_score = log_p_struct[r1][r1][q1] + log_p_struct[r1][r2][q2]
            r2_score = log_p_struct[r2][r1][q1] + log_p_struct[r2][r2][q2]
            if r1_score > r2_score:
                ML_bases.append(r1)
            else:
                ML_bases.append(r2)
        elif r1 in bases and q1 > 2:
            ML_bases.append(r1)
        elif r2 in bases and q2 > 2:
            ML_bases.append(r2)
        else:
            ML_bases.append(r1)
    return ''.join(ML_bases)
    
def parse_fastq_lines(gzipped_filename):
    with gzip.open(gzipped_filename) as fh:
        for record in SeqIO.parse(fh, 'fastq'): 
            yield record
            
def rand_seq(seq_len):
    return ''.join(random.choice('ACGT') for _ in xrange(seq_len))
    
def write_read_names_by_sequence(read_names_given_seq, out_file_path):
    with open(out_file_path, 'w') as out:
        for seq, read_names in sorted(read_names_given_seq.items()):
            out.write('{}\t{}\n'.format(seq, '\t'.join(read_names)))