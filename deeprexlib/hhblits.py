import subprocess
import logging
import re
from . import deeprexconfig as cfg
from . import conservation as cons

def a3m_to_aln(a3m_file, aln_file):
    of = open(aln_file, 'w')
    iif = open(a3m_file)
    for line in iif.readlines():
        if not re.match('^>', line):
            of.write(re.sub('[a-z]', '', line))
    iif.close()
    of.close()

def aln_to_faln(aln_file, faln_file):
    aln_f = open(aln_file)
    o_fal = open(faln_file, 'w')
    seq_c = 0
    for line in aln_f.readlines():
        print(">seq_%d" % seq_c, file=o_fal)
        print(line.rstrip(), file=o_fal)
        seq_c = seq_c + 1
    o_fal.close()
    aln_f.close()
    return seq_c

def run_hhblits(acc, db_prefix, fasta_file, we, gap_cutoff=0.7, cpus=1, data_cache=None):
    hhblits_a3m_out = we.createFile(acc+".hhblits.", ".a3m")
    hhblits_aln_out = we.createFile(acc+".hhblits.", ".aln")
    hhblits_hhm_out = we.createFile(acc+".hhblits.", ".hhm")
    hhblits_fal_out = we.createFile(acc+".hhblits.", ".aln.fa")
    hhblits_stdout = we.createFile(acc+".hhblits.stdout.", ".log")
    hhblits_stderr = we.createFile(acc+".hhblits.stderr.", ".log")

    try:
        exec_hhblits = True
        sequence = "".join([x.strip() for x in open(fasta_file).readlines()[1:]])
        if data_cache is not None:
            if data_cache.lookup(sequence, 'hhblits.aln'):
                if data_cache.lookup(sequence, 'hhblits.hhm'):
                    if data_cache.lookup(sequence, 'hhblits.a3m'):
                        exec_hhblits = False
        if exec_hhblits:
            subprocess.check_output(['hhblits', '-i', fasta_file,
                                    '-d', db_prefix,
                                    '-n', "2",
                                    '-cpu', str(cpus),
                                    '-oa3m', hhblits_a3m_out,
                                    '-ohhm', hhblits_hhm_out,
                                    '-o', hhblits_stdout],
                                    stderr=open(hhblits_stderr, 'w'))
            a3m_to_aln(hhblits_a3m_out, hhblits_aln_out)
            if data_cache is not None:
                data_cache.store(hhblits_aln_out, sequence, 'hhblits.aln')
                data_cache.store(hhblits_hhm_out, sequence, 'hhblits.hhm')
                data_cache.store(hhblits_a3m_out, sequence, 'hhblits.a3m')
        else:
            data_cache.retrieve(sequence, 'hhblits.aln', hhblits_aln_out)
            data_cache.retrieve(sequence, 'hhblits.hhm', hhblits_hhm_out)
        seq_c = aln_to_faln(hhblits_aln_out, hhblits_fal_out)
        if seq_c > 1:
            msa_conservation = cons.score_conservation(hhblits_fal_out, gap_cutoff=gap_cutoff)
        else:
            msa_conservation = [0] * len(sequence)
    except:
        logging.error("HHblits failed. For details, please see stderr file %s" % hhblits_stderr)
        raise
    return hhblits_aln_out, hhblits_hhm_out, msa_conservation
