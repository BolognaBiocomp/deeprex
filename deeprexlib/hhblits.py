import subprocess
import logging
import re
from . import deeprexconfig as cfg

def a3m_to_aln(a3m_file, aln_file):
    of = open(aln_file, 'w')
    iif = open(a3m_file)
    for line in iif.readlines():
        if not re.match('^>', line):
            of.write(re.sub('[a-z]', '', line))
    iif.close()
    of.close()

def run_hhblits(acc, db_prefix, fasta_file, we, cpus=1, data_cache=None):
    hhblits_a3m_out = we.createFile(acc+".hhblits.", ".a3m")
    hhblits_aln_out = we.createFile(acc+".hhblits.", ".aln")
    hhblits_hhm_out = we.createFile(acc+".hhblits.", ".hhm")
    hhblits_stdout = we.createFile(acc+".hhblits.stdout.", ".log")
    hhblits_stderr = we.createFile(acc+".hhblits.stderr.", ".log")

    try:
        exec_hhblits = True
        if data_cache is not None:
            sequence = "".join([x.strip() for x in open(fasta_file).readlines()[1:]])
            if data_cache.lookup(sequence, 'hhblits.aln'):
                if data_cache.lookup(sequence, 'hhblits.hhm'):
                    if data_cache.lookup(sequence, 'hhblits.a3m'):
                        exec_hhblits = False
        if exec_hhblits:
            subprocess.check_output(['hhblits', '-i', fasta_file,
                                    '-d', db_prefix,
                                    '-n', str(cpus),
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
    except:
        logging.error("HHblits failed. For details, please see stderr file %s" % hhblits_stderr)
        raise
    return hhblits_aln_out, hhblits_hhm_out
