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
    

def run_hhblits(acc, db_prefix, fasta_file, we, cpus=1):
    #"/home/cas/software/hh-suite/build/bin/hhblits -i temp/"+id+".fasta
    #-d /mnt/fat1/databases/HHBlits/UniRef30_2020_02/UniRef30_2020_02
    #-n 2 -oa3m temp/"+id+".a3m -ohhm temp/"+id+".hhm -o /dev/null -cpu 2"
    hhblits_a3m_out = we.createFile(acc+".hhblits.", ".a3m")
    hhblits_aln_out = we.createFile(acc+".hhblits.", ".aln")
    hhblits_hhm_out = we.createFile(acc+".hhblits.", ".hhm")
    hhblits_stdout = we.createFile(acc+".hhblits.stdout.", ".log")
    hhblits_stderr = we.createFile(acc+".hhblits.stderr.", ".log")

    #sequence = "".join([x.strip() for x in open(fastaFile).readlines()[1:]])

    try:
        subprocess.check_output(['hhblits', '-i', fasta_file,
                                 '-d', db_prefix,
                                 '-n', str(cpus),
                                 '-oa3m', hhblits_a3m_out,
                                 '-ohhm', hhblits_hhm_out,
                                 '-o', hhblits_stdout],
                                 stderr=open(hhblits_stderr, 'w'))
        a3m_to_aln(hhblits_a3m_out, hhblits_aln_out)
    except:
        logging.error("HHblits failed. For details, please see stderr file %s" % hhblits_stderr)
        raise
    return hhblits_aln_out, hhblits_hhm_out
