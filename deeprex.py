#!/usr/bin/env python
import sys
import os
import argparse
import logging
if 'DEEPREX_ROOT' in os.environ:
    sys.path.append(os.environ['DEEPREX_ROOT'])
else:
    logging.error("DEEPREX_ROOT environment varible is not set")
    logging.error("Please, set and export DEEPREX_ROOT to point to deeprex root folder:")
    logging.error("$ export DEEPREX_ROOT=/path/to/deeprex")
    sys.exit(1)

from Bio import SeqIO

from deeprexlib import workenv
import deeprexlib.deeprexconfig as cfg
from deeprexlib import hhblits
from deeprexlib import utils

def main():
    DESC="DeepREx: Deep learning-based predictor of Residue EXposure"
    parser = argparse.ArgumentParser(description=DESC)
    parser.add_argument("-f", "--fasta",
                        help = "The input multi-FASTA file name",
                        dest = "fasta", required = True)
    parser.add_argument("-d", "--hhblits-database",
                        help = "The HHBlits database file",
                        dest = "hhblits_db", required = True)
    parser.add_argument("-o", "--outf",
                        help = "The output TSV file",
                        dest = "outf", required = True)
    parser.add_argument("-a", "--cpus",
                        help = "Number of CPUs to use",
                        dest = "cpus", type = int, default = 1)
    ns = parser.parse_args()
    try:
        we = workenv.TemporaryEnv()
        ofs = open(ns.outf, 'w')
        for record in SeqIO.parse(ns.fasta, 'fasta'):
            prefix = record.id.replace("|","_")
            fasta_seq  = we.createFile(prefix+".", ".fasta")
            SeqIO.write([record], fasta_seq, 'fasta')
            hhblits_aln_out, hhblits_hhm_out = hhblits.run_hhblits(prefix, ns.hhblits_db, fasta_seq, we, cpus=ns.cpus)
            sequence_profile_file = utils.build_sequence_profile(prefix, hhblits_aln_out, we)
            sequence = str(record.seq)
            protein = utils.encode_protein(sequence, sequence_profile_file, hhblits_hhm_out)
            predictions = utils.predict(protein, cfg.DEEPREX_MODEL_FILE)
            utils.write_tsv_output(record.id, sequence, predictions, ofs)
    except:
        ofs.close()
        logging.exception("Errors occurred:")
        sys.exit(1)
    else:
        ofs.close()
        we.destroy()
        sys.exit(0)

if __name__ == "__main__":
    main()
