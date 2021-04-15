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

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

from Bio import SeqIO

from deeprexlib import workenv
import deeprexlib.deeprexconfig as cfg
from deeprexlib import hhblits
from deeprexlib import utils
from deeprexlib import conservation

def run_fasta(ns):
    try:
        we = workenv.TemporaryEnv()
        data_cache = utils.get_data_cache()
        ofs = open(ns.outf, 'w')
        for record in SeqIO.parse(ns.fasta, 'fasta'):
            prefix = record.id.replace("|","_")
            fasta_seq  = we.createFile(prefix+".", ".fasta")
            SeqIO.write([record], fasta_seq, 'fasta')
            utils.print_date("Running HHBlits and building sequence profile [protein=%s]" % record.id)
            hhblits_aln_out, hhblits_hhm_out, msa_conservation = hhblits.run_hhblits(prefix, ns.hhblits_db, fasta_seq, we, gap_cutoff=ns.gapth, cpus=ns.cpus, data_cache=data_cache)
            sequence_profile_file = utils.build_sequence_profile(prefix, hhblits_aln_out, we)
            sequence = str(record.seq)
            hydrophobicity = utils.score_hp(sequence, ns.hpwin)
            utils.print_date("Encode protein sequence [protein=%s]" % record.id)
            protein = utils.encode_protein(sequence, sequence_profile_file, hhblits_hhm_out)
            utils.print_date("Predict residue solvent exposure [protein=%s]" % record.id)
            predictions = utils.predict(protein, cfg.DEEPREX_MODEL_FILE)
            utils.print_date("Writing predictions to TSV file [protein=%s]" % record.id)
            utils.write_tsv_output(record.id, sequence, predictions,
                                   hydrophobicity, msa_conservation, ofs)
    except:
        ofs.close()
        logging.exception("Errors occurred:")
        sys.exit(1)
    else:
        utils.print_date("Cleaning temporary environment end exit")
        ofs.close()
        we.destroy()
        sys.exit(0)

def run_aln(ns):
    try:
        we = workenv.TemporaryEnv()
        record = SeqIO.read(ns.fasta, 'fasta')
        prefix = record.id.replace("|","_")
        fasta_seq  = we.createFile(prefix+".", ".fasta")
        SeqIO.write([record], fasta_seq, 'fasta')
        aln_file = we.createFile(prefix+".", ".aln")
        faln_file = we.createFile(prefix+".", ".aln.fa")
        hhblits.a3m_to_aln(ns.a3m, aln_file)
        hhblits.aln_to_faln(aln_file, faln_file)
        msa_conservation = conservation.score_conservation(faln_file, gap_cutoff=ns.gapth)
        sequence_profile_file = utils.build_sequence_profile(prefix, aln_file, we)
        sequence = str(record.seq)
        hydrophobicity = utils.score_hp(sequence, ns.hpwin)
        utils.print_date("Encode protein sequence [protein=%s]" % record.id)
        protein = utils.encode_protein(sequence, sequence_profile_file, ns.hhm)
        utils.print_date("Predict residue solvent exposure [protein=%s]" % record.id)
        predictions = utils.predict(protein, cfg.DEEPREX_MODEL_FILE)
        utils.print_date("Writing predictions to TSV file [protein=%s]" % record.id)
        ofs = open(ns.outf, 'w')
        utils.write_tsv_output(record.id, sequence, predictions,
                               hydrophobicity, msa_conservation, ofs)
    except:
        logging.exception("Errors occurred:")
        sys.exit(1)
    else:
        utils.print_date("Cleaning temporary environment end exit")
        ofs.close()
        we.destroy()
        sys.exit(0)

def run_ss(ns):
    try:
        we = workenv.TemporaryEnv()
        ofs = open(ns.outf, 'w')
        for record in SeqIO.parse(ns.fasta, 'fasta'):
            sequence = str(record.seq)
            hydrophobicity = utils.score_hp(sequence, ns.hpwin)
            utils.print_date("Encode protein sequence [protein=%s]" % record.id)
            protein = utils.encode_protein_single_seq(sequence)
            utils.print_date("Predict residue solvent exposure [protein=%s]" % record.id)
            predictions = utils.predict(protein, cfg.DEEPREX_MODEL_FILE)
            utils.print_date("Writing predictions to TSV file [protein=%s]" % record.id)
            
            utils.write_tsv_output(record.id, sequence, predictions,
                                   hydrophobicity, [0.0]*len(sequence), ofs)
    except:
        logging.exception("Errors occurred:")
        sys.exit(1)
    else:
        utils.print_date("Cleaning temporary environment end exit")
        ofs.close()
        we.destroy()
        sys.exit(0)

def main():
    DESC="DeepREx: Deep learning-based predictor of Residue EXposure"
    parser = argparse.ArgumentParser(description=DESC)
    subparsers   = parser.add_subparsers(title = "subcommands",
                                         description = "valid subcommands",
                                         help = "additional help",
                                         required = True)
    multifasta  = subparsers.add_parser("fasta",
                                          help = "Multi-FASTA input module",
                                          description = "DeepREx: Multi-FASTA input module.")
    aln  = subparsers.add_parser("aln", help = "Alignment input module (one sequence at a time)",
                                    description = "DeepREx: Alignment input module.")
    singless = subparsers.add_parser("singleseq", help = "Single sequence input module (one sequence at a time)",
                                    description = "DeepREx: Single sequence input module.")
    singless.add_argument("-f", "--fasta",
                        help = "The input multi-FASTA file name",
                        dest = "fasta", required = True)

    singless.add_argument("-o", "--outf",
                        help = "The output TSV file",
                        dest = "outf", required = True)

    singless.add_argument("-w", "--hp-window",
                        help = "Window size for hydrophobicity computation (default: 5)",
                        dest = "hpwin", required = False, type = int, default= 5)
    singless.add_argument("-t", "--cpus",
                        help = "Number of CPUs to use",
                        dest = "cpus", type = int, default = 1)
    singless.set_defaults(func=run_ss)
    multifasta.add_argument("-f", "--fasta",
                        help = "The input multi-FASTA file name",
                        dest = "fasta", required = True)
    multifasta.add_argument("-d", "--hhblits-database",
                        help = "The HHBlits database file",
                        dest = "hhblits_db", required = True)
    multifasta.add_argument("-o", "--outf",
                        help = "The output TSV file",
                        dest = "outf", required = True)
    multifasta.add_argument("-g", "--gap-th",
                        help = "Score conservation of MSA columns with less than this gap threshold (default: 0.7)",
                        dest = "gapth", required = False, type = float, default= 0.7)
    multifasta.add_argument("-w", "--hp-window",
                        help = "Window size for hydrophobicity computation (default: 5)",
                        dest = "hpwin", required = False, type = int, default= 5)
    multifasta.add_argument("-t", "--cpus",
                        help = "Number of CPUs to use",
                        dest = "cpus", type = int, default = 1)
    multifasta.set_defaults(func=run_fasta)
    aln.add_argument("-f", "--fasta",
                        help = "The input FASTA file name (single sequence)",
                        dest = "fasta", required = True)
    aln.add_argument("-a", "--a3m",
                        help = "The input multiple sequence alignment in A3M format",
                        dest = "a3m", required = True)
    aln.add_argument("-m", "--hhm",
                        help = "The input HHM file from HHblits",
                        dest = "hhm", required = True)
    aln.add_argument("-o", "--outf",
                        help = "The output TSV file",
                        dest = "outf", required = True)
    aln.add_argument("-g", "--gap-th",
                        help = "Score conservation of MSA columns with less than this gap threshold (default: 0.7)",
                        dest = "gapth", required = False, type = float, default= 0.7)
    aln.add_argument("-w", "--hp-window",
                        help = "Window size for hydrophobicity computation (default: 5)",
                        dest = "hpwin", required = False, type = int, default= 5)
    aln.add_argument("-t", "--cpus",
                        help = "Number of CPUs to use",
                        dest = "cpus", type = int, default = 1)
    aln.set_defaults(func=run_aln)
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        ns = parser.parse_args()
        ns.func(ns)

if __name__ == "__main__":
    main()
