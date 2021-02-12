################################################################################
# score_conservation.py - Copyright Tony Capra 2007 - Last Update: 03/09/11
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# -----------------------------------------------------------------------------

import math

PSEUDOCOUNT = .0000001

amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
iupac_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V", "W", "Y", "Z", "X", "*", "-"]

# dictionary to map from amino acid to its row/column in a similarity matrix
aa_to_index = {}
for i, aa in enumerate(amino_acids):
    aa_to_index[aa] = i

def weighted_gap_penalty(col, seq_weights):
    """ Calculate the simple gap penalty multiplier for the column. If the
    sequences are weighted, the gaps, when penalized, are weighted
    accordingly. """

    # if the weights do not match, use equal weight
    if len(seq_weights) != len(col):
        seq_weights = [1.] * len(col)
    gap_sum = 0.
    for i in range(len(col)):
        if col[i] == '-':
            gap_sum += seq_weights[i]
    return 1 - (gap_sum / sum(seq_weights))

def gap_percentage(col):
    """Return the percentage of gaps in col."""
    num_gaps = 0.
    for aa in col:
        if aa == '-': num_gaps += 1
    return num_gaps / len(col)

def weighted_freq_count_pseudocount(col, seq_weights, pc_amount):
    """ Return the weighted frequency count for a column--with pseudocount."""

    # if the weights do not match, use equal weight
    if len(seq_weights) != len(col):
        seq_weights = [1.] * len(col)

    aa_num = 0
    freq_counts = len(amino_acids)*[pc_amount] # in order defined by amino_acids

    for aa in amino_acids:
        for j in range(len(col)):
            if col[j] == aa:
                freq_counts[aa_num] += 1 * seq_weights[j]

        aa_num += 1

    for j in range(len(freq_counts)):
        freq_counts[j] = freq_counts[j] / (sum(seq_weights) + len(amino_acids) * pc_amount)

    return freq_counts

def shannon_entropy(col, seq_weights, gap_penalty=1):
    """Calculates the Shannon entropy of the column col. sim_matrix  and
    bg_distr are ignored. If gap_penalty == 1, then gaps are penalized. The
    entropy will be between zero and one because of its base. See p.13 of
    Valdar 02 for details. The information score 1 - h is returned for the sake
    of consistency with other scores."""

    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

    h = 0.
    for i in range(len(fc)):
        if fc[i] != 0:
            h += fc[i] * math.log(fc[i])

    h /= math.log(min(len(fc), len(col)))

    inf_score = 1 - (-1 * h)

    if gap_penalty == 1:
        return inf_score * weighted_gap_penalty(col, seq_weights)
    else:
        return inf_score

def calculate_sequence_weights(msa):
    """ Calculate the sequence weights using the Henikoff '94 method
    for the given msa. """

    seq_weights = [0.] * len(msa)
    for i in range(len(msa[0])):
        freq_counts = [0] * len(amino_acids)
        col = []
        for j in range(len(msa)):
            if msa[j][i] != '-' and msa[j][i] in amino_acids: # ignore gaps
                freq_counts[aa_to_index[msa[j][i]]] += 1
        num_observed_types = 0
        for j in range(len(freq_counts)):
            if freq_counts[j] > 0: num_observed_types +=1
        for j in range(len(msa)):
            if msa[j][i] in amino_acids:
                d = freq_counts[aa_to_index[msa[j][i]]] * num_observed_types
                if d > 0:
                    seq_weights[j] += 1. / d
    for w in range(len(seq_weights)):
            seq_weights[w] /= len(msa[0])
    return seq_weights

def read_fasta_alignment(filename):
    """ Read in the alignment stored in the FASTA file, filename. Return two
    lists: the identifiers and sequences. """

    f = open(filename)

    names = []
    alignment = []
    cur_seq = ''

    for line in f:
        line = line[:-1]
        if len(line) == 0: continue

        if line[0] == ';': continue
        if line[0] == '>':
            names.append(line[1:].replace('\r', ''))

            if cur_seq != '':
                cur_seq = cur_seq.upper()
                for i, aa in enumerate(cur_seq):
                    if aa not in iupac_alphabet:
                        cur_seq = cur_seq.replace(aa, '-')
                alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
                cur_seq = ''
        elif line[0] in iupac_alphabet:
            cur_seq += line.replace('\r', '')

    # add the last sequence
    cur_seq = cur_seq.upper()
    for i, aa in enumerate(cur_seq):
        if aa not in iupac_alphabet:
            cur_seq = cur_seq.replace(aa, '-')
    alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))

    return names, alignment

def get_column(col_num, alignment):
    """Return the col_num column of alignment as a list."""
    col = []
    for seq in alignment:
        if col_num < len(seq): col.append(seq[col_num])
    return col

def window_score(scores, window_len, lam=.5):
    """ This function takes a list of scores and a length and transforms them
    so that each position is a weighted average of the surrounding positions.
    Positions with scores less than zero are not changed and are ignored in the
    calculation. Here window_len is interpreted to mean window_len residues on
    either side of the current residue. """

    w_scores = scores[:]
    for i in range(window_len, len(scores) - window_len):
        if scores[i] < 0:
            continue
        sum = 0.
        num_terms = 0.
        for j in range(i - window_len, i + window_len + 1):
            if i != j and scores[j] >= 0:
                num_terms += 1
                sum += scores[j]
        if num_terms > 0:
            w_scores[i] = (1 - lam) * (sum / num_terms) + lam * scores[i]
    return w_scores

def score_conservation(align_file, gap_cutoff=0.7, window_size=1,
                       use_gap_penalty=1, win_lam=.5):
    names, alignment = read_fasta_alignment(align_file)
    seq_len = len(alignment[0])
    seq_weights = calculate_sequence_weights(alignment)
    scores = []
    for i in range(seq_len):
        col = get_column(i, alignment)
        if gap_percentage(col) <= gap_cutoff:
            scores.append(shannon_entropy(col, seq_weights, use_gap_penalty))
        else:
            scores.append(0.0)
    if window_size > 0:
        scores = window_score(scores, window_size, win_lam)
    return scores
