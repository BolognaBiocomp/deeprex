import numpy
from keras.models import load_model
from keras import backend as K

def build_sequence_profile(acc, aln_file, we):
    aa_order = '-ARNDCQEGHILKMFPSTWYV'
    sequence_profile_file = we.createFile(acc+".profile.", ".prof")
    matrix = []
    #Parse aln file
    with open(aln_file) as iif:
        seq_all = [line.rstrip() for line in iif.readlines() if len(line)>1]
        l = len(seq_all[0])
        iif.close()

    #For each column:
    for i in range(l):
        counts = numpy.zeros(len(aa_order))
        for seq in seq_all:
            try:
                counts[aa_order.index(seq[i])] += 1.0
            except:
                pass
        n = numpy.sum(counts)
        counts[1:] /= numpy.sum(counts[1:])
        counts[0] /= n
        matrix.append(counts)
    matrix = numpy.array(matrix)
    numpy.savetxt(sequence_profile_file, matrix, fmt="%.2f")
    return sequence_profile_file

def encode_protein(sequence, sequence_profile_file, hhblits_hhm_out):
    aa_order = '-ARNDCQEGHILKMFPSTWYV'
    profile = numpy.loadtxt(sequence_profile_file)
    hhm_file  = open(hhblits_hhm_out)
    hhm_line = hhm_file.readline()
    while hhm_line[:4] != "HMM ":
        hhm_line = hhm_file.readline()
    hhm_file.readline()
    hhm_file.readline()
    prot = []
    for (i, residue) in enumerate(sequence):
        one_hot = [0.0] * 20
        try:
            one_hot[aa_order.index(residue)-1] = 1.0
        except:
            pass
        one_hot = numpy.array(one_hot)

        hhm_line = hhm_file.readline().strip()
        hhm_line = hhm_line.split()

        #Parse hhm profile
        hhm_profile = []
        for prob in hhm_line[2:-1]:
            try:
                hhm_profile.append(2**(float(prob)/-1000))
            except ValueError:
                hhm_profile.append(0.0)
        hhm_profile = numpy.array(hhm_profile)
        if hhm_profile.sum() == 0:
            hhm_profile = one_hot
        else:
            hhm_profile /= hhm_profile.sum()

        #Parse transitions and neff from hhm file
        tra_neff = hhm_file.readline().strip().split('\t')
        transitions = []
        for tra in tra_neff[:7]:
            try:
                transitions.append(2**(float(tra)/-1000))
            except ValueError:
                transitions.append(0.0)
        transitions = numpy.array(transitions)
        neff = numpy.array(tra_neff[7:],float)/1000
        hhm_file.readline()
        vect = numpy.concatenate((one_hot, profile[i], hhm_profile, transitions, neff))
        prot.append(vect)
    prot = numpy.array([prot])
    return prot

def predict(protein, model_file):
    model = load_model(model_file)
    xs = K.constant(protein)
    predictions = model.predict_on_batch(xs).tolist()[0]
    return predictions

def write_tsv_output(acc, sequence, predictions, out_file):
    for i in range(len(predictions)):
        label = 0
        if predictions[i][0] > 0.5:
            label = 1
        reliability = round(2.0 * numpy.absolute(predictions[i][0] - 0.5), 2)
        print(acc, i+1, sequence[i], label, reliability, sep="\t", file=out_file)