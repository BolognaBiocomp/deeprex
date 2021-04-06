#!/usr/bin/env python
import sys
import numpy as np

def read_data(filename):
    ypred = []
    ytrue = []
    cls_map = {"Buried": 0, "Exposed": 1}
    for line in open(filename).readlines()[1:]:
        line = line.split()
        ypred.append(cls_map[line[3]])
        ytrue.append(cls_map[line[7]])
    return ytrue, ypred

def confusion_matrix(ytrue, ypred):
    cm = np.zeros((2,2))
    for i in range(len(ytrue)):
        cm[ytrue[i],ypred[i]]+=1
    return cm

def mcc(cm):
    n = cm[0,0]*cm[1,1]-cm[0,1]*cm[1,0]
    d = (cm[0,0]+cm[0,1])*(cm[0,0]+cm[1,0])*(cm[1,1]+cm[0,1])*(cm[1,1]+cm[1,0])
    if d > 0:
        mcc = n / np.sqrt(d)
    else:
        mcc = 0.0
    return mcc

def prec(cm):
    d = cm[1,1]+cm[0,1]
    if d > 0:
        prec = cm[1,1]/d
    else:
        prec = 0.0
    return prec

def sens(cm):
    return cm[1,1]/(cm[1,1]+cm[1,0])

def q2(cm):
    return (cm[0,0]+cm[1,1])/np.sum(cm)

def f1(cm):
    sens = sens(cm)
    prec = prec(cm)
    if sens+prec > 0:
        return 2*sens*prec/(sens+prec)
    else:
        return 0.0

cm = read_data(sys.argv[1])
print("SENSITIVITY:", sens(cm), sep="\t")
print("PRECISION:", prec(cm), sep="\t")
print("F1:", f1(cm), sep="\t")
print("Q2:", q2(cm), sep="\t")
print("MCC:", mcc(cm), sep="\t")
