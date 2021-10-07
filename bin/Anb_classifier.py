#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 15:15:29 2020

@author: dulab19
"""
import sys 
sys.path.append("..") 
import utils
import os
from processfeature import ProcessFeature
import processdata as pda
import numpy as np
import pickle
import argparse
from collections import Counter
import stacking
from sklearn import metrics


parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, help="IMGT HighV-QUEST results file dir." +
                    "\n--input dir\n\tresult file1\n\tresult file2...\n")
parser.add_argument('-o', type=str, help="Output path")
parser.add_argument('--remove_redandunt', type=bool, default=True,
                    help="True or False.\n"+
                    "It is faster when running on non-redandunt data.")
parser.add_argument('--format', type=str, help="Output file format of prediction", choices=['csv', 'tsv'], default='csv')
parser.add_argument('--cpu', type=str, help="Specify how many cores to use", default=1)
parser.add_argument('--olga_path', type=str, help="The absolute path of python script: compute_pgen.py. /xxx/xx/olga/compute_pgen.py")
parser.add_argument('--pssm_path', type=str, help="The absolute path of PSSM result files")
parser.add_argument('--mode', type=str, help="validate or predict", default="predict")
parser.add_argument('--antigen_names', type=str, help="The folder names of IMGT HighV-QUEST results files, split by ','", default="hiv,flu,pps,acpa,tt,hcv,hbv")
args = parser.parse_args()

# process input files
if os.path.exists(args.i) is False:
    raise FileNotFoundError("input file path does not exist")
    exit(0)
print("processing IMGT HighV-QUEST results")
# process IMGT result files
names =  args.antigen_names.split(',')
name_index = {i: name for i, name in enumerate(names)}
for top, dirs, _ in os.walk(args.i):
    for dir_ in names:
        tmp_path = os.path.join(args.i, dir_)
        print(tmp_path)
        os.system(r"python3 ../joint_imgt.py -i {} -o ./data".format(tmp_path))
        
    break
os.system(r"python3 ../imgt_shm.py -i {} -o ./data".format(args.i))
ac_regions, CDR3, unqualifyed_ac, vj = pda.process_data(name_index)
n_acs, n_sqs, n_labels = pda.remove_unqualified(name_index, CDR3, unqualifyed_ac, True)
# remove and shuffle redundant sqs
if args.remove_redandunt is True:
    n_acs, n_sqs, n_labels = pda.remove_redundant(n_acs, n_sqs, n_labels)
print(Counter(n_labels))
print("processing done")
# cat fasta files of input
cmd = ""
for name in names:
    cmd += "./data/fasta/" + name + "_p.fasta "
os.system(r"cat {} > ./data/fasta/all.fasta".format(cmd))
# use olga to generate probability
olga_path = args.olga_path
olga = []
for x,y in zip(n_acs, n_sqs):
    tmp = ['C'+y+'W'] + vj[x].split('|')
    olga.append(tmp)
np.savetxt("./data/olga/olga_cdr3s.tsv", olga, fmt="%s", delimiter='\t')
del olga
os.system(r"python {} -i {} --humanIGH --display_off off --v_in 1 --j_in 2 -o ./data/olga/out_pgen.tsv".format(olga_path, "./data/olga/olga_cdr3s.tsv")) 
np.savetxt("./data/olga/all_acs.csv", n_acs, fmt="%s", delimiter='\n')

# load stacking model
with open("./model/stacking_model.txt", 'rb') as f:
    stack_model = pickle.load(f)
# count feature vector for each antibody
pf = ProcessFeature(ac_regions, args.pssm_path)
pf.read_external_files(n_acs)
matrix = pf.count_features(n_sqs, n_acs)
matrix = np.array(matrix)
y_train = np.loadtxt("./data/y_train.csv", delimiter=',')
# the index of stack_model predictions
ag_names = ['Anti-HIV-1 Ab', 'Anti-FLU Ab', 'Anti-PPS Ab', 'ACPA', 'Anti-TT Ab', 'Anti-HBV Ab']
stack_p = stack_model.predict(matrix, y_train)
stack_p = [ag_names[int(label)] for label in stack_p]
stack_p = np.array(stack_p)
y_score = stack_model.predict_proba(y_train)
print(y_score.shape)
output = np.concatenate((stack_p.reshape((-1, 1)), y_score), axis = 1)
# save output
header = "predicted_label," + ",".join(ag_names)
if args.format == 'csv':
    np.savetxt(os.path.join(args.o, 'prediction.csv'), output, fmt='%s', delimiter=',', header=header)
else:
    np.savetxt(os.path.join(args.o, 'prediction.tsv'), output, fmt='%s', delimiter='\t', header=header)
if args.mode == 'validate':
    print(metrics.classification_report(n_labels, stack_p, digits=4))
else:
    print(Counter(stack_p))
