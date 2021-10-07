#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 10:10:50 2020

@author: dulab19
"""

import re
import pickle
import argparse
import os

def count_mutation(path):
    """
    read AA mutation number of V region
    Parameters
    ----------
    path : str
        file path of V-REGION-AA-change-statistics.txt

    Returns
    -------
    cdr_shm : dic
        shm number of each cdr region.
        {ac:[cdr shm nunber]}
    V_shm : dic
        shm number of V-region.
    fr_shm : dic
        shm number of each fr region.
    """
    regex = re.compile('[Ss]imilar')
    regex1 = re.compile('[0-9]{1,2}')
    header = []
    cdr_shm = {}
    cdr_ind = [[], [], []]
    fr_ind = []
    ind_len = 0
    V_shm = {}
    fr_shm = {}
    v_ind = 0
    error_num = 0
    mark =  False
    if 'hcv' in path:
        mark =  True
    with open(path, 'r') as f:
        for line in f:
            line = line.strip('\n').split('\t')
            if header == []:
                header = line
                ind_len = len(header)
                for i in range(ind_len):
                    col_name = header[i]
                    if 'CDR1' in col_name and regex.search(col_name):
                        cdr_ind[0].append(i)
                    elif 'CDR2' in col_name and regex.search(col_name):
                        cdr_ind[1].append(i)
                    elif 'CDR3' in col_name and regex.search(col_name):
                        cdr_ind[2].append(i)
                    elif 'V-REGION' in col_name and 'AA changes' in col_name:
                        v_ind = i
                    elif 'AA changes' in col_name and "CDR" not in col_name:
                        fr_ind.append(i)
                continue
            sq_number = line[1].split(' ')[0]
            if mark:
                sq_number = sq_number.replace(',', '_')
            if len(line) < ind_len:
                error_num += 1
                continue
            cdr_shm[sq_number] = []
            fr_shm[sq_number] = []
            num_list = regex1.findall(line[v_ind])  # shm number for V-region
            V_shm[sq_number] = max([int(num) for num in num_list])
            # record cdr1 cdr2 cdr3 similar and dissimilar AA changes
            for inds in cdr_ind:
                # inds : similar index and dissimilar index
                tmp = [line[i] for i in range(ind_len) if i in inds]
                out_vet = []
                # check () and +
                for ele in tmp:
                    num_list = regex1.findall(ele)
                    if num_list != []:
                        num_list = [int(num) for num in num_list]
                        out_vet.append(max(num_list))
                    else:
                        out_vet.append(0)
                cdr_shm[sq_number].append(out_vet) 
            for ind in fr_ind:
                num_list = regex1.findall(line[ind])
                if num_list == []:
                    num_list = [0]
                fr_shm[sq_number].append(max([int(num) for num in num_list])) 
    return cdr_shm, V_shm, fr_shm

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, help="IMGT HighV-QUEST results path")
parser.add_argument('-o', type=str, default="./original_model/data", help="output path")
args = parser.parse_args()
cdr_dic = {}  # record cdr shm
# get all subfolder names
for top, dirs, _ in os.walk(args.i):
    names = dirs
    break  # do not search in subfoldes
v_dic = {}  # record v-region shm
fr_dic = {}  # record fr region shm
for name in names:
    path = os.path.join(args.i, "{}/9_V-REGION-AA-change-statistics.txt".format(name))
    sub_cdr, sub_v, sub_fr = count_mutation(path)
    cdr_dic.update(sub_cdr)
    v_dic.update(sub_v)
    fr_dic.update(sub_fr)
path = os.path.join(args.o, "cdr_shm_dic.txt")
with open(path, 'wb') as f:
    pickle.dump(cdr_dic, f)
path = os.path.join(args.o, "v_shm_dic.txt")
with open(path, 'wb') as f:
    pickle.dump(v_dic, f)
path = os.path.join(args.o, "fr_shm_dic.txt")
with open(path, 'wb') as f:
    pickle.dump(fr_dic, f)