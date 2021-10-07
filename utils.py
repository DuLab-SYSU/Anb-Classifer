#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 16:17:15 2020

@author: dulab19
"""
import os
import random
import numpy as np
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pickle


# functions for write and read files
# def add(i, out_list):
#     """
#     add ',' in list`s element
#     """
#     if i+10 >= len(out_list):
#         j = len(out_list)
#     else:
#         j = i+10
#     out_list = [str(ele) for ele in out_list]
#     tmp = out_list[i:j]
#     for k in range(0,len(tmp)):
#         tmp[k] += ','
#     tmp[-1] = tmp[-1][:-1]
#     return tmp

def write_out_file(path, out_list, explain="", m='w'):
    with open(path, m) as f:
        temp = [",".join(map(str, out_list[i:i+10])) for i in range(0, len(out_list), 10)]
        if explain != "":
            f.write(explain+'\n')
        for line in temp:
           f.write(line+'\n')
        f.write('\n')

def read_as_list(path, list_num):
    """
    read created csv file  .
    """
    list_out = []
    with open(path,'r') as f:
        marker = -1
        for line in f:
            line = line.strip('\n')
            if ":" in line :
                marker += 1
                list_out.append([])
            elif line == '':
                pass
            elif marker >= 0:
                list_out[marker] += line.strip(',').split(',')
    if list_num == 2:
        return list_out[0], list_out[1]
    elif list_num == 1:
        return list_out[0]
    elif list_num == 3:
        return list_out[0], list_out[1], list_out[2]
    else:
        return list_out

def write_fasta(path, ac_list, sq_list, b_list, m='w', mode='normal'):
    if mode == 'normal':
        with open(path, m) as f: 
            for ac,sq in zip(ac_list, sq_list):
                f.write(">{}\n".format(ac, b_list))
                f.write(sq+'\n')
    else:
        with open(path, m) as f: 
            for ac,sq in zip(ac_list,sq_list):
                f.write(">{} (({}))\n".format(ac, b_list))
                f.write(sq+'\n')
    return None

def read_dic_title(path = r"dic_title_1221.txt"):
    dic_title = {}
    with open(path,'r') as f:
        for line in f:
            line = line.strip('\n')
            if "REF" in line:
                ref = line.split("REF ")[-1]
                dic_title[ref] = []
            elif '//' in line:
                pass
            else:
                tmp = line.strip(',').split(',')
                dic_title[ref] += tmp
    return dic_title

def read_fasta(path):
    sq_dic = {}
    with open(path, 'r') as f:
        for line in f:
            if '>' in line:
                AC = line.strip('>').split('|')[0].strip('\n')
                if AC not in sq_dic:
                    sq_dic[AC] = ""
            else:
                sq_dic[AC] += line.strip('\n')
    return sq_dic

def read_config(path):
    out = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.rstrip().split(':')
            name = line[0].split(' ')[0]
            out[name] = line[-1]
    return out


def one_hot_feature(feature, threshold=0):
    vet = [0]*len(feature)#ag0 1 2 3…… |no only max
    max_num = max(feature)   
    index = [feature.index(ft) for ft in feature if ft == max_num]
    if len(index) == 1 and max_num >= threshold:#max score > threshold
        vet[index[0]] = 1 
    return vet


def one_hot_label(y_list, num):
    y_new = []
    for y in y_list:
        vet = [0]*num
        vet[y]= 1
        y_new.append(vet)
    return y_new


def sq_pi(sq):
    """
    count isoelectric point score. 
    """
    pep = ProteinAnalysis(sq)
    instability = pep.instability_index()
    return [instability]

def read_pssm():
    with open("./pssm_dict", 'rb') as f:
        pssm = pickle.load(f)
    return pssm

def extract_pssm(file_dir):
    row = re.compile("\s*\d{1,3}\s[A-Z]")
    number = re.compile('[-]*\d[.]*\d*')
    out = {}
    count = 0
    error = []
    for root, dir_, files in os.walk(file_dir):
        for f_name in files:
            cur_path = os.path.join(root, f_name)
            matrix = []
            count += 1
            with open(cur_path, 'r') as f:
                sq = ""
                for line in f:
                    match = re.match(row, line)
                    if match:
                        header = match.group(0)
                        tmp = line[len(header): ]
                        nums = re.findall(number, tmp)
                        pssm_20 = list(map(int, nums[:20]))
                        matrix.append(pssm_20)
                        sq += header[-1]
            matrix = np.array(matrix)
            out[sq] = matrix
    with open("./pssm_dict", 'wb') as f:
        pickle.dump(out, f)



def position_w_matrix(cdr3_list, length, column = 22):
    pwm = np.zeros((column, length))  # 22*n
    aa_index = {'A':0,'R':1,'N':2,"D":3,'C':4,'Q':5,'E':6,\
           'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,\
           'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'X':20,'-':21}
    for i in range(length):
        vet = [cdr3[i] for cdr3 in cdr3_list if i < len(cdr3)]
        for aa in vet:
            pwm[aa_index[aa]][i] += 1
        sum_n = np.sum(pwm[:,i])
        pwm[:,i] = (pwm[:,i]/sum_n).reshape(column)
    return pwm
        
def count_global_similarity(sq, pwms, sizes, length=15, st=0, num=7):
    aa_index = {'A':0,'R':1,'N':2,"D":3,'C':4,'Q':5,'E':6,\
           'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,\
           'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'X':20,'-':21}
    if sq == ['-']*40:
        return [0]*num
    sq = [aa_index[aa] for aa in sq]
    vet = [0]*num
    for i in range(num):
        weight = 0.005/np.log(sizes[i])
        for j in range(st,len(sq)):
            index = sq[j]
            if j == length:
                break
            elif pwms[i][index][j] == 0 or index == 20:
                vet[i] += -np.log2(weight)  # when an AA is never observed in pwm.
            else:
                vet[i] += -np.log2(pwms[i][index][j])
    sum_n = sum(vet)
    vet = [1-tmp/sum_n for tmp in vet]
    sum_n = sum(vet)
    vet = [tmp/sum_n for tmp in vet]
    return vet

def gly_position(sq):
    regex = re.compile("N[^P][S|T]")#N-glycosylation position
    match = re.findall(regex, sq)
    feature = 0
    if match != []:
        feature = len(match)
    return feature

def read_aaindex(path):
    pc_dic = {}
    tmp = []
    AA = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X']
    with open(path) as f:
        for line in f:
            line = line.strip().split()
            tmp += [float(line[i]) for i in range(len(line))]
    tmp.append(0.0)  # X: 0.0
    for x,y in zip(AA,tmp):
        pc_dic[x] = y
    return pc_dic


def AAC_PSSM(region_sqs, pssm, ac):
    vector = []
    st = len(region_sqs[0])
    for sq in region_sqs[1:-1]:
        length = len(sq)
        if length <= 1:
            print(ac, sq)
            if length == 1:
                aac = pssm[st]
            else:
                aac = [0]*20
            vector += list(aac)
            st += length
            continue
        aac = pssm[st: st+length].mean(axis=0) # [st: end]
        vector += list(aac)
        st += length
    end = st
    st = len(region_sqs[0])
    vector += list(pssm[st: end].mean(axis=0))
    return vector


def count_physiochemical(ac, sq, aaindexs, region):
    l1 = lambda x, y:0 if y==0 else x/y
    vet = []
    length_lim = {'cdr3':40, 'cdr2':10, 'cdr1':12, 'fr1':25, 'fr2':17,\
                  'fr3':38, 'fr4':11}
    sq_lim = length_lim[region]
    sq_length = len(sq)
    f_names =  ["helix probability", "hydrophobicity", "isoelectric point",\
               "sheet probability", "steric parameter", "volume", "polarisability"]
    # physiochemial encoding
    val_sum = [0]*7
    val_mean = [0]*7
    for aa in sq[:sq_lim]:
        for (f_n, i) in zip(f_names, range(7)):
            vet.append(aaindexs[f_n][aa])
            val_sum[i] += aaindexs[f_n][aa]
            val_mean[i] += aaindexs[f_n][aa]
    if sq_length < sq_lim:  # fill with 0
        vet += [-10]*7*(sq_lim - len(sq))
    vet += val_sum
    val_mean = [l1(val, sq_length) for val in val_mean]
    return vet+val_mean

def fasta_weblog(path, imgt_dic_full, acs):
    ac_list = []
    sq_list = []
    for ac in acs:
        if imgt_dic_full.get(ac):
            ac_list.append(ac)
            sq_list.append("".join(imgt_dic_full[ac][0]))
    write_fasta(path, ac_list, sq_list, "")
    return None


def read_olga(f_dir):
    pgen = []
    with open(f_dir + "out_pgen.tsv", 'r') as f:
        for line in f:
            tmp = line.rstrip('\n').split('\t')[1]
            if tmp == '0.0':
                tmp = '1e-70'
            pgen.append(float(tmp))
    acs = np.loadtxt(f_dir +"all_acs.csv",  delimiter='\n', dtype=str)
    pgen_dic = { x : y for x, y in zip(acs, pgen)}
    return pgen_dic


def read_imgt_index(mode):
    path = "./params/imgt_index_{}.csv".format(mode)
    header = []
    index_dic = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.strip('\n')
            if "header" in line:
                header = line.split(',')[1:]
            else:
                line = line.split(',')
                length = int(line[0])
                index_dic[length] = line[1:]
    return index_dic, header

def align_cdrs(ac_list, cdr_list, mode = 'cdr3', mark = 'skip'):
    """
    align cdrs in imgt format. 
    ----------
    ac_list : list
        serial numbers of sqs
    cdr_list : list
        cdr sqs
    mode : str, optional
        ['cdr3', 'cdr2', 'cdr1']. One of the three CDRs. The default is 'cdr3'.
    mark : str, optional
        'skip': skip '-' in sqs. output imgt index for each AA. The default is 'skip'.
        'full': output standard aligned imgt cdrs with '-' in sqs. 
    Returns
    -------
    imgt_dic : dictionary
        DESCRIPTION.
    header : list
        DESCRIPTION.

    """
    imgt_dic = {}
    index_dic, header = read_imgt_index(mode)
    num = len(header)
    sq_limit = {'cdr3':[5, 40], 'cdr1':[5,12], 'cdr2':[0, 10]}
    u_lim = sq_limit[mode][1]  # upper limitation
    l_lim = sq_limit[mode][0]  # lower limitation
    cdr_length = len(header)
    error_num = 0
    error_ac = []
    for sq, ac in zip(cdr_list, ac_list):
        sq_length = len(sq)
        if sq_length > u_lim or sq_length < l_lim or '*' in sq:
            error_num += 1
            error_ac.append(ac)
            continue
        inds = index_dic[sq_length]
        if mark == 'skip':  # skip. will not align sq with '-'
            sq_index = [header[i] for i in range(cdr_length) if inds[i] == '1']
            tmp = [aa for aa in sq]
            imgt_dic[ac] = [tmp, sq_index]
        else:  # converts a sq to a imgt-standard sq with '-'
            tmp = []
            ed = 0
            for i in range(sq_length):
                if inds[i] == '1':
                    tmp += [sq[i]]
                else:
                    ed = i
                    break
            if sq_length < u_lim:
                tmp += ['-']*(num - sq_length)
                tmp += [aa for aa in sq[ed:]]
            imgt_dic[ac] = [tmp, header]
    print(error_num, "{} sqs out of length range".format(mode))
    return imgt_dic, header, error_ac

def wrong_predicted_sequence(ag_index, sq_test, ac_test, y_p, y_t):
    """
    return FP sequences

    Parameters
    ----------
    ag_index : dictionary
        DESCRIPTION.
    sq_test : array-like
        DESCRIPTION.
    ac_test : array-like
        DESCRIPTION.
    y_p : array-like
        predicted labels.
    y_t : array-like
        true labels.

    Returns
    -------
    wrong_sq_and_ac : dic
        DESCRIPTION.

    """
    num = len(ag_index)
    wrong_sq_and_ac = {'{}{}'.format(i, j):[[],[],[]] for i in range(num) for j in range(num)}
    k = 0
    for i,j,x,y in zip(y_t, y_p, sq_test, ac_test):
        if i != j:
            key = str(i)+str(j)
            wrong_sq_and_ac[key][0].append(x)
            wrong_sq_and_ac[key][1].append(y)
            wrong_sq_and_ac[key][2].append(k)
        k += 1
    return wrong_sq_and_ac

def shuffle_data(acs, sqs, labels, times):
    """
    repeat shuffle data.
    ----------
    acs : TYPE
        DESCRIPTION.
    sqs : TYPE
        DESCRIPTION.
    labels : TYPE
        DESCRIPTION.
    times : TYPE
        DESCRIPTION.

    Returns
    -------
    out_acs : TYPE
        DESCRIPTION.
    out_sqs : TYPE
        DESCRIPTION.
    out_labels : TYPE
        DESCRIPTION.

    """
    
    out_acs, out_sqs, out_labels = [], [], []
    data = [ [x,y,z] for x,y,z in zip(acs, sqs, labels)]
    data = np.array(data)
    for i in range(times):
        inds = list(range(len(acs)))
        random.shuffle(inds)
        data = data[inds]
    for line in data:
        out_acs.append(line[0])
        out_sqs.append(line[1])
        out_labels.append(line[2])
    out_acs = np.array(out_acs)
    out_sqs = np.array(out_sqs)
    out_labels = np.array(out_labels).astype(int)
    return out_acs, out_sqs, out_labels