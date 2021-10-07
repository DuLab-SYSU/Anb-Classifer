#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import utils
import numpy as np
import os
# read all sqs and align sqs by imgt scheme


def process_data(name_index):
    imgt_regions = {'cdr1':1,'cdr2':3,'cdr3':5,'fr1':0,'fr2':2,'fr3':4,'fr4':6}
    ac_regions = {}
    vj = {}
    CDR3 = {}
    all_acs = []
    all_cdr3s = []
    unqualifyed_ac = []
    for i in name_index.keys():
        tmp = []
        name = name_index[i]
        acs, v_gene, j_gene = utils.read_as_list("./data/{}_vj.csv".format(name), 3)
        if len(v_gene) != len(j_gene):
            print(name, 'lost part of genes')
        # store AC numbers, v genes, j genes
        for x, y, z in zip(acs, v_gene, j_gene):
            vj[x] = y.split('*')[0] + '|' + z.split('*')[0]
        CDR3[name+'_ac'], CDR3[name+"_cdr3"] = utils.read_as_list("./data/{}_cdr3.csv".format(name), 2)
        # remove redundant cdr3s
        nonrd_ac, nonrd_sq, _ = remove_redundant(CDR3[name+'_ac'], CDR3[name+"_cdr3"])
        del _
        utils.write_out_file("./data/{}_nonrd_ac.csv".format(name), nonrd_ac, "acs:")
        # when individual analysis lose 104 'C', add 'C' to CDR3 sq
        for sq in CDR3[name+'_cdr3']:
            if sq[0] == 'C':
                tmp += [sq[1:]]  # remove 'C' at the beginning of cdr3
            else:
                tmp += [sq]
        # update cdr3 sqs
        CDR3[name+'_cdr3'] = tmp
        # record sqs from other regions
        for rg in imgt_regions.keys():
            acs, sqs = utils.read_as_list("./data/{}_{}.csv".format(name, rg), 2)
            for x, y in zip(acs, sqs):
                if ac_regions.get(x) is None:
                    ac_regions[x] = [""]*7
                ac_regions[x][imgt_regions[rg]] = y
        all_acs += CDR3[name + '_ac']
        all_cdr3s += CDR3[name + '_cdr3']
    all_cdr1s = [ac_regions.get(ac)[1] for ac in all_acs]
    all_cdr2s = [ac_regions.get(ac)[3] for ac in all_acs]
    # align cdr1 cdr2 sqs
    _, _, e_cdr1_ac = utils.align_cdrs(all_acs, all_cdr1s, mode='cdr1', mark='full')
    _, _, e_cdr2_ac = utils.align_cdrs(all_acs, all_cdr2s, mode='cdr2', mark='full')
    del all_cdr1s
    del all_cdr2s
    # align cdr3 sq
    _, _, e_cdr3_ac = utils.align_cdrs(all_acs, all_cdr3s, 'cdr3', 'full')
    unqualifyed_ac += e_cdr1_ac + e_cdr2_ac + e_cdr3_ac
    return ac_regions, CDR3, unqualifyed_ac, vj


def remove_unqualified(name_index, CDR3, unqualifyed_ac, mark=False):
    n_acs, n_sqs, n_labels = [], [], []  # nonredundant data
    sizes = []
    b_list = []
    b_list += unqualifyed_ac
    # l1 = lambda x, y, z, n: y.get(x)[n] if y.get(x) != None else ['-']*z
    for i in name_index.keys():
        # remove redundant sq and ac for each subset
        name = name_index[i]
        cdr3s = CDR3[name+'_cdr3']
        acs = CDR3[name+'_ac']
        # non-redundant acs
        nonrd_ac = utils.read_as_list("./data/{}_nonrd_ac.csv".format(name), 1)
        tmp = set(b_list)  # list to sest
        # remove ac in black list
        nonrd_ac = [ac for ac in nonrd_ac if ac not in tmp]
        tmp = set(nonrd_ac)
        nonrd_cdr3 = [y for x, y in zip(acs, cdr3s) if x in tmp]
        num = len(nonrd_cdr3)
        if num != len(tmp):
            print("Error! ac number != cdr3 number")
            print("Please check your data")
            break
        sizes.append(num)  
        n_acs += nonrd_ac
        n_sqs += nonrd_cdr3
        if mark is True:
            n_labels += [i]*num
    if mark is True:
        return n_acs, n_sqs, n_labels
    return n_acs, n_sqs


def remove_redundant(n_acs, n_sqs, n_labels=None):
    if n_labels is None:
        n_labels = [0]*len(n_acs)
    sq_set = set(n_sqs)
    if len(sq_set) != len(n_acs):
        print("redundant sq exists in dataset", len(sq_set))
        tmp_acs, tmp_sqs, tmp_labels = [], [], []
        for x, y, z in zip(n_sqs, n_acs, n_labels):
            if x in sq_set:
                sq_set.remove(x)
                tmp_acs.append(y)
                tmp_sqs.append(x)
                tmp_labels.append(z)
        n_acs = np.array(tmp_acs)
        n_sqs = np.array(tmp_sqs)
        n_labels = np.array(tmp_labels)
    else:
        print("no overlap between different specific abs")
    n_acs, n_sqs, n_labels = utils.shuffle_data(n_acs, n_sqs, n_labels, 5000)
    return n_acs, n_sqs, n_labels
