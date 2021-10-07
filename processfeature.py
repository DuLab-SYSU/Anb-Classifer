#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 19:13:35 2020

@author: dulab19
"""
import utils
from collections import Counter
import numpy as np
import pickle
from sklearn.model_selection import train_test_split


class ProcessFeature():

    def __init__(self, ac_regions, pssm_path):
        # self.out_feature_names = None
        self.sizes = None
        self.ac_regions = ac_regions
        utils.extract_pssm(pssm_path)
        


    def read_external_files(self, acs):
        with open('./data/cdr_shm_dic.txt', 'rb') as f:
            self.cdr_shm_dic = pickle.load(f)
        with open('./data/v_shm_dic.txt', 'rb') as f:
            self.v_shm_dic = pickle.load(f)
        with open('./data/fr_shm_dic.txt', 'rb') as f:
            self.fr_shm_dic = pickle.load(f)
        self.pgen_dic = utils.read_olga("./data/olga/")
        fats = utils.read_fasta("./data/fasta/all.fasta")
        cur_pssm = utils.read_pssm()
        pssm = {}
        for key, value in fats.items():
            if value in cur_pssm:
                pssm[key] = cur_pssm[value]
        self.pssm = pssm


    def count_features(self, data_list, ac_list):
        import math
        features = []
        num = len(data_list)
        l1 = lambda x: 1 if x == "" else len(x)
        length_lim = {'cdr3':40, 'cdr2':10, 'cdr1':12, 'fr1':25, 'fr2':17,\
                    'fr3':38, 'fr4':11}
        l2 = lambda x, part: x+[0]*(length_lim[part]-len(x)) if len(x) < length_lim[part] else x[: length_lim[part]]
        aaindexs = {}
        index_names = ["helix probability", "hydrophobicity", "isoelectric point",\
                "sheet probability", "steric parameter", "volume", "polarisability"]
        for f_n in index_names:
            aaindexs[f_n] = utils.read_aaindex("./aaindex/{}.txt".format(f_n))
        for i in range(num):
            feature = []
            cdr3, ac = data_list[i], ac_list[i].strip(',') # cdr3 
            # get full v-region sq
            v_region = "".join(self.ac_regions[ac])
            len_v = len(v_region)-len(cdr3)-len(self.ac_regions[ac][-1])
            region_sqs = self.ac_regions[ac]
            cdr_sqs = [region_sqs[1], region_sqs[3], cdr3]
            # shm number of cdr1, cdr2
            shm = []
            for k in range(2):
                shm_num = self.cdr_shm_dic[ac][k]
                shm += [sum(shm_num)/len(cdr_sqs[k])]
            # shm number of fr1, fr2, fr3
            shm += [self.fr_shm_dic[ac][0]/l1(region_sqs[0]), 
                    self.fr_shm_dic[ac][1]/l1(region_sqs[2]), 
                    self.fr_shm_dic[ac][2]/l1(region_sqs[4])]
            feature += shm
            feature.append(self.v_shm_dic[ac]/len_v)  # fr1-fr3 region change rate
            # add aaindex
            tmp = []
            region_order = ['cdr1', 'fr2', 'cdr2', 'fr3', 'cdr3']
            for x, y in zip(region_sqs[1: -1], region_order):
                tmp += utils.count_physiochemical(ac, x, aaindexs, y)
            feature += tmp
            for x in region_sqs[1:-1]:
                if x == "":
                    feature += [0.0]
                else:
                    feature += utils.sq_pi(x)
            feature += [min(-math.log(self.pgen_dic[ac], 10), 50)]
            feature += utils.AAC_PSSM(region_sqs, self.pssm[ac], ac)
            features.append(feature)
            print('\r  cdr3 {} total {}     '.format(i+1,num),end = '')
        return features
