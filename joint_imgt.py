# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 15:28:25 2020

@author: ainer
"""
import os
import utils as ut
import re
import argparse

def joint_file(filename, ac_list, export_path):
    """
    Parameters
    ----------
    filename : str
    ac_list : list
    export_path : str

    Returns
    -------
    None.
    conbine individual files into a file

    """
    if os.path.exists(export_path) == True:
        os.remove(export_path)
    for ac in ac_list:
        path = filename + ac.replace(' ', '_')
        if os.path.isfile(path):
            content = []
            with open(path, 'r') as f:
                 for line in f:
                     content.append(line)
            with open(export_path, 'a+') as f:
                f.writelines(content)

def check_part_region(path, mode='CDR3'):
    mark = 0
    mark1 = 0
    ac = []
    region = []
    v_gene = []
    j_gene = []
    regex = re.compile("[ ]+[A-Z]+")
    regex1 = re.compile("IGHV[\d][-]?[\d]{0,2}[-]?[\d]{0,2}[*]?[\d]{0,2}")
    regex2 = re.compile("IGHJ[\d][*]?[\d]{0,2}")
    regex3 = re.compile("[0-9]{1,3}[.]{2,3}[0-9]{1,3}")
    with open(path,"r") as f:
        for line in f:
            if "---------------------" in line:
                mark = 1
                sq = ""
            elif mark == 1 and "Sequence number" in line:
                if "no results" in line:
                    mark = 0
                else:
                    AC = line.split(": ")[1]
                    AC = AC.strip('\n').split(' ')[0]
            elif mark == 1 and ("V-GENE and allele;" in line) and ("IMGT Label" not in line):
                search1 = re.search(regex1, line)
                if search1:
                    v_gene.append(search1.group(0))
            elif mark == 1 and ("J-GENE and allele;" in line) and ("IMGT Label" not in line):
                search2 = re.search(regex2, line)
                if search2:
                    j_gene.append(search2.group(0))    
            elif mark == 1 and "V-D-J-REGION " in line:
                mark1 = 1
            elif mark1 == 1 and "/CDR_length=" in line and "X" not in line and ".0" not in line:  # CDR complex
                mark1 = 2
            elif mark1 == 2 and "{}-IMGT".format(mode) in line:
                mark1 = 3
            elif mark1 == 3 and "/translation" in line:
                mark1 = 4
            elif mark1 == 4 and regex.match(line):  # AA
                tmp = line.strip().split()[0]
                sq += tmp
            elif mark1 == 4 and regex3.search(line):  # xx..xx
                region.append(sq)
                ac.append(AC)
                mark1 = 0
            elif mark1 == 4 and mode == 'FR4':
                region.append(sq)
                ac.append(AC)
                mark1 = 0
            else:
                pass
    return ac, region, v_gene, j_gene

def read_summary(folder_dir, save_path, ag=None):
    num = 0
    ac_list = []
    path = r"{}/1_Summary.txt".format(folder_dir)
    with open(path,'r') as f:
        for line in f:
            line= line.strip('\n')
            if "No results" in line or "unproductive" in line:
                pass
            else:
                if "IGH" in line:
                    tmp = line.split('\t')
                    ac_list.append(tmp[1]+"_"+tmp[0])
                    num += 1
    ut.write_out_file(save_path, ac_list, "{} productive IGH:".format(ag))
    return num, ac_list

def get_sq(path):
    mark = 0
    mark1 = 0
    ac_list = []
    sqs = []
    regex = re.compile("[ ]+[A-Z]+")
    error_ac = []
    with open(path,"r") as f:
        for line in f:
            if "---------------------" in line:
                mark = 1
            elif mark == 1 and "Sequence number" in line:
                if "no results" in line:
                    mark = 0
                else:
                    AC = line.split(": ")[1]
                    AC = AC.strip('\n').split(' ')[0]
                    sq = ""
            elif mark == 1 and "V-D-J-REGION " in line:
                mark1 = 1
            elif mark1 == 1 and "/translation" in line:
                mark1 = 2
            elif mark1 == 2 and regex.match(line):
                sq += line.strip().split()[0]
            elif "V-REGION" in line and mark1 == 2:
                mark1 = 0
                if ('*' not in sq) and ('X' not in sq):
                    sqs.append(sq)
                    ac_list.append(AC)
                else:
                    error_ac.append(AC)
            else:
                pass
    return ac_list, sqs, error_ac

def split_fasta(path, ac, sq_list, max_num):
    total = len(ac)
    acs = [ac[k: k + max_num] for k in range(0, total, max_num)]
    sqs = [sq_list[k: k + max_num] for k in range(0, total, max_num)]
    sub_num = len(acs)
    for i in range(sub_num):
        ut.write_fasta(path + str(i) + '.fasta', acs[i], sqs[i], "")
    return None

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, help="IMGT HighV-QUEST result folder")
parser.add_argument('-o', type=str, default='./original_model/data', help="output path for cdrs and v-region sqs")
args = parser.parse_args()
error_ac = []
# output sqs of each region
if '/' in args.i:
    # linux path
    name = args.i.strip('/').split('/')[-1]
elif "\\" in args.i:
    # windows path
    name = args.i.strip('\\').split('\\')[-1]
# read summary file in folder parser.path
num, ac_list = read_summary(args.i, os.path.join(args.o, '{}_h.csv'.format(name)))
# individual files path
path = os.path.join(args.i, "IMGT_HighV-QUEST_individual_files_folder/")
# output path for integrated individual files
igh_path = os.path.join(args.i, "IGH.txt")
# joint individual files
joint_file(path , ac_list, igh_path)
# 7 regions
imgt_regions = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'FR4', 'CDR3']
# get v-region protein sq
ac, sq_list, error_ac = get_sq(igh_path)
for x,y in zip(ac, sq_list):
    if (len(y) <= 100):  # the lower limit for ab sq length
        error_ac.append(x)
print("unqualified ACs:", len(error_ac))
# get cdr3 sqs and v gene j gene
ac, region, v_ge, j_ge =  check_part_region(igh_path, 'CDR3')
# 
ac = list(set(ac)-set(error_ac))
print("number of ACs:", len(ac))
tmp = []  # store qualified acs
for ele1 in ac:
    for ele2 in ac_list:
        if ele1 in ele2:
            tmp.append(ele2)
            break
# remove unqualified sqs in IGH
joint_file(path, tmp, igh_path)
for rg in imgt_regions:
    ac, region, v_ge, j_ge = check_part_region(igh_path, rg)
    ac = [ele.replace(",", '_') for ele in ac]
    if len(ac) == len(region):
        ut.write_out_file(
            os.path.join(args.o, "{}_{}.csv".format(name, rg.lower())),\
            ac, "{} ac:".format(name)
            )
        ut.write_out_file(
            os.path.join(args.o, "{}_{}.csv".format(name, rg.lower())),\
            region, "{} {}:".format(name, rg), 'a+'
            )
        if rg == 'CDR3':
            print("unique cdr3:", len(set(region)))
            ut.write_out_file(
                os.path.join(args.o, "{}_vj.csv".format(name)),\
                ac, "{} ac:".format(name)
                )
            ut.write_out_file(
                os.path.join(args.o, "{}_vj.csv".format(name)),
                v_ge, "{} gene:".format(name), 'a'
                )
            ut.write_out_file(
                os.path.join(args.o, "{}_vj.csv".format(name)),
                j_ge, "{} j gene:".format(name), 'a'
                )
# get and output v-region sqs
ac, sq_list, error_ac = get_sq(igh_path)
ac = [ele.replace(",", '_') for ele in ac]
# check fasta folder
if os.path.exists(os.path.join(args.o, "fasta/")) is False:
    os.mkdir(os.path.join(args.o, "fasta/"))
ut.write_fasta(
                os.path.join("./data", "fasta/{}_p.fasta".format(name)),
                ac, sq_list, "", "w"
                )