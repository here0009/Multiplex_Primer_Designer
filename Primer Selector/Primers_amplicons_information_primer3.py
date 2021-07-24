
# coding: utf-8

# # Primers Amplicons Interaction
# the aim of this script is to check the interactions between primers and amplicons, the basic idea is :
# **12 or < 12 bases between the 3’-end and any amplicon sequence**
# 
# For most of the primers used in viruses detection is ambigous primer, so I should modify to suite this case.
# 
# ## Moudles
# 
# ### Primers and Consensus Sequences to Amplicons
# **Input Files :**
# 1. Primers
# 2. Consensus Sequences
# 
# **Process :**
# 1. Find the postion of primers in consensus seq, the module of string manipulation in Python can be used.
# 2. If can not find primers in consensus seq, or there is no pair of primer in the given consenus seq, show a warning message.
# 
# **Output Files :**
# 1. The Amplicon Sequences
# 2. The position of amplicons in consensus sequences
# 
# 
# ### Import the needed module, and read the test primers into memory
# using primer3-py replace mpprimer for the calculation of primer tm and delta # G, it can also calculate hairpin and homodimer structure, NOT done yet. edited at 20170913
import random
import os
import sys
import re
import itertools
import csv
import pandas as pd
import numpy as np
from collections import defaultdict
import operator
import copy
import TmDeltaG as DG
import math
import ThermodynamicsParameters as TP
import primer3

#Parameters to calculate delatG
mono = 50 #mono-valent ion concentration, mM
diva = 1.5 #divalent ion concentraion, mM
oligo = 50 #oligo concentraion to calculate Tm
dntp = 0.25 #dntp concentraction, mM
degenerate_max_tm = 68 #The max Tm of degenerate primers, Celsius Degree
degenerate_min_tm = 42 #The min Tm of degenerate primers, Celsius Degree
max_degenerate_number = 24 #The max degeneration number of primer


ATGC = {'A','T','G','C'}
IUPAC_DICT = {'R':{'A','G'}, 'Y':{'C','T'},'S':{'G','C'},'W':{'A','T'},
              'K':{'G','T'},'M':{'A','C'},'B':{'C','G','T'},'D':{'A','G','T'},
              'H':{"A","C","T"},'V':{"A","C","G"},'N':{"A","T","C","G"},
             "X":{"A","T","G","C"}, 'A':{'A'}, 'G':{'G'}, 'T':{'T'}, 'C':{'C'}}
class Primer:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        # self.tm = primer3.calcTm(primer_seq)
        self.tm = DG.calTm(primer_seq, mono, diva, oligo, dntp)
        # self.self_hairpin = primer3.calcHairpin(primer_seq)
        self.len = len(sequence)

class DegeneratePrimer:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.atcg_sequence = seq_conversion(self.sequence)
        self.len = len(sequence)
        self.degenerate_number = len(self.atcg_sequence)
    def thermol_result(self): #return the tm range of primers in the degenerate primer
        self.max_tm = 0
        self.min_tm = 1000
        # self.hairpin_structure_number = 0
        # hairpin_structure_number = 0
        for seq in self.atcg_sequence:
            # hairpin = primer3.calcHairpin(seq)
            tm = DG.calTm(seq, mono, diva, oligo, dntp)
            if tm>=self.max_tm:
                self.max_tm = tm
            if tm<=self.min_tm:
                self.min_tm = tm
            # if hairpin.structure_found == True:
            #     self.hairpin_structure_number+=1
        return self.max_tm,self.min_tm

class Primer_Ref_Match:
    def __init__(self, primer_seq, ref_seq, match_score, strand, pos_start, pos_end):
        self.primer_seq = primer_seq
        self.ref_seq = ref_seq
        self.match_score = match_score
        self.strand = strand
        self.pos_start = pos_start
        self.pos_end = pos_end

def seq_split(seq):
    seq_atgc = ''
    seq_other = ''
    seq = str.upper(seq)
    for i in seq:
        if i in ATGC:
            seq_atgc += i
        elif i in IUPAC_DICT:
            seq_atgc += '?'
            seq_other += i
        else:
            seq_atgc += '!'
            print("Warning :There is an unrecognizable characters %s." %i)
    return seq_atgc,seq_other

def seq_generation(seq_atgc, seq_other):
    seqs_converted = []
    for amb_seq in itertools.product(*[IUPAC_DICT[i] for i in seq_other]):
    #generate the all combinations of ambiguous letters
        seq_converted = '' 
        amb_seq = ''.join(list(amb_seq))
        amb_index = 0
        for seq_index, letter in enumerate(seq_atgc):
            if letter == '?':
                seq_converted += amb_seq[amb_index]
                amb_index += 1
            else:
                seq_converted += seq_atgc[seq_index]
        seqs_converted.append(seq_converted)
    return tuple(seqs_converted) #change the returned atcg sequences to tuple

def seq_conversion(seq):
    seq_atgc, seq_other = seq_split(seq)
    primers = seq_generation(seq_atgc,seq_other)
    return primers


def read_fasta_file(file_name) :
    '''
    Function to read fasta file, store the names and sequences in a dictionary
    '''
    fasta_dict = dict() # a dict store the primers
    number = 1
    with open (file_name) as fhand :
        for line in fhand:
            if line.startswith('>'):
                name = line.strip()[1:]
                seq = primer_seq = fhand.readline().strip()
                if name in fasta_dict:
                    print("Already exsiting fasta name of {}".format(name))
                    name = name + str(number)
                    number+=1
                    print("Change it to the name of {}".format(name))
                fasta_dict[name] = seq
    return fasta_dict

def print_dict(test_dict):
    '''
    Function to print a dictionary, so it is easy to read
    '''
    for key, value in test_dict.items():
        print('{}\t{}'.format(key, value))



def cpl(seq):
    '''
    return a complementary sequence of a given sequence
    '''
    complementary_dict ={'A':'T', 'T':'A', 'G':'C', 'C':'G','Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'H':'D', 'V':'B', 'B':'V', 'X':'X', 'N':'N'}
    try:
        cpl_seq = [complementary_dict[i] for i in seq]
        cpl_seq = ''.join(cpl_seq)
        return cpl_seq
    except KeyError:
        print('undefiend character in' , seq)
        pass

def rev_cpl(seq):
    '''
    return a reverse complementary sequence of a given sequence
    '''
    cpl_seq = cpl(seq)
    rev_cpl_seq = cpl_seq[::-1]
    return rev_cpl_seq


def primer_information_check(primer_name, primer_seq,file_hand):
    flag = True
    msg_deg = '{:35}\t{:40}\t{:^10d}\t{:^10d}\t{:^5.2f}\t{:^5.2f}\t{:^10d}\n' #for output of degenerate primers
    primer = DegeneratePrimer(primer_name, primer_seq)
    max_tm, min_tm = primer.thermol_result()
    if max_tm >degenerate_max_tm or min_tm<degenerate_min_tm or primer.degenerate_number>max_degenerate_number :
        flag = False
    file_hand.write(msg_deg.format(primer_name,primer_seq,primer.len,primer.degenerate_number,max_tm,min_tm, flag))
    return flag
#The find substring process, remember there must be two and only two 
#primers should be found in every amplicon
 
# Rule No. 6
# 12 or < 12 bases between the 3’-end and any amplicon sequence
# I think this is simple, just return the 3’-end of a primer and check it with the whole amplicon seqs
# 
# 
# ## output file
# 
# 1. the amplicons 
# 2. the positions of primers on consensus sequences


def match_score_letter(primer_letter, amplicon_letter):
    '''
    Return match score of a single letter
    '''
    primer_letter_converted = IUPAC_DICT[str.upper(primer_letter)]
    amplicon_letter_converted = IUPAC_DICT[str.upper(amplicon_letter)]
    letter_in_common = primer_letter_converted & amplicon_letter_converted #两者的共同碱基
    letter_only_in_primer = primer_letter_converted - amplicon_letter_converted #primer多出amplicon的部分
    letter_only_in_amplicon = amplicon_letter_converted - primer_letter_converted #amplicon多出primer的部分
    if  len(letter_in_common) == 0: #无共同碱基，match值为0
        match_score_letter = 0
    elif len(letter_only_in_primer) >= 0 and len(letter_only_in_amplicon) == 0 :
    #primer的碱基包括amplicon,match值为1
        match_score_letter =1
    else : #其他情况：amplicon包括primer，或二者有交集,match_score为交集除以amplicon
        match_score_letter = len(primer_letter_converted&amplicon_letter_converted)/len(amplicon_letter_converted)
    return match_score_letter



def match_score_seq(primer_seq, amplicon_sub_seq):
    '''
    Return the match score of two equal length sequence
    '''
    len_seq = len(primer_seq)
    if len_seq != len(amplicon_sub_seq):
        print("two seq length are not equal")
        return
    match_score_seq = 0
    for i in range(len_seq):
        match_score_each_letter = match_score_letter(primer_seq[i],amplicon_sub_seq[i])
        match_score_seq += match_score_each_letter
    return match_score_seq, primer_seq, amplicon_sub_seq




def primer_amplicon_match_score(primer, ampliseq, threshold, strand): 

    '''
    Test the interaction between 3' end of primer and the amplicon, the parameter sense is (f)orward or (r)reverse
    thershold is percent of matches
    '''
    #确保primer长度大于length
    # primer_amplicon_match_score = 0
    primer_ref_match_list = list()
    length = len(primer)
    # if len_primer-length < 0:
    #     print('The thershold is longer than the length of primer')
    #     return 
    #确认序列方向为正向或反向
    if strand == 'forward':
        primer_3_end = primer[-length:]
    elif strand == 'reverse':
        primer_3_end = primer[:length]
    else :
        print("Wrong Parameter, Please give the strand information of primer")
        return
    
    for i in range(len(ampliseq)-length+1):
        match_len, primer_seq, ref_seq = match_score_seq(primer_3_end, ampliseq[i:i+length])
        match_score = match_len / length
        if match_score >= threshold:
            # fhand_primer_positon.write("headline.format(match_score_sub_seq[1],strand, match_score_sub_seq[2]),match_score_sub_seq[0], i, i+length)
            primer_pos_info = Primer_Ref_Match(primer_seq,ref_seq, match_score, strand, i, i+length)
            primer_ref_match_list.append(primer_pos_info)
            # primer_amplicon_match_score +=1

    return primer_ref_match_list


def primer_amplicon_check(primer_dict, consensus_seq_dict, file_tag):
    '''
    Check the interaction between the primers and amplicons, \
    the primers and amplicons are stored in different dicionary.
    '''
    primer_amplicon_match_dict = dict()
    for primer_name in primer_dict.keys():
        primer_amplicon_match_dict[primer_name] = 0

    # primer_postion_file_name = primer_file.split(".")[0]+'_'+ "position.txt"
    primer_postion_file_name = ".\\{0}\\{0}_position.txt".format(file_tag)
    fhand_primer_positon = open(primer_postion_file_name, 'w')
    primer_pos_headline = '{:35}\t{:35}\t{:35}\t{:35}\t{:^5}\t{:^5}\t{:^10}\t{:^10}\t{:^10}\n'
    fhand_primer_positon.write(primer_pos_headline.format('Primer_Name', 'Ref_Name' ,'Primer_Seq','Ref_Seq','Score','Strand','Pos_Start','Pos_End','Identical'))
    for ref_name, ref_seq in consensus_seq_dict.items():
        for primer_name, primer_seq in primer_dict.items():
            
            threshold = 0.75 #if the score is larger than 0.8*length, return true
            rev_cpl_primer = rev_cpl(primer_seq)
            match_list_f = primer_amplicon_match_score(primer_seq, ref_seq,  threshold, 'forward')
            match_list_r = primer_amplicon_match_score(rev_cpl_primer, ref_seq,threshold, 'reverse')
            match_list = match_list_f + match_list_r
            #正向引物与amplcion的match score
            if len(match_list) > 0:
                for i in range(len(match_list)):
                    if match_list[i].primer_seq == match_list[i].ref_seq:
                        identical = True
                    else:
                        identical = False
                    fhand_primer_positon.write(primer_pos_headline.format(primer_name, ref_name, match_list[i].primer_seq, match_list[i].ref_seq, match_list[i].match_score, match_list[i].strand, match_list[i].pos_start, match_list[i].pos_end, identical))

            match_score = len(match_list)
            if match_score >= 1:
                primer_amplicon_match_dict[primer_name] = primer_amplicon_match_dict[primer_name] + match_score

    
    print_dict(primer_amplicon_match_dict)
    return primer_amplicon_match_dict



#read the primers in memory
#primer_file = input("Enter the file name of primers :")


def primer_info_doc(file_tag):
    #primer information check
    # outputfile_primer_information = primer_file.split(".")[0] + "_primer_information.txt"
    outputfile_primer_information = ".\\{0}\\{0}_primer_information.txt".format(file_tag)
    fhand_primer_information = open(outputfile_primer_information,'w')
    headline = '{:35}\t{:40}\t{:^10}\t{:^10}\t{:^5}\t{:^5}\t{:^10}\n'
    fhand_primer_information.write(headline.format('Name','Sequence','Len','DegNumber','MaxTm','MinTm','Flag'))
    for primer_name, primer_seq in primer_dict.items():
        primer_information_check(primer_name, primer_seq, fhand_primer_information)
    fhand_primer_information.close()

def primer_pos_doc(file_tag):
    # output_file_name = primer_file.split(".")[0]+"_info_report.fasta"
    output_file_name = ".\\{0}\\{0}_info_report.fasta".format(file_tag)
    sys.stdout=open(output_file_name,"w")
    print("The match score between the primer and amplicons are :")
    print("===============================================")
    primer_amplicon_match_dict = primer_amplicon_check(primer_dict, consensus_seq_dict, file_tag)
    print("These primers got more than one match or no match:")
    print("===============================================")
    for primer, match_score in primer_amplicon_match_dict.items():
        if match_score != 1 :
            print('{} \t {}'.format(primer, match_score) )
    sys.stdout.close()


primer_file = sys.argv[1]
primer_dict = read_fasta_file(primer_file) 

#read the consensus sequences in memory
consensus_seq_file = sys.argv[2]
consensus_seq_dict = read_fasta_file(consensus_seq_file)
file_tag = primer_file.split(".")[0]
try:
    os.stat(".\\{}".format(file_tag))
except:
    os.mkdir(".\\{}".format(file_tag))
primer_info_doc(file_tag)
primer_pos_doc(file_tag)
