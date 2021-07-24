#!/usr/bin/python
#coding:utf-8
#modification : merge the function of primer check and primer amplicon check. The primers are forward and reverse primers are allocated to different group.
# Add hairpinCheck() and homodimerCheck() function to probe, if the  tm of the hairpin structure  is greater than 50, hairpin_number +=1, if hairpin number great than 25% of the total degenerate number, the probe is unqualified.
# If the delta G value of  homodimer structure is less than -9kcal*mol-1, the probe is unqualified
# Import the needed module, and read the test primers into memory
import random
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

#Parameters to calculate primer interactions
rule_1_len = 4 #3' ends len
rule_1_cutoff = 4 #3' ends cutoff score
rule_2_len = 8 #3' ends len(with one mismatch)
rule_2_cutoff = 7 #3' ends cutoff score(with one mismatch)
rule_3_len = 10 #primer len
rule_3_cutoff = 10 #primer cutoff score
rule_4_len = 13 #primer len(with one mismatch)
rule_4_cutoff = 12 #primer cutoff score(with one mismatch)
rule_5_len = 14 #match percent calculation min len
rule_5_cutoff = 0.75 #match percent cutoff score
rule_6_len = 9 #alignment score len (min length)
rule_6_cutoff = 9 #alignment cutoff score
rule_7_len = 7 #delta G calculation min len
rule_7_cutoff = -9 #delta G cutoff score, kcal/mol
rule_7_cutoff = -1*rule_7_cutoff #change cut-off score to positive value


amplicon_primer_check_len = 14 #The 3' ends of primer that will be checked with amplicon sequence
amplicon_primer_cutoff = 12 #The score of primer-amplicon interaction
amplicon_len_max = 250 #The maximum length of amplicon
amplicon_len_min = 70 #The minimum length of amplicon
amplicon_lendiff_min = 0 #The minimum length differences between amplicons

confliction_table_threshold = 1 #The minimum confliction table, default value is 1.
elminate_self_interaction = 1 #Wether eliminate primer interaction, default is 1.

threshold_position = 0.9 #the threshold for postion check, if primer can overlap the ref seq to this ratio, it will not used as a candidate primer
#Parameters for the probe
primer_probe_gap = 1 #the nucleotides between primer and probe, default set to be 1.
probe_len_min = 23
probe_len_max = 30
probe_tm_min = 60
probe_tm_max = 80
probe_init_noG_check = True #If set to true, the 5' nucleotide of the primer should not be G
probe_more_CthanG = True #If ste to true, select the probe has more C than G
probe_max_deg_threshold = 32 #The maximum degeneration number allowed for probe
probe_ambigous_letter_threshold = 4 #The max ambigous letter allowed in probe.
hairpin_number_cutoff = 0.25 #The max allowed ratio of  (probe got hairpin struct / degenerate number)
homodimer_number_cutoff = 0.25 #The max allowed ration of (homodimers that got larger deltaG value than homodimer_dg_threshold/degenerate number)
hairpin_tm_threshold = 50 #The max tm allowed for a hairpin structure
homodimer_dg_threshold = -9000 #The min deltaG value allowed for homodimer structure, cal/mol

#Dictionary of ATGC and Degenerate Bases
ATGC = ['A', 'T', 'G', 'C']
IUPAC_DICT = {'R':{'A','G'}, 'Y':{'C','T'},'S':{'G','C'},'W':{'A','T'},
              'K':{'G','T'},'M':{'A','C'},'B':{'C','G','T'},'D':{'A','G','T'},
              'H':{"A","C","T"},'V':{"A","C","G"},'N':{"A","T","C","G"},
             "X":{"A","T","G","C"}, 'A':{'A'}, 'G':{'G'}, 'T':{'T'}, 'C':{'C'}}

#dictonary and  functions for probe design - begin
degenerate_number_dict = {'R':2, 'Y':2,'S':2,'W':2,'K':2,'M':2,'B':3,'D':3,'H':3,'V':3,'N':4,"X":4}

def seq_ambigous_number(seq):
    '''
    return the ambigous letter number in a given sequence
    '''
    ambigous_number = 0
    for letter in seq:
        if letter not in ATGC:
            ambigous_number += 1
    return ambigous_number

def probe_init_G(seq):
    '''
    Check if the initiation of the probe is G or ambigous letter inculde G
    return True if the initiation letter is G
    '''
    init_letter = seq[0]
    if 'G' in IUPAC_DICT[init_letter]:
        return True
    else:
        return False

def seq_degenerate_number(seq):
    '''
    return the degenerate_number of the given seq
    '''
    degenerate_number = 1
    for letter in seq:
        if letter in degenerate_number_dict:
            degenerate_number = degenerate_number_dict[letter]*degenerate_number

    return degenerate_number

#dictonary and  functions for probe design - end

def read_fasta_file(file_name) :
    '''
    Function to read fasta file, store the names and sequences in a dictionary
    '''
    fasta_dict = dict() # a dict store the primers
    with open (file_name) as fhand :
        for line in fhand:
            if line.startswith('>'):
                name = line.strip()[1:]
                seq = primer_seq = fhand.readline().strip()
                fasta_dict[name] = seq
    return fasta_dict

def print_dict(test_dict):
    '''
    Function to print a dictionary, so it is easy to read
    '''
    for key, value in test_dict.items():
        print(key)
        print(value)

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

# ### Funcion *cpl* to return the value of the complementary seq of a seq
def cpl(seq):
    complementary_dict ={'A':'T', 'T':'A', 'G':'C', 'C':'G','Y':'R', 'R':'Y', 'W':'W', 'S':'S', 'K':'M', 'M':'K', 'D':'H', 'H':'D', 'V':'B', 'B':'V', 'X':'X', 'N':'N',' ':' '}
    try:
        cpl_seq = [complementary_dict[i] for i in seq]
        cpl_seq = ''.join(cpl_seq)
        return cpl_seq
    except KeyError:
        print('undefiend character in' , seq)
        pass

def rev_cpl(seq):
    '''
    Return the complementary sequence of a given sequence
    '''
    cpl_seq = cpl(seq)
    rev_cpl_seq = cpl_seq[::-1]
    return rev_cpl_seq

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
        match_score_letter = len(primer_letter_converted&amplicon_letter_converted)/        len(amplicon_letter_converted)
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

def primer_amplicon_match_score(primer, ampliseq, length, threshold, strand): 
    #根据文献primer 3'端与amplicon最多有12个碱基一样
    #因为有简并引物，引入一个thershold值控制match_score,先尝试length为13，threshold为11
    '''
    Test the interaction between 3' end of primer and the amplicon, the parameter sense is (f)orward or (r)reverse
    '''
    #确保primer长度大于length
    primer_amplicon_match_score = 0
    len_primer = len(primer)
    if len_primer-length <= 0:
        print('The thershold is longer than the length of primer')
        return 
    #确认序列方向为正向或反向
    if strand == 'forward':
        primer_3_end = primer[-length:]
    elif strand == 'reverse':
        primer_3_end = primer[:length]
    else :
        print("Wrong Parameter, Please give the strand information of primer")
        return
    
    for i in range(len(ampliseq)-length):
        match_score_sub_seq = match_score_seq(primer_3_end, ampliseq[i:i+length])
        if match_score_sub_seq[0] >= threshold:
            #print("Found one match between the primer and amplicon above threshold, the match score is {}.\nThe primer 3 end seq is {} it is on {} strand.\nThe amplicon seq is {}".format(match_score_sub_seq[0],match_score_sub_seq[1],strand, match_score_sub_seq[2]))
            primer_amplicon_match_score +=1
    return primer_amplicon_match_score

def primer_amplicon_check(primer_dict, ref_dict, length, threshold):
    '''
    Check the interaction between the primers and amplicons, \
    the primers and amplicons are stored in different dicionary.
    '''
    primer_amplicon_match_dict = dict()
    for amplicon_name, amplicon_seq in ref_dict.items():
        for primer_name, primer_seq in primer_dict.items():
            rev_cpl_primer = rev_cpl(primer_seq)
            match_score_f = primer_amplicon_match_score(primer_seq, amplicon_seq, length,threshold, 'forward') 
            #正向引物与amplcion的match score
            match_score_r = primer_amplicon_match_score(rev_cpl_primer, amplicon_seq, length,threshold, 'reverse')
            #反向引物与amplicon的match score,反向引物应该为引物前半段与序列的比较
            match_score = match_score_f + match_score_r
            if match_score >=1:
                # print("Find {} matches of primer {} in amplicon {}".format(match_score,primer_name, amplicon_name))
                # print('\n')
                primer_amplicon_match_dict[primer_name] = primer_amplicon_match_dict.get(primer_name,0)+match_score

    return primer_amplicon_match_dict

def amplicon_length_check(amplicon_min,amplicon_max,amplicon_diff,primer_tuple, position_dict, primer_amplicon_match_dict):
    '''
    Return True if the amplicon len is between min and max and the lenght diffrent between amplicons is above threshold.
    '''
    flag = True
    amplicon_len_list = []
    for comb in itertools.combinations(primer_tuple,2):
        if position_dict[comb[0]].ref_name == position_dict[comb[1]].ref_name: #the same organism
            if primer_amplicon_match_dict[comb[0]] + primer_amplicon_match_dict[comb[1]] >=4:
                flag = False
            if position_dict[comb[0]].strand == 'reverse':
                reverse_position = position_dict[comb[0]].pos_end
                forward_position = position_dict[comb[1]].pos_start
            else:
                reverse_position = position_dict[comb[1]].pos_end
                forward_position = position_dict[comb[0]].pos_start
            amplicon_len = reverse_position - forward_position
            amplicon_len_list.append(amplicon_len)
            if amplicon_len<amplicon_min or amplicon_len>amplicon_max:
                flag = False
    amplicon_len_list = sorted(amplicon_len_list)
    for i in range(len(amplicon_len_list)-1):
        if amplicon_len_list[i+1] - amplicon_len_list[i] < amplicon_diff: #amplicon_diff is the minimum lenghth diff between different amplicons
            flag = False
    return flag,amplicon_len_list

class Primer:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        # self.tm = primer3.calcTm(primer_seq)
        self.tm = DG.calTm(primer_seq, mono, diva, oligo, dntp)
        # self.self_hairpin = primer3.calcHairpin(primer_seq)
        self.len = len(sequence)
    def self_hairpin(self):
        return primer3.calcHairpin(primer_seq)

class DegeneratePrimer:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.len = len(sequence)
        # self.degenerate_number = len(self.atcg_sequence)
    def thermol_result(self): #return the tm range of primers in the degenerate primer
        self.max_tm = 0
        self.min_tm = 1000
        for seq in seq_conversion(self.sequence):
            # hairpin = primer3.calcHairpin(seq)
            tm = DG.calTm(seq, mono, diva, oligo, dntp)
            if tm>=self.max_tm:
                self.max_tm = tm
            if tm<=self.min_tm:
                self.min_tm = tm
            # if hairpin.structure_found == True:
            #     self.hairpin_structure_number+=1
        return self.max_tm,self.min_tm
    def atcg_sequence(self):
        return seq_conversion(self.sequence)
    def degenerate_number(self):
        return seq_degenerate_number(self.sequence)
    def hairpin_check(self):
        '''
        Check the hairpin structure of the given degenerate primer
        '''
        # max_hairpin_dg = 0
        # min_hairpin_dg = 0
        max_hairpin_tm = 0
        hairpin_structure_number = 0 #the number primer got hairpin tm larger than hairpin_tm_threshold
        for seq in self.atcg_sequence():
            hairpin = primer3.calcHairpin(seq)
            if hairpin.tm >= hairpin_tm_threshold:
                hairpin_structure_number+=1
                if hairpin.tm >= max_hairpin_tm:
                    max_hairpin_tm = hairpin.tm
                # if hairpin.dg < min_hairpin_dg:
                #     min_hairpin_dg = hairpin.dg
        return hairpin_structure_number, max_hairpin_tm

    def homo_dimer_check(self):
        '''
        Check the homo dimer of the given degenerate primer
        '''
        # max_homo_dimer_tm = 0
        min_homo_dimer_dg = 0
        homodimer_number = 0 #the number of homo_dimer_dg smaller than homodimer_dg_threshold
        for seq in self.atcg_sequence():
            homo_dimer = primer3.calcHomodimer(seq)
            if homo_dimer.dg < min_homo_dimer_dg:
                min_homo_dimer_dg = homo_dimer.dg
            if homo_dimer.dg <= homodimer_dg_threshold:
                homodimer_number+=1
                # if homo_dimer.tm > max_homo_dimer_tm:
                #     max_homo_dimer_tm = homo_dimer.tm
        # return homodimer_number, max_homo_dimer_tm
        return homodimer_number, min_homo_dimer_dg




class Primer_Ref_Match:
    def __init__(self, primer_name, ref_name, primer_seq, ref_seq, match_score, strand, pos_start, pos_end):
        self.primer_name = primer_name
        self.ref_name = ref_name
        self.primer_seq = primer_seq
        self.ref_seq = ref_seq
        self.match_score = match_score
        self.strand = strand
        self.pos_start = pos_start
        self.pos_end = pos_end


def primer_information_check(primer_name, primer_seq,file_hand):
    flag = True
    msg_deg = '{:35}\t{:40}\t{:^10d}\t{:^10d}\t{:^5.2f}\t{:^5.2f}\t{:^10d}\n' #for output of degenerate primers
    primer = DegeneratePrimer(primer_name, primer_seq)
    max_tm, min_tm = primer.thermol_result()
    if max_tm >degenerate_max_tm or min_tm<degenerate_min_tm or primer.degenerate_number()>max_degenerate_number :
        flag = False
    file_hand.write(msg_deg.format(primer_name,primer_seq,primer.len,primer.degenerate_number(),max_tm,min_tm, flag))
    return flag

# In[30]:
def primer_alignment_print(primer_1, primer_2, partial_or_full_len, parameter):
    '''
    A function return the alignment of two primers.
    Input: two primer sequences, 
        partial_or_full_len : means the state of primer_1 and primer_2 alignment
        parameter: 
        if partial_alignment, 
            parameter is the overlap_len between two primers, 
                overlap_len > 0 means the first n character of primer_2 will be represented as ' '.
                overlap_len < 0 means the first n character of primer_1 will be represented as ' '.
        if full_alignment,
            since primer_1 is the longer primer, then parameter is the stagger_value of primer_2. 
            The first n character of primer_2 will be represented as ' '.
    Output: alignment prints of the two primers
    '''
    primer_2 = primer_2[::-1]
    if partial_or_full_len == 'partial':
        #Partial alignment
        overlap_len = parameter
        len_seq = len(primer_1) + len(primer_2) - abs(overlap_len)
        #the initialization of three strings that will be printed
        middle_line = [' ']*len_seq
        #Put the sequence into the corresponding position of the list
        if overlap_len >= 0:
            primer_1_print = primer_1 + (len_seq - len(primer_1))*' '
            primer_2_print = (len_seq - len(primer_2))*' ' + primer_2
        if overlap_len <0:
            primer_2_print = primer_2 + (len_seq - len(primer_2))*' '
            primer_1_print = (len_seq - len(primer_1))*' ' + primer_1
    elif partial_or_full_len == 'full_len':
        #Full Length alignment
        stagger_value = parameter
        len_seq = len(primer_1)
        len_diff = len_seq - len(primer_2)
        middle_line = [' ']*len_seq
        primer_1_print = primer_1
        primer_2_print = ' '*stagger_value + primer_2 + (len_diff-stagger_value)*' '
    else:
        print('Warning : Wrong Input value of partial_or_full_len')
        return
    #Convert the list to string
    for i in range(len_seq):
        if primer_1_print[i] == cpl(primer_2_print[i]):
            middle_line[i] = '|'
    middle_line = ''.join(middle_line)

    print (primer_1_print)
    print (middle_line)
    print (primer_2_print)

# ### Filter of Rule No.1 and No.2, they are usable

# '''
# Rule No.1
# If seq1 and seq2 got 4 complementary seq in 3 end, 
# return True,
# else return False
# '''
# def primer_check(primer1, primer2):
#     if primer1[-4:] == primer2[-4:]
#         return True
#     else:
#         return False


def three_ends_match_score(primer_1, primer_2, length):
    '''
    Useful
    Modify primer_check
    given seq1, seq2 and length 3 parameters
    return the complemantry score of seq1 and seq2
    This can be applied to Rule1, and 2
    '''
    len_1 = len(primer_1)
    len_2 = len(primer_2)
    if len_2 > len_1: # make primer1 to be the longer one
        primer_1, primer_2 = primer_2, primer_1
    score = 0
    seq1 = primer_1[-1*length:] #seq1 is the 3' ends of seq1
    seq2 = primer_2[:length] #modify, the seq2 is 5' ends of the cpl_rev(seq)
    for i in range(length):
        if seq1[i] == seq2[i]:
            score +=1          
    return score


# ### Filter of Rule No.3 and No.4 Version 2 
# ### Done
# 
# Main Idea: 
# 1. one funtion match_score(seq1, seq2, panelty) return the match score of seq1 and seq2. seq1 and seq2 are equal length.
# 2. primer1 and primer2 have limited match styles, use match_score() to get the largest match score of primer1 and primer2
# 3. may be these function need to be modified, so it is suitable for insertion and deletion(do it later)

# In[31]:


def match_score(seq1, seq2, max_mismatch):
    '''
    Return the match_score of two equal length seq1 and seq2, allow mismatch
    Useful
    '''
    len_seq = len(seq1)
    if len_seq != len(seq2):
        print("two seq length are not equal")
        return
    match = 0
    mis_match = 0
    match_score = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            match+=1
            if match > match_score:
                match_score = match
        else:
            mis_match+=1
            if mis_match > max_mismatch:
                match = 0
                mis_match = 0
            continue
    return match_score

def alignment_score(seq1, seq2):

    '''
    A function given back the max alignment score of two equal length sequence
    match : score + 1
    mismatch open : score - 2
    mismatch extension : score -1
    '''
    max_alginment_score = 0
    temp_alginment_score = 0
    gap_open_flag = False #a parameter record the gap opening status
    len_seq = len(seq1)
    if len_seq != len(seq2):
        print("two seq length are not equal")
        return
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            temp_alginment_score+=1
            if temp_alginment_score > max_alginment_score:
                max_alginment_score = temp_alginment_score
                gap_open_flag = True
        else:
            if gap_open_flag == True : #判断是否出现gap open的状况，如果之前一个是配对，则为出现
                temp_alginment_score-=1 #score-1, 如果之前一个为gap，则未出现。
                gap_open_flag = False
            temp_alginment_score-=1 #出现gap, score -1 
            if temp_alginment_score < 0:
                temp_alginment_score = 0 #最小值为0，不存在负值
    return max_alginment_score

def match_percents(seq1, seq2):
    '''
    Return the match_percents of two equal length seq
    '''
    len_seq = len(seq1)
    if len_seq != len(seq2):
        print("two seq length are not equal")
        return
    match = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            match+=1
        match_percent = match/len_seq
    return match_percent

def deltaG_score(seq1, seq2): 
    '''
    calculate the deltaG socre of seq1 and seq2 edited: 20170420
    for incorprate into primer_score , change deltaG_score to positive vaule
    '''
    max_deltaG_score = DG.calDeltaG(seq1, seq2, mono, diva, dntp)

    return -1*max_deltaG_score


def primer_score(primer_1, primer_2, match_method, *extra_parameter):
    '''
    Use different match method, the function match_score, alignment_score, and match_percents deltaG_score to check two primers
    '''
    len_1 = len(primer_1)
    len_2 = len(primer_2)
    if len_2 > len_1: # make primer1 to be the longer one
        primer_1, primer_2 = primer_2, primer_1     
    match_len = min(len_1, len_2)
    len_diff = abs(len_1 - len_2) 
    final_score = -np.inf
     
    if match_method == match_score:
        max_mismatch = extra_parameter[0]
        min_match_len = extra_parameter[1]
        # add min len for match_score
    else:
        #add min len for deltaG_score and alignment_score, match_precent
        min_match_len = extra_parameter[0]
    match_range = range(min_match_len, match_len)

    if match_method == match_score:
        for i in match_range:# the partial match score of two primer
            #文献中为3‘端和any frame, 3'端配对情况
            partial_score_3 = match_method(primer_1[-i:], primer_2[:i], max_mismatch)
            
            if partial_score_3 > final_score:
                final_score = partial_score_3
                flag = 'partial'
                parameter = i #parameter used for the function primer_alignment_print 
            #文献中为3‘端和any frame, 5’端配对情况
        
            partial_score_5 = match_method(primer_1[:i], primer_2[-i:], max_mismatch)
            if partial_score_5 > final_score:
                final_score = partial_score_5
                flag = 'partial'
                parameter = -i 
        for i in range(len_diff+1): # the full match score of two primer
            full_len_score = match_method(primer_1[i:i+match_len], primer_2[:],max_mismatch)
            if full_len_score > final_score:
                final_score = full_len_score
                flag = 'full_len'
                parameter = i
    else:
        for i in match_range:# the partial match score of two primer
            #文献中为3‘端和any frame, 3'端配对情况
            partial_score_3 = match_method(primer_1[-i:], primer_2[:i])
            
            if partial_score_3 > final_score:
                final_score = partial_score_3
                flag = 'partial'
                parameter = i #parameter used for the function primer_alignment_print 
            #文献中为3‘端和any frame, 5’端配对情况
        
            partial_score_5 = match_method(primer_1[:i], primer_2[-i:])
            if partial_score_5 > final_score:
                final_score = partial_score_5
                flag = 'partial'
                parameter = -i
             
        for i in range(len_diff+1): # the full match score of two primer
            full_len_score = match_method(primer_1[i:i+match_len], primer_2[:])
            if full_len_score > final_score:
                final_score = full_len_score
                flag = 'full_len'
                parameter = i
    
    return final_score, flag, parameter



# ### Rule No. 6
# 
# 12 or < 12 ATGC between the 3’-end and any amplicon sequence
# 
# I think this is simple, just return the 3’-end of a primer and check it with the whole amplicon seqs

# In[33]:

# ## Other concerns about the primer design
# I can write my own function to filter rule 2,3,4,5 and I can also use repeatmasker or other tools to check it, I will find some database to check rule no.1
# 
# 1. No SNPs in the primer regions ( within 50 ATGC on any side of target sequences).
# 2. No ≧ 10 consecutive mononucleotides the primer regions.
# 3. No ≧ 9 dinucleotide repeating units the primer regions.
# 4. No ≧ 5 trinucleotide repeating units in the primer regions.
# 5. CG content: 25%–75%.

def simple_repeats_check(primer_dict):
    '''
    Check if the degenerate primers got simple repeats
    '''
    mono_10 = set()
    for base in ATGC:
        mono_10.add(base*10)

    #print(mono_10)

    tri_words = set()
    ATGC2 = ATGC*2
    tri_combinations = itertools.permutations(ATGC2, 3)
    for i in tri_combinations:
        tri_word = ''.join(i)*5
        tri_words.add(tri_word)
    #print (tri_words)

    di_words = set()
    di_combinations = itertools.permutations(ATGC, 2)
    for i in di_combinations:
        di_word = ''.join(i)*9
        di_words.add(di_word)
    #print(di_words)

    repeat_elements = mono_10|di_words|tri_words
    print("\n")
    print('Primers that have simple repeats are :')
    simple_repeats_primers = set()
    for primer_name, primer_seq in primer_dict.items():
        flag = True
        atgc_primers = seq_conversion(primer_seq)
        for primer in atgc_primers:
            if any (x in primer for x in repeat_elements):
                print(">%s" %primer_name)
                print(primer_seq)
                flag = False
                simple_repeats_primers.add(primer_name)
                break

    print ("\n")
    return simple_repeats_primers


# So the basic steps of the function is complete, the next step is put them together to filter the unquanlified primer

# In[37]:

'''
Put the rules together, find the primers that do not meet the needs, print 
which rule they have violated
If the primer violate one rule the return False.
''' 
def primer_check(primer_1, primer_2):
    msg = '{} and {} violate rule No.{:d} the score is {}'
    msg_float = '{} and {} violate rule No.{:d} the score is {:.2f}'
    primer_1 = str.upper(primer_1)
    primer_2 = str.upper(primer_2)
    if len(primer_1) < len(primer_2):
        primer_1, primer_2 = primer_2, primer_1 #保证primer_1是长的那一段  
    rev_cpl_primer_2 = rev_cpl(primer_2) #把primer2转变为反向互补序列进行比较
    panelty_score = 0
    
    #match_score, alignment_score, and match_percents

    rule_1 = three_ends_match_score(primer_1, rev_cpl_primer_2, rule_1_len)
    if rule_1 == rule_1_cutoff:
        print (msg.format(primer_1, primer_2, 1, rule_1))
        primer_alignment_print(primer_1, primer_2, 'partial', rule_1_len)
        panelty_score += 1

    rule_2 = three_ends_match_score(primer_1, rev_cpl_primer_2, rule_2_len)
    if rule_2 >=rule_2_cutoff  :
        print (msg.format(primer_1, primer_2, 2, rule_2))
        primer_alignment_print(primer_1, primer_2, 'partial', rule_2_len)
        panelty_score += 1

    rule_3, flag, parameter = primer_score(primer_1, rev_cpl_primer_2, match_score,0,rule_3_len) #0 means no gap, 7 is the min len
    if rule_3 >= rule_3_cutoff:
        print (msg.format(primer_1, primer_2, 3, rule_3 ))
        primer_alignment_print(primer_1, primer_2, flag, parameter)
        panelty_score += 1

    rule_4, flag, parameter = primer_score(primer_1, rev_cpl_primer_2, match_score,1,rule_4_len)#1 means 1 gap, 7 is the min len
    if rule_4 >= rule_4_cutoff:
        print (msg.format(primer_1, primer_2, 4, rule_4))
        primer_alignment_print(primer_1, primer_2, flag, parameter)
        panelty_score += 1

    rule_5, flag, parameter = primer_score(primer_1, rev_cpl_primer_2,match_percents, rule_5_len)
    if rule_5 >= rule_5_cutoff:
        print (msg_float.format(primer_1, primer_2, 5, rule_5))
        primer_alignment_print(primer_1, primer_2, flag, parameter)
        panelty_score += 1
    
    #Rule No.6 alignment score can not greater than 8
    rule_6, flag, parameter = primer_score(primer_1, rev_cpl_primer_2, alignment_score,rule_6_len) #8 is the min len
    if rule_6 >= rule_6_cutoff:
        print (msg.format(primer_1, primer_2, 6, rule_6))
        primer_alignment_print(primer_1, primer_2, flag, parameter)
        panelty_score += 1
    
    #Rule No.7 deltaG score can not smaller than -9, the min match len is 5. edited: 20170420
    rule_7, flag, parameter = primer_score(primer_1, rev_cpl_primer_2, deltaG_score, rule_7_len)
    #the script in TmDeltaG.py will change to complement sequence, 7 is the min len
    if rule_7 >= rule_7_cutoff:
        rule_7 = -1*rule_7 #change the score back to negative value.
        print (msg_float.format(primer_1, primer_2, 7, rule_7))
        primer_alignment_print(primer_1, primer_2, flag, parameter)
        panelty_score += 1

    
    if panelty_score >= 1:
        print ('The panelty_score is ', panelty_score)
        return False
    else:
        return True

# **for primer_comb in itertools.combinations(primer_dict_atgc.values(), 2):primer_check (primer_comb[0]], primer_comb[1])**   
def primer_amplicon_match_score_max(primer_name, ref_name, primer, ampliseq, strand): 

    '''
    Find the postion of primer in the ampliseq, the largest score, return value is a Primer_Ref_Match
    '''
    #确保primer长度大于length
    # primer_amplicon_match_score = 0
    primer_ref_match_list = list()
    length = len(primer)
    max_score = 0
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
        if match_score > max_score:
            max_score = match_score
            # fhand_primer_positon.write("headline.format(match_score_sub_seq[1],strand, match_score_sub_seq[2]),match_score_sub_seq[0], i, i+length)
            primer_pos_info = Primer_Ref_Match(primer_name, ref_name, primer_seq,ref_seq, max_score, strand, i, i+length)
            # primer_amplicon_match_score +=1

    return primer_pos_info


def primer_position_check(primer_dict, ref_dict, threshold, outputfile_primer_position):
    #Find the postion information of primer, if can not align primer to ref, return 0, else return the alignment object.
    #If the alignment score larger than threshold, write to a file
    fhand_primer_positon = open(outputfile_primer_position, 'w')
    primer_pos_headline = '{:35}\t{:35}\t{:35}\t{:35}\t{:^5}\t{:^5}\t{:^10}\t{:^10}\t{:^10}\n'
    fhand_primer_positon.write(primer_pos_headline.format('Primer_Name', 'Ref_Name' ,'Primer_Seq','Ref_Seq','Score','Strand','Pos_Start','Pos_End','Identical'))
    
    position_dict = dict()
    
    for primer_name, primer_seq in primer_dict.items():
        rev_cpl_primer = rev_cpl(primer_seq)
        primer_max_score = threshold
        for ref_name, ref_seq in ref_dict.items():
            #match object f, if larger than threhold, write to file
            match_object_f = primer_amplicon_match_score_max(primer_name, ref_name, primer_seq, ref_seq, 'forward')
            if match_object_f.match_score > threshold:
                if match_object_f.primer_seq == match_object_f.ref_seq:
                        identical = True
                else:
                    identical = False
                fhand_primer_positon.write(primer_pos_headline.format(match_object_f.primer_name, match_object_f.ref_name, match_object_f.primer_seq, match_object_f.ref_seq, match_object_f.match_score, match_object_f.strand, match_object_f.pos_start, match_object_f.pos_end, identical))

            #match object r, if larger than threhold, write to file
            match_object_r = primer_amplicon_match_score_max(primer_name, ref_name, rev_cpl_primer, ref_seq, 'reverse')
            if match_object_r.match_score > threshold:
                if match_object_r.primer_seq == match_object_r.ref_seq:
                        identical = True
                else:
                    identical = False
                fhand_primer_positon.write(primer_pos_headline.format(match_object_r.primer_name, match_object_r.ref_name, match_object_r.primer_seq, match_object_r.ref_seq, match_object_r.match_score, match_object_r.strand, match_object_r.pos_start, match_object_r.pos_end, identical))

            #if match score f or r larger than pervious max score, it is the new max score, store the max score in primer_max_score
            if match_object_f.match_score > primer_max_score:
                primer_max_score = match_object_f.match_score
                max_match_object = match_object_f
            if match_object_r.match_score > primer_max_score:
                primer_max_score = match_object_r.match_score
                max_match_object = match_object_r

        #if primer_max score larger than threshold,store it in position_dict 
        if primer_max_score > threshold:
            position_dict[primer_name] = max_match_object
        else:
            position_dict[primer_name] = 0
            fhand_primer_positon.write('{}\t{}\n'.format(primer_name, 0))
                
    fhand_primer_positon.close()
    return position_dict

def remove_primer_pos_wrong(primer_dict, ref_dict, position_dict):
    #Eliminate the primer if position_dict value is 0
    primer_dict_copy = copy.deepcopy(primer_dict)
    for primer_name in primer_dict_copy:
        if position_dict[primer_name] == 0:
            primer_dict.pop(primer_name, None)
    return primer_dict


def remove_primer_info_wrong(primer_dict, outputfile_primer_information, outputfile_primer_atgc):
    '''
    remove the primer degenerate number larger than threshold, or tm out of range. Write the primer information to file outputfile_primer_information. Write the primer atgc file to outputfile_primer_atgc
    '''
    group_name_dict = dict() #A dict contain all the group names
    primer_index_dict = dict() #A dict contain the index of primer_name in group_name_dict
    fhand_primer_atgc = open (outputfile_primer_atgc,'w')
    fhand_primer_information = open(outputfile_primer_information,'w')
    headline = '{:35}\t{:40}\t{:^10}\t{:^10}\t{:^5}\t{:^5}\t{:^10}\n'
    fhand_primer_information.write(headline.format('Name','Sequence','Len','DegNumber','MaxTm','MinTm','Flag'))
    #Remove the primer information not quanlified
    primer_dict_copy = copy.deepcopy(primer_dict)
    for primer_name, primer_seq in primer_dict_copy.items():
        if primer_information_check(primer_name,primer_seq,fhand_primer_information) == True:
        #check if the primer tm in range, and degenerate number <= threshold value
            # primer_name_set_degenerate.add(primer_name)
            group_name = '#'.join(primer_name.split('#')[:-1])
            # primer_index = primer_name.split('#')[-1] 
            #if primer_name is >Norovirus_GI#R#3, primer_index is 3
            # modified in 20170802, using number to represent the primer_index, so the minimize the limitation to primer name.
            if group_name not in group_name_dict:
                group_name_dict[group_name] = dict()
                primer_index_dict[group_name] = 0
                primer_index = primer_index_dict[group_name]
                group_name_dict[group_name][primer_index] = primer_name
            else:
                primer_index_dict[group_name] += 1
                primer_index = primer_index_dict[group_name]
                group_name_dict[group_name][primer_index] = primer_name
            primer_dict[primer_name] = primer_seq

            primers = seq_conversion(primer_seq)
            for index, primer in enumerate(primers):
                primer_name_index = '>'+primer_name+'#'+str(index)
                fhand_primer_atgc.write(primer_name_index+'\n')
                fhand_primer_atgc.write(primer+'\n')
        else:
            primer_dict.pop(primer_name, None)
    fhand_primer_atgc.close()
    fhand_primer_information.close()
    return primer_dict, group_name_dict


def read_primer_atgc(outputfile_primer_atgc):
    primer_dict_atgc ={} # A dictionary contains all the atgc primer sequences
    with open(outputfile_primer_atgc) as fhand:
        for line in fhand:
            if line.startswith('>'): #e.g: >adenovirus#2#F#0
                primer_name = line[1:].strip() # adenovirus#2#F#0
            else:
                primer_seq = line.strip()
                primer_dict_atgc[primer_name] = primer_seq
    return primer_dict_atgc

def primers_confliction_report(outputfile_confliction_report,outputfile_confliction_csv, primer_dict, group_name_dict, confliction_table):
    '''
    Write the conflication report to file
    '''
    sys.stdout=open(outputfile_confliction_report,"w")
    simple_repeats_primer_set = simple_repeats_check(primer_dict)
    for primer_name in simple_repeats_primer_set:
        confliction_table.loc[primer_name,primer_name] += 1 #set the confliction table of simple_repeats_primers to 1

    print ("The confliction of primers are : \n")
    msg_6 = '{} and {} may have interaction.\n'
    degenerate_unusable_primers_dict = dict()
    for primer_comb in itertools.combinations_with_replacement(primer_dict.keys(), 2): #use the primer_name of each primer, check if the two primers got an interaction
    #delete the primers have simple repeats
        degenerate_primers_0 = DegeneratePrimer(primer_comb[0],primer_dict[primer_comb[0]]).atcg_sequence()
        degenerate_primers_1 = DegeneratePrimer(primer_comb[1],primer_dict[primer_comb[1]]).atcg_sequence()
        for primer_atcg_comb in itertools.product(degenerate_primers_0,degenerate_primers_1):
            if primer_check(primer_atcg_comb[0],primer_atcg_comb[1]) == False:
                degenerate_unusable_primers_dict[primer_comb[0]]= degenerate_unusable_primers_dict.get(primer_comb[0],0)+1
                degenerate_unusable_primers_dict[primer_comb[1]]= degenerate_unusable_primers_dict.get(primer_comb[1],0)+1
                print(msg_6.format(primer_comb[0],primer_comb[1]))
                #Convert the primer name to degenerate primer name
                #If there is a confliction between primers, write 1 to the corresponding position in confliction table
                confliction_table.loc[primer_comb[0],primer_comb[1]] += 1
                confliction_table.loc[primer_comb[1],primer_comb[0]] += 1
                break         
    print ('\n')
    print ("The panelty score of degenerate primers are :")
    msg_5 = '{:35}\t{:^15}'
    for name,value in degenerate_unusable_primers_dict.items():
        print(msg_5.format(name, value))

    with open(outputfile_confliction_csv,'w') as fhand:
            confliction_table.to_csv(fhand)

        #Eliminate the primers with self interaction
    self_interaction_primer = set()
    if elminate_self_interaction :
        group_name_dict_copy = copy.deepcopy(group_name_dict)
        for key in group_name_dict_copy.keys():
            for primer_index,primer_name in group_name_dict_copy[key].items():
                if confliction_table.loc[primer_name,primer_name]>=1:
                    self_interaction_primer.add(primer_name)
                    group_name_dict[key].pop(primer_index,None)

    
    print("==================================================")
    print("The group name dict is: ")
    print_dict(group_name_dict)

    print('The self interaction primers are:')
    for primer_name in self_interaction_primer:
        print('>'+ primer_name)
        print(primer_dict[primer_name])
    sys.stdout.close()

def probe_region_generation(primer_list, ref_dict, position_dict):
    #strip the primer sequence and primer probe gap sequence from ref sequence, then it is the sequence for probe design. 
    probe_region_dict = copy.deepcopy(ref_dict)

    for primer_name in primer_list:
        #the key of position_dict  primer_name, the value of position_dict is Primer_Ref_Match.
        position_info = position_dict[primer_name]
        if position_info.strand == 'forward':
            probe_start_index = len(position_info.primer_seq) + primer_probe_gap
            #Change the probe region, trime the primer_seq + primer_probe_gap from ref_seq
            probe_region_dict[position_info.ref_name] = probe_region_dict[position_info.ref_name][probe_start_index:]
        elif position_info.strand == 'reverse':
            probe_end_index = len(position_info.primer_seq) + primer_probe_gap
            probe_region_dict[position_info.ref_name] = probe_region_dict[position_info.ref_name][:-1*probe_end_index]
    return probe_region_dict

def probe_candi_dict_generation(probe_region_dict, probe_info_file_name):
    '''
    Check the probe in probe_region_seq, which are store in probe_region_dict, if it meets the criteria in probe_check, store in the probe_candi_dict, the key is ref_name
    '''
    probe_candi_dict = defaultdict(list)
    probe_info = '{:35}\t{:35}\t{:^4}\t{:^4}\t{:^.2f}\t{:^.2f}\t{:^4}\t{:^.2f}\t{:^.2%}\t{:^4}\t{:^.2f}\t{:^.2%}\n'
    probe_info_header = '{:35}\t{:35}\t{:^4}\t{:^4}\t{:^4}\t{:^4}\t{:^4}\t{:^4}\t{:^4}\t{:^4}\t{:^4}\t{:^4}\n'
    probe_info_fhand = open(probe_info_file_name, 'w')
    probe_info_fhand.write(probe_info_header.format('Probe_Name', 'Probe_Seq', 'Ambigous_Number', 'Deg_Number', 'Tm_min', 'Tm_max', 'Hairpin_Number', 'Max_Hairpin_Tm', 'Hairpin_Ratio','Homodimer_number', 'Min_HomoDimer_dg','HomoDimer_Ratio'))
    for ref_name, ref_seq in probe_region_dict.items():
        for probe_len in range(probe_len_min, probe_len_max+1):
            for i in range(len(ref_seq) - probe_len):
                probe_name = ref_name + '@' + str(i)
                probe_seq = ref_seq[i:i+probe_len]
                if probe_check(probe_name, probe_seq, probe_info_fhand, probe_info):
                    probe_candi_dict[ref_name].append(probe_seq)

    probe_info_fhand.close()
    return probe_candi_dict


def probe_check(probe_name, probe_seq, fhand, probe_info):
    '''
    Check if the probe_seq meets certain criteria
    '''
    if probe_init_noG_check: #this parameter is true, check if the first letter of probe is G
        if probe_init_G(probe_seq): #probe initiate with G, return false
            return False
    ambigous_number = seq_ambigous_number(probe_seq)
    if  ambigous_number > probe_ambigous_letter_threshold:
        # check if the ambigous number of the probe is larger than the threshold
        return False
    degenerate_probe = DegeneratePrimer(probe_name, probe_seq)
    degenerate_number = degenerate_probe.degenerate_number()
    if  degenerate_number > probe_max_deg_threshold:
        # check if the degenerate number of the probe is larger than the threshold
        return False
    tm_max, tm_min = degenerate_probe.thermol_result()
    if tm_max > probe_tm_max or tm_min < probe_tm_min:
        return False

    homodimer_number, min_homo_dimer_dg = degenerate_probe.homo_dimer_check()
    homodimer_ratio = homodimer_number/degenerate_number
    if homodimer_ratio > homodimer_number_cutoff : #Degenerate Probe is unqualified only if larger than homodimer_number_cutoff (default set ot 25%) of the atgc probes got homo_dimer (dg smaller than homodimer_dg_threshold)
        return False

    hairpin_number, max_hairpin_tm = degenerate_probe.hairpin_check()
    hairpin_ratio = hairpin_number/degenerate_number
    if hairpin_ratio > hairpin_number_cutoff : #Degenerate Probe is unqualified only if larger than hairpin_number_cutoff (default set ot 25%) of the atgc probes got hairpin structure (tm larger than hairpin_tm_threshold)
        return False

    else:
        fhand.write(probe_info.format(probe_name, probe_seq, ambigous_number, degenerate_number, tm_min, tm_max, hairpin_number, max_hairpin_tm, hairpin_ratio, homodimer_number, min_homo_dimer_dg, homodimer_ratio))
        return True

def primer_combinations(primer_combinations_file,outputfile_probe_info, group_name_dict, confliction_table, primer_amplicon_match_dict, position_dict, primer_dict, ref_dict):
    '''
    Check the primer combinations according to the DataFrame file , and Write the Primer combinations to file
    '''
    
    #file that store the probe information
    combination_number = 0
    primer_without_interaction = dict()
    combination_without_ineraction_number = 0
    with open(primer_combinations_file,'w') as fhand:
        fhand.write('Combinations of primers are:\n')
        for product in itertools.product(*group_name_dict.values()):
            #print(product)
            primer_list = []
            flag = False #no interation
            combination_number +=1
            for index, group_name in enumerate(group_name_dict.keys()):
                primer_list.append(group_name_dict[group_name][product[index]])
            primer_tuple = tuple(primer_list)
            amplicon_len_flag, amplicon_len_list = amplicon_length_check(amplicon_len_min,amplicon_len_max,amplicon_lendiff_min,primer_tuple, position_dict, primer_amplicon_match_dict) 
            if amplicon_len_flag == False:
                flag = True #Amplicon length is inappropriate
            if flag == False:
                for comb in itertools.combinations(primer_list,2):
                    # if confliction_table.loc[comb[0],comb[1]] ==confliction_table_threshold or confliction_table.loc[comb[1],comb[0]] ==confliction_table_threshold: 
                    #include the self interaction primers
                    if confliction_table.loc[comb[0],comb[1]] >=confliction_table_threshold or confliction_table.loc[comb[1],comb[0]] >=confliction_table_threshold:
                        flag = True #different primer interaton spotted
                        break
            if flag == False: #no interaction, write the primers to the file
                combination_without_ineraction_number +=1
                fhand.write('Primer Combinations Without Interaction:\n')
                fhand.write('Amplicon Length are : %s\n'%str(amplicon_len_list) )
                #Add the function to design probe according to the given primer_lsit
                probe_region_dict = probe_region_generation(primer_list, ref_dict, position_dict)
                probe_candi_dict = probe_candi_dict_generation(probe_region_dict, outputfile_probe_info)
                for primer_name in primer_list:
                    primer_without_interaction[primer_name] = primer_without_interaction.get(primer_name,0)+1
                    fhand.write('>'+primer_name+'\n')
                    fhand.write(primer_dict[primer_name]+'\n')
                fhand.write('\n')
                fhand.write("The probe can be used for the primers are: \n")
                for ref_name, probe_seqs in probe_candi_dict.items():
                    fhand.write('>'+ref_name+'\n')
                    for probe_seq in probe_seqs:
                        fhand.write(probe_seq + '\n')
                fhand.write('\n')
        
        fhand.write('\nThe total combination number is %d, the combination number without interaction is %d\n'%(combination_number, combination_without_ineraction_number))
        fhand.write('The primers can be used are:\n')
        for key in primer_without_interaction.keys():
            fhand.write('>%s\n'%key)
            fhand.write('%s\n'%primer_dict[key])
        fhand.write('\n')
        fhand.write('The primers can be used and there appearance value are:\n')
        for key, value in primer_without_interaction.items():
            fhand.write('{:35}\t{:d}\n'.format(key,value))

        #Add the report of primer amplicon interaction
        fhand.write('\nPrimer_Amplcion interaction number are:\n')
        msg_primer_amplicon = '{:30}\t{:d}\n'
        for key, value in primer_amplicon_match_dict.items():
            fhand.write(msg_primer_amplicon.format(key,value))
        fhand.write('\nPrimers with unspecific interactions are:\n')
        for key, value in primer_amplicon_match_dict.items():
            if value>1:
                fhand.write('>%s\n'%key)
                fhand.write('%s\n'%primer_dict[key])


def main():
    primer_file = sys.argv[1]
    ref_file = sys.argv[2]

    primer_dict = read_fasta_file(primer_file)  
    ref_dict = read_fasta_file(ref_file)

    outputfile_primer_position = primer_file.split(".")[0] + "_primer_position.txt"
    position_dict = primer_position_check(primer_dict, ref_dict, threshold_position, outputfile_primer_position)
    primer_dict = remove_primer_pos_wrong(primer_dict, ref_dict, position_dict)
 

    outputfile_primer_atgc = primer_file.split(".")[0] + "_atgc.fasta"
    outputfile_primer_information = primer_file.split(".")[0] + "_primer_information.txt"
    primer_dict, group_name_dict = remove_primer_info_wrong(primer_dict, outputfile_primer_information, outputfile_primer_atgc)

    primer_dict_atgc = read_primer_atgc(outputfile_primer_atgc)
    
    #the output file of primer amplicon interaction report
    primer_amplicon_outputfile = primer_file.split(".")[0]+'_'+ ref_file.split(".")[0]+"_confliction.fasta"
    confliction_table_len = len(primer_dict.keys())
    confliction_table = pd.DataFrame(np.zeros((confliction_table_len, confliction_table_len)), index = primer_dict.keys(), columns = primer_dict.keys())


    #Find the position of primers and record in position_dict


    #Add the function to eliminat the primer pairs if both f and r had unspecific priming
    primer_amplicon_match_dict = primer_amplicon_check(primer_dict, ref_dict, amplicon_primer_check_len, amplicon_primer_cutoff)

    #Write the primers confliction report 
    outputfile_confliction_report = primer_file.split(".")[0]+"_confliction_report.fasta"
    outputfile_confliction_csv = primer_file.split(".")[0]+'_confliction_table.csv'
    primers_confliction_report(outputfile_confliction_report, outputfile_confliction_csv, primer_dict, group_name_dict,confliction_table)
    
    #Write the primer combinations
    outputfile_probe_info = ref_file.split(".")[0] + "_probes_info.txt"
    primer_combination_file = primer_file.split(".")[0]+'_primer_probe_combinations.txt'
    primer_combinations(primer_combination_file, outputfile_probe_info,group_name_dict, confliction_table, primer_amplicon_match_dict, position_dict, primer_dict, ref_dict)

    
main()

