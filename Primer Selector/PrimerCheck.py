#!/usr/bin/python
#coding:utf-8
#modification : merge the function of primer check and primer amplicon check. The primers are forward and reverse primers are allocated to different group.
#Disable the function of amplicon check in this version.
# Test of Git
# ### Import the needed module, and read the test primers into memory
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

#Parameters to calculate delatG
mono = 50 #mono-valent ion concentration, mM
diva = 1.5 #divalent ion concentraion, mM
oligo = 50 #oligo concentraion to calculate Tm
dntp = 0.25 #dntp concentraction, mM

#Parameters to calculate primer interactions
rule_1_len = 5 #3' ends len
rule_1_cutoff = 5 #3' ends cutoff score
rule_2_len = 8 #3' ends len(with one mismatch)
rule_2_cutoff = 7 #3' ends cutoff score(with one mismatch)
rule_3_len = 10 #primer len
rule_3_cutoff = 10 #primer cutoff score
rule_4_len = 13 #primer len(with one mismatch)
rule_4_cutoff = 12 #primer cutoff score(with one mismatch)
rule_5_len = 14 #match percent calculation min len
rule_5_cutoff = 0.75 #match percent cutoff score
rule_6_len = 9 #alignment score len
rule_6_cutoff = 9 #alignment cutoff score
rule_7_len = 7 #delta G calculation min len
rule_7_cutoff = -9 #delta G cutoff score, kcal/mol
rule_7_cutoff = -1*rule_7_cutoff #change cut-off score to positive value

#Dictionary of ATGC and Degenerate Bases
ATGC = ['A', 'T', 'G', 'C']
IUPAC_DICT = {'R':{'A','G'}, 'Y':{'C','T'},'S':{'G','C'},'W':{'A','T'},
              'K':{'G','T'},'M':{'A','C'},'B':{'C','G','T'},'D':{'A','G','T'},
              'H':{"A","C","T"},'V':{"A","C","G"},'N':{"A","T","C","G"},
             "X":{"A","T","G","C"}}



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


# In[9]:

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
    return seqs_converted

# In[10]:

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

# inputfile = input("Enter the primer design results : ") #fasta format primer name example: >adenovirus#2#F
inputfile = sys.argv[1]
primer_dict_degenerate = dict()

primer_name_set_degenerate = set() #A set contain all the degenerate primer names
outputfile_2 = inputfile.split(".")[0] + "_atgc.fasta"
fhand_output = open (outputfile_2,'w')
with open(inputfile) as fhand_input:
    for line in fhand_input:
        if line.startswith('>'):
            primer_name = line.strip()[1:] #e.g: adenovirus#F#2
            primer_name_set_degenerate.add(primer_name)
            primer_seq = fhand_input.readline().strip()
            primer_dict_degenerate[primer_name] = primer_seq
            primers = seq_conversion(primer_seq)
            for index, primer in enumerate(primers):
                primer_name_index = '>'+primer_name+'@'+str(index)
                fhand_output.write(primer_name_index+'\n')
                fhand_output.write(primer+'\n')
fhand_output.close()


primer_dict_atgc ={} # A dictionary contains all the atgc primer sequences



with open(outputfile_2) as fhand:
    for line in fhand:
        if line.startswith('>'): #e.g: >adenovirus#2#F#0
            primer_name = line[1:].strip() # adenovirus#2#F#0
        else:
            primer_seq = line.strip()
            primer_dict_atgc[primer_name] = primer_seq


confliction_table_len = len(primer_name_set_degenerate)
confliction_table = pd.DataFrame(np.zeros((confliction_table_len, confliction_table_len)), index = primer_name_set_degenerate, columns = primer_name_set_degenerate)


def amplicon_length_check(amplicon_min,amplicon_max,amplicon_diff,primer_tuple):
    '''
    Return True if the amplicon len is between min and max and the lenght diffrent between amplicons is above threshold.
    '''
    flag = True
    amplicon_len_list = []
    for comb in itertools.combinations(primer_tuple,2):
        if position_dict[comb[0]][0] == position_dict[comb[1]][0]: #the same organism
            if position_dict[comb[0]][1] == 'reverse':
                reverse_position = position_dict[comb[0]][2]
                forward_position = position_dict[comb[1]][2]
            else:
                reverse_position = position_dict[comb[1]][2]
                forward_position = position_dict[comb[0]][2]
            amplicon_len = reverse_position - forward_position
            amplicon_len_list.append(amplicon_len)
            if amplicon_len<amplicon_min or amplicon_len>amplicon_max:
                flag = False
    amplicon_len_list = sorted(amplicon_len_list)
    for i in range(len(amplicon_len_list)-1):
        if amplicon_len_list[i+1] - amplicon_len_list[i] < amplicon_diff: #amplicon_diff is the minimum lenghth diff between different amplicons
            flag = False
    return flag


outputfile_3 = inputfile.split(".")[0]+"_confliction_report.fasta"
sys.stdout=open(outputfile_3,"w")

# ### Filter of Rule No.1 and No.2, they are usable

# '''
# Rule No.1
# If seq1 and seq2 got 4 complementary seq in 3 end, 
# return True,
# else return False
# '''
# def primer_check(primer1, primer2):
#     if primer1[-4:] == primer2[-4:]:
#         return True
#     else:
#         return False

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

'''
Useful
Modify primer_check
given seq1, seq2 and length 3 parameters
return the complemantry score of seq1 and seq2
This can be applied to Rule1, and 2
'''
def three_ends_match_score(primer_1, primer_2, length):
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
    hits : score + 1
    gap open : score - 2
    gap extension : score -1
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

# In[34]:
def simple_repeats_check(primer_dict_degenerate):
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
    for primer_name, primer_seq in primer_dict_degenerate.items():
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
If the primer violate one rule the panelty_score will plus one
Finaly return the panelty_score, if panelty_score >= 1 return false
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

def main():
    # In[27]:
    simple_repeats_set = simple_repeats_check(primer_dict_degenerate)
    for primer_name in simple_repeats_set:
        confliction_table.loc[primer_name,primer_name]+=1
        
    print ("The confliction of primers are : \n")
    msg_6 = '{} and {} may have interaction.\n'
    unusable_primers_dict = dict()
    for primer_comb in itertools.combinations(primer_dict_atgc.keys(), 2): #use the primer_name of each primer, check if the two primers got an interaction
        if primer_check (primer_dict_atgc[primer_comb[0]], primer_dict_atgc[primer_comb[1]]) == False:
            unusable_primers_dict[primer_comb[0]]= unusable_primers_dict.get(primer_comb[0],0)+1
            unusable_primers_dict[primer_comb[1]]= unusable_primers_dict.get(primer_comb[1],0)+1
            print(msg_6.format(primer_comb[0],primer_comb[1]))
            #Convert the primer name to degenerate primer name
            primer_name_degenerate_0 = primer_comb[0].split('@')[0]
            primer_name_degenerate_1 = primer_comb[1].split('@')[0]
            #If there is a confliction between primers, write 1 to the corresponding position in confliction table
            confliction_table.loc[primer_name_degenerate_0,primer_name_degenerate_1] += 1
            confliction_table.loc[primer_name_degenerate_1,primer_name_degenerate_0] += 1
        

    print ('\n')
    print ("The panelty score of primers are :")
    msg_4 = '{:35}\t{:^15}\t{:35}'
    print(msg_4.format("Primer Name", "Penalty Value", "Primer Seq"))
    for name, value in primer_dict_atgc.items():
        if name in unusable_primers_dict.keys():
            print(msg_4.format(name, unusable_primers_dict[name], value))


    degenerate_unusable_primers_dict = dict()
    for name in unusable_primers_dict.keys():
        degenerate_name = name.split('@')[0]
        degenerate_unusable_primers_dict[degenerate_name] = degenerate_unusable_primers_dict.get(degenerate_name,0)+ unusable_primers_dict[name]

    print ('\n')
    print ("The panelty score of degenerate primers are :")
    msg_5 = '{:35}\t{:^15}'
    for name,value in degenerate_unusable_primers_dict.items():
        print(msg_5.format(name, value))
            
    #Eliminate the primers with self interaction
    self_interaction_primer = set()
    for primer_name in degenerate_unusable_primers_dict.keys():
        if confliction_table.loc[primer_name,primer_name]>=1:
            self_interaction_primer.add(primer_name)
    
    print('\n')
    print('The self interaction primers are:')
    for primer_name in self_interaction_primer:
        print('>%s'%primer_name)
        print('%s'%primer_dict_degenerate[primer_name])

    sys.stdout.close()

    with open(inputfile.split(".")[0]+'_confliction_table.csv','w') as fhand:
        confliction_table.to_csv(fhand)

    

main()
