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
#modification : change to the whole length of template, not only the amplicons, use the consensus sequence to update amplicon sequence.
# In[186]:

import sys
ATGC = {'A','T','G','C'}
IUPAC_DICT = {'R':{'A','G'}, 'Y':{'C','T'},'S':{'G','C'},'W':{'A','T'},
              'K':{'G','T'},'M':{'A','C'},'B':{'C','G','T'},'D':{'A','G','T'},
              'H':{"A","C","T"},'V':{"A","C","G"},'N':{"A","T","C","G"},
             "X":{"A","T","G","C"}, 'A':{'A'}, 'G':{'G'}, 'T':{'T'}, 'C':{'C'}}


# In[172]:

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


# In[27]:

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


# In[152]:

#read the primers in memory
primer_file = sys.argv[1]
#primer_file = input("Enter the file name of primers :")
if len(primer_file) == 0:
    primer_file = 'rotavirus_primers.fasta'
primer_dict = read_fasta_file(primer_file) 

#read the consensus sequences in memory
consensus_seq_file = sys.argv[2]
#consensus_seq_file = input("Enter the file name of consensus sequences :")
if len(consensus_seq_file) == 0:
    consensus_seq_file = 'rotavirus_a_consensus_seq.fasta'
consensus_seq_dict = read_fasta_file(consensus_seq_file) 


# In[ ]:

output_file_name = primer_file.split(".")[0]+'_'+consensus_seq_file.split(".")[0]+"_confliction.fasta"
sys.stdout=open(output_file_name,"w")




print("===============================================")
print ("The consensus seq dict is:")
print("===============================================")
print_dict(consensus_seq_dict)

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
    #primer_3_end = primer[-length:]
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
            print("Found one match between the primer and amplicon above threshold, the match score is {}.\nThe primer 3 end seq is {} it is on {} strand.\nThe amplicon seq is {}".                  format(match_score_sub_seq[0],match_score_sub_seq[1],strand, match_score_sub_seq[2]))
            primer_amplicon_match_score +=1
    return primer_amplicon_match_score

def primer_amplicon_check(primer_dict, consensus_seq_dict, length, threshold):
    '''
    Check the interaction between the primers and amplicons, \
    the primers and amplicons are stored in different dicionary.
    '''
    primer_amplicon_match_dict = dict()
    for amplicon_name, amplicon_seq in consensus_seq_dict.items():
        for primer_name, primer_seq in primer_dict.items():
            rev_cpl_primer = rev_cpl(primer_seq)
            match_score_f = primer_amplicon_match_score(primer_seq, amplicon_seq, length,threshold, 'forward') 
            #正向引物与amplcion的match score
            match_score_r = primer_amplicon_match_score(rev_cpl_primer, amplicon_seq, length,threshold, 'reverse')
            #反向引物与amplicon的match score,反向引物应该为引物前半段与序列的比较
            match_score = match_score_f + match_score_r
            if match_score >=1:
                print("Find {} matches of primer {} in amplicon {}".format(match_score,primer_name, amplicon_name))
                print('\n')
                primer_amplicon_match_dict[primer_name] = primer_amplicon_match_dict.get(primer_name,0)+match_score
    
    for primer, value in primer_amplicon_match_dict.items():
        if value>1:
            print(primer,value)

    return primer_amplicon_match_dict


# In[184]:

print("===============================================")
print("The match score between the primer and amplicons are :")
print("===============================================")
primer_amplicon_check(primer_dict, consensus_seq_dict, 14, 12)


# In[ ]:

sys.stdout.close()

