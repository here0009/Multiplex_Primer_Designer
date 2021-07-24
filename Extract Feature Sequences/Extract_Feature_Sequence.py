
# coding: utf-8

# In[1]:

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
import os
import sys
#Parameters that needed tob be modified

search_terms = {"tcdB","toxin b","toxin B"} #The CDS name you want to extract
min_seq, max_seq = (7000,7200) #The CDS length specified
# search_terms = {'protein utxA', 'protein UtxA'}
# min_seq, max_seq = (400,600) #The CDS length specified
output_file = "output_file.fasta" #the output file name
outrange_file = "outrange.fasta" #stores the sequence records out of range

# ### Note
# 
# The main purpose of this script is to extract feature sequences of given genbank files.
# 
# **input:** 
# 
# Genbank files:the file name are stored in input_file_list(commented out), or it can read all the files in the current directory(using now). 
# 
# The feature you want to extrat : search_terms
# can change  the parameter of  find_feature_seq(record, "CDS", "product",search_terms) to specify the features you want to extrat
#     
# **process:**
# 
# Extract feature and written them to an output file. Mainly using BioPython package to do the job.
# 
# **output:**
# 
# fasta_headline will inform the record name|feature location|feature length.

fasta_headline_format = ">{}|{}|{}|{:d}\n"
record_number = 0
find_item = 0
written_item = 0
outrange_item = 0
outputfile_index = 0
output_fhand = open(output_file,'w')
outrange_fhand = open(outrange_file,'w')

def find_feature_seq(seqrecord, feature_type, qualifier_key, search_terms):
    '''
    Find the feature seq in a given genbank record stored in a SeqRecord object
    The parameters are :
    record: the SeqRecord object
    feauter_type: CDS or gene (feature.type)
    qualifier_keys: the key of (feature.qualifiers["product"])
    search_term : the one you want to search, store in a set ({"tcdB","toxin b","toxin B"})
    '''
    flag = False
    for feature in seqrecord.features:
        if feature.type == feature_type:
            if qualifier_key in feature.qualifiers:
                qualifier_name = feature.qualifiers[qualifier_key][0]
                #print(feature.qualifiers)
                if qualifier_name in search_terms:
#                     print(seqrecord.name)
#                     print(product_name)
#                     print(feature.location)
                    feature_seq = feature.location.extract(seqrecord).seq
                    #output_seq = Seq(feature_seq)
                    feature_seq_record = SeqRecord(feature_seq)
                    feature_seq_record.name = seqrecord.name
                    feature_seq_record.id = seqrecord.id
                    feature_seq_record.description = qualifier_name
                    feature_seq_record.features.append(feature)
#                     print(len(seq))
                    flag = True
                    break #find one and stop, in case find multiple matches
    if flag == False:
        feature_seq_record = None
    return flag, feature_seq_record


# In[37]:

def write_record_to_file(record, fasta_headline_format, fhand):
    '''
    Write the record to fhand
    '''
    feature_location = record.features[0].location
    fasta_headline = fasta_headline_format.format(record.id, record.description,                                                   feature_location, len(record.seq))
    fhand.write(fasta_headline)
    fhand.write(str(record.seq)+'\n')


# In[38]:

directory_name = sys.argv[1].replace('\\','\\')
input_file_list = list()
for file_name in os.listdir(directory_name):
    if file_name.endswith(".gb"):
        input_file_list.append(file_name)

# In[47]:

#Specify the parameters



# In[48]:

# write the output to file
for input_file in input_file_list:
    input_file = '\\'.join([directory_name, input_file])
    with open(input_file) as fhand:
            records = SeqIO.parse(fhand,"genbank")
            for record in records:
                record_number += 1
                flag,feature_seq_record = find_feature_seq(record, "CDS", "product",search_terms)
                if flag:
                    len_seq = len(feature_seq_record.seq.strip(chars="N")) 
                    #some of the contig or wgs files do not have sequence, they are
                    #make up of Ns, eliminate these sequences
#                     len_seq = len(feature_seq_record.seq)
                    find_item+=1
                    if len_seq in range(min_seq, max_seq):
                    #if there are too many records, split the output file, currently 
                    #inactivated
#                         if written_item%record_per_file == 0:
#                             outputfile_index +=1
#                             output_fhand.close()
#                             file_name, appendix_name = output_file.split(".")
#                             output_file_name = file_name+'_'+str(outputfile_index)+'.'+appendix_name
# #                             output_file_name = change_output_file(output_fhand, output_file,outputfile_index )
#                             output_fhand = open(output_file_name,'w')
#                         else:
                        write_record_to_file(feature_seq_record, fasta_headline_format,output_fhand)
                        written_item += 1
                    else:
                        write_record_to_file(feature_seq_record,fasta_headline_format, outrange_fhand)
                        outrange_item +=1
output_fhand.close()
outrange_fhand.close()
print("There are total %d records in files."%record_number)
print("There are total %d featured sequence generated."%find_item)
print("There are total %d sequences have written to the file."%written_item)
print("There are total %d sequences are out range."%outrange_item)
