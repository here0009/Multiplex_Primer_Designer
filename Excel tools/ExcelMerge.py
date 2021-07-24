# coding: utf-8
# Merge Primer Files, if the primerName and primerSeq are identical, deduplicate. If the primerName is the same while the primerSeq is different, change the primerName of the second primer.
import os
import sys
import copy

def read_fasta_file(file_name, fasta_dict, uniquePrimer) :
    '''
    Function to read fasta file, store the names and sequences in a dictionary
    '''
     # a dict store the primers
    number = 1
    with open (file_name) as fhand :
        for line in fhand:
            if line.startswith('>'):
                name = line.strip()[1:]
                seq = fhand.readline().strip()
                # print(name, seq)
                if seq in uniquePrimer :
                    print("Delete primer name of {}".format(name))
                    print("Already exsiting primer sequence of {}".format(seq))
                else:
                    uniquePrimer.add(seq)
                    if name in fasta_dict:
                        print("Already exsiting fasta name of {}".format(name))
                        if seq == fasta_dict[name]:
                            print("Duplicated primer of {}".format(name))
                        else:
                            name = name + str(number)
                            number+=1
                            print("Change it to the name of {}".format(name))
                    fasta_dict[name] = seq
    # print (uniquePrimer)
    return fasta_dict

def load_files(directory_name, fasta_dict, uniquePrimer):
    '''
    read the files in directory_name, store the sequences information in fasta_dict
    '''
    input_file_list = list()
    for file_name in os.listdir(directory_name):
        if file_name.endswith(".fasta"):
            input_file_list.append(file_name)
            print(file_name) #for test

    print("Read the following files into dict:")
    print("===================================")
    for input_file in input_file_list:
        input_file = '\\'.join([directory_name, input_file])
        print(input_file) 
        with open(input_file) as fhand:
                read_fasta_file(input_file, fasta_dict, uniquePrimer)

def write_dict(test_dict, file_name):
    '''
    Function to print a dictionary, so it is easy to read
    '''
    with open(file_name, 'w') as fhand:
        for key, value in test_dict.items():
            fhand.write('>{}\n'.format(key))
            fhand.write('{}\n'.format(value))


def main():
    directory_name = sys.argv[1].replace('\\','\\')
    fasta_dict = dict()
    uniquePrimer = set()
    load_files(directory_name, fasta_dict, uniquePrimer)
    outputfile = "Merged_sequences.txt"
    write_dict(fasta_dict, outputfile)
    
main()