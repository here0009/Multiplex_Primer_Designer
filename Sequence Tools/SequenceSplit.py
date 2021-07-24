# coding: utf-8
# Split a fasta file contains multiple sequences to different fasta files, contain one sequence each.
import sys
def read_fasta_file(file_name, fasta_dict) :
    '''
    Function to read fasta file, store the names and sequences in a dictionary
    '''
     # a dict store the primers
    number = 1
    with open (file_name) as fhand :
        for line in fhand:
            if line.startswith('>'):
                name = line.strip()[1:]
                seq = primer_seq = fhand.readline().strip()
                # print(name, seq)
                if name in fasta_dict:
                    print("Already exsiting fasta name of {}".format(name))
                    if seq == fasta_dict[name]:
                        print("Duplicated primer of {}".format(name))
                    else:
                        name = name + '&' + str(number)
                        number+=1
                        print("Change it to the name of {}".format(name))
                fasta_dict[name] = seq
    return fasta_dict

def write_fasta(fasta_dict):
    '''
    Function to print a dictionary, so it is easy to read
    '''
    n = len(fasta_dict)
    for key, value in fasta_dict.items():
        with open('{}.fasta'.format(key), 'w') as fhand:
            fhand.write('>{}\n'.format(key))
            fhand.write('{}\n'.format(value))

# def print_dict(test_dict):
#     '''
#     Function to print a dictionary, so it is easy to read
#     '''
#     for key, value in test_dict.items():
#         print('{}\t{}'.format(key, value))


def main():
    file_name = sys.argv[1]
    fasta_dict = dict()
    read_fasta_file(file_name, fasta_dict)

    write_fasta(fasta_dict)

main()