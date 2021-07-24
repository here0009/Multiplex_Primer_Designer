# coding: utf-8
# Input: the sequence name and sequence value, the length of the output region, the length of the most conserved region.
# Process: Find the most conserved region of the input sequences, if the sequence is shorter than conserved_len, return the sequence directly.
# Output : the most conserved region of the sequence
import string
import sys

conserved_len = 150 #The length of the most conserved region
ajacent_len = 100 #The length next to the most conserved region
ouput_len = conserved_len + 2*ajacent_len

ATGC = {'A','T','G','C'}

input_file = sys.argv[1]
# input_file = input("Enter the file name :")


#Store the sequence information into a dictionary
seq_dict = dict()
with open(input_file) as fhand:
    for line in fhand:
        if line.startswith('>'):
            seq_name = line.strip()
        else:
            seq_value = line.strip()
            seq_dict[seq_name] = seq_value

#print(seq_dict)

def ambiguous_frequency(seq):
    seq = str.upper(seq)
    atgc_number = 0
    for letter in seq:
        if letter in ATGC:
            atgc_number+=1
    frequency = 1 - atgc_number/len(seq)
    return frequency

def conserved_region(seq):
	conserved_frequency = 1
	if len(seq)<= conserved_len:
		region = seq
		conserved_pos = 0
	else:
		for i in range(len(seq)-conserved_len):
			frequency = ambiguous_frequency(seq[i:i+conserved_len])
			if frequency < conserved_frequency:
				conserved_frequency = frequency
				region = seq[i:i+conserved_len]
				conserved_pos = i
	return conserved_pos

#Write the conserved region to a file
output_file = input_file.split('.')[0]+'_conserved_region.fasta'
sys.stdout = open(output_file, 'w')
output_region_dict = dict() #store the output region used to design primers
conserved_dict = dict()   #store the conserved region used for report
for seq_name,seq_value in seq_dict.items():
	conserved_pos = conserved_region(seq_value)
	conserved_dict[seq_name] = seq_value[conserved_pos:conserved_pos+conserved_len]
	if len(seq_value)<= ouput_len:
		output_region_dict[seq_name] = seq_value
	else:
		left_margin = conserved_pos - ajacent_len
		right_margin = conserved_pos+conserved_len+ajacent_len
		if left_margin<0: #exceed the left margin
			compensation_len = abs(left_margin)
			output_region_dict[seq_name] = seq_value[:right_margin+compensation_len]
		elif right_margin > len(seq_value): #exceed the right margin
			compensation_len = abs(right_margin-len(seq_value))
			output_region_dict[seq_name] = seq_value[left_margin-compensation_len:]
		else:
			output_region_dict[seq_name] = seq_value[left_margin:right_margin]

for seq_name,seq_value in output_region_dict.items():
	print(seq_name)
	print(seq_value)

sys.stdout.close()

#Write the conserved region to a file, with ambiguous_frequency
output_report = input_file.split('.')[0]+'_conserved_report.fasta'
sys.stdout = open(output_report, 'w')  
for seq_name,region in output_region_dict.items():
	frequency_conserved = ambiguous_frequency(conserved_dict[seq_name])
	frequency_all = ambiguous_frequency(region)
	print('The conserved region is: ')
	print('%s %.2f'%(seq_name,frequency_conserved))
	print(conserved_dict[seq_name])
	print('The output region is: ')
	print('%s %.2f'%(seq_name, frequency_all))
	print(region)
	print()
sys.stdout.close()
