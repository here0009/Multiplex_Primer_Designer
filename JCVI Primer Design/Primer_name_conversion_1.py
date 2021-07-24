#modification:export primer name are changed from group_name#primer_index#F/R to goup_name#F/R#index, notice: the primers got the same index does not mean they are a pair.
import sys
def jcvi_primer_name_conversion(inputfile): #e.g. >Adenovirus_F~G_Hexon_00002_0002.r /begin=752 /end=775 /orientation=-1 /length=23
    outputfile = inputfile.split(".")[0] + "_name_converted.fasta" 
    fhand_output = open (outputfile,'w')
    group_name_set = set()
    with open(inputfile) as fhand_input:
        for line in fhand_input:
            if line.startswith('>'):
                words = line.split("/")
                primer_l_r = words[0].split(".")[1].strip()
                sub_words = words[0].split('_')
                group_name = '_'.join([sub_words[0][1:],sub_words[1]])
                if group_name not in group_name_set:
                    index_f = -1 #first apperance of group_name, set index to -1
                    index_r = -1
                    group_name_set.add(group_name)
                if primer_l_r == "l":
                    primer_F_R = "F"
                    index_f +=1
                    primer_name = '#'.join([group_name,primer_F_R,str(index_f)])
                else:
                    primer_F_R = "R"
                    index_r +=1
                    primer_name = '#'.join([group_name,primer_F_R,str(index_r)]) # primer index is composed of primer name ,f or r and  index letter
                primer_seq = fhand_input.readline().strip()
                fhand_output.write('>'+primer_name+'\n')
                fhand_output.write(primer_seq+'\n')
    fhand_output.close()

def primer3_plus_name_conversion(inputfile): #e.g.>Norovirus_GII_1_F or >SapoVirus_GI~IV_R
    outputfile = inputfile.split(".")[0] + "_name_converted.fasta" 
    fhand_output = open (outputfile,'w')
    with open(inputfile) as fhand_input:
        index = 0
        for line in fhand_input:
            if line.startswith('>'):
                words = line.split("_")
                index = words[2].strip()
                primer_F_R = words[-1].strip() #The first output primer of primer3plus got not index e.g. Norovirus_GII_F
                if index == primer_F_R:
                    index = str(0)
            group_name = '_'.join([words[0][1:],words[1]])
            primer_name = '#'.join([group_name,primer_F_R,index]).strip()
            primer_seq = fhand_input.readline().strip()
            fhand_output.write('>'+primer_name+'\n')
            fhand_output.write(primer_seq+'\n')
    fhand_output.close()

# inputfile = input("Enter the primer design results : ")
inputfile = sys.argv[1]
if sys.argv[2] == 'jcvi':
    jcvi_primer_name_conversion(inputfile)
elif sys.argv[2] == 'primer3_plus':
    primer3_plus_name_conversion(inputfile)
else:
    print('The input primer name format is not supported!')






