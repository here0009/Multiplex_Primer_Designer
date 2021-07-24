# coding: utf-8
# Read the combination report, return the information on demand, store the already order primer with a primer_name in the given primer_dict. If the all of the primer names are in the given set, Print the combination out, and Add the primer name apperance to combination_number_dict. Check the appearance value of primers in all the combinations.
import sys
report_name = sys.argv[1]
primer_each_group = 12
primer_set = {"Adenovirus_f-g#F#4o","Adenovirus_f-g#F#7","Adenovirus_f-g#F#7o","Adenovirus_f-g#R#1","Adenovirus_f-g#R#2","Adenovirus_f-g#R#6","Astrovirus_1-8#F#4n","Astrovirus_1-8#F#2o","Astrovirus_1-8#F#4o","Astrovirus_1-8#R#0","Astrovirus_1-8#R#1n","Astrovirus_1-8#R#2o","Norovirus_GI#F#1","Norovirus_GI#F#5o","Norovirus_GI#R#1o","Norovirus_GII#F#3","Norovirus_GII#R#1","Rotavirus_A#F#12o","Rotavirus_A#F#2o","Rotavirus_A#R#0","Sapovirus_GI#F#0","Sapovirus_GI#F#16o","Sapovirus_GI#F#2","Sapovirus_GI#R#4","Sapovirus_GI#R#1","Sapovirus_GI#R#2n"}
print(report_name)
with open(report_name) as fhand:
    combination_number_dict = dict() #store the primer and there appearance number in primer combinations.
    total_combination = 0
    for line in fhand:
        if line.startswith("Primer Combinations Without Interaction:"):
            length_info = fhand.readline().strip()
            primer_dict = dict()
            flag = True
            for i in range(primer_each_group):
                primer_name = fhand.readline().strip()[1:]
                if primer_name not in primer_set:
                    flag = False
                    break
                else:
                    primer_seq = fhand.readline().strip()
                    primer_dict[primer_name] = primer_seq
            if flag:
                total_combination += 1
                print(length_info)    
                for primer_name, primer_seq in primer_dict.items():
                    combination_number_dict[primer_name] = combination_number_dict.get(primer_name, 0) + 1
                    print('>' + primer_name)
                    print(primer_seq)
                print()

    print("The toatal combinations are :", total_combination)
    print("The primer appearance in combinations are :")
    for primer_name, combination_number in combination_number_dict.items():
        print(primer_name + '\t' + str(combination_number))
