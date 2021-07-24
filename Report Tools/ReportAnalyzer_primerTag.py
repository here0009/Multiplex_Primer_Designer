# coding: utf-8
# Read the combination report, return the information on demand, store the already order primer with a primer_name endswith o. Check the group contains the most already ordered primers.
import sys
report_name = sys.argv[1]
primer_each_group = 12

print(report_name)
with open(report_name) as fhand:
    max_order_primers = 7 #previously tested to be 7
    combination_number_dict = dict() #store the primer and there appearance number in primer combinations.

    for line in fhand:
        if line.startswith("Primer Combinations Without Interaction:"):
            order_primer_num = 0
            length_info = fhand.readline().strip()
            primer_dict = dict()
            for i in range(primer_each_group):
                primer_name = fhand.readline().strip()
                if primer_name.endswith('o'):
                    order_primer_num += 1
                primer_seq = fhand.readline().strip()
                primer_dict[primer_name] = primer_seq
            if order_primer_num >= max_order_primers:
                max_order_primers = order_primer_num
                print(length_info)
                print("The number of order primer is {}".format(max_order_primers))
                for primer_name, primer_seq in primer_dict.items():
                    combination_number_dict[primer_name] = combination_number_dict.get(primer_name, 0) + 1
                    print(primer_name)
                    print(primer_seq)
                print()

    print("The primer appearance in combinations are :")
    for primer_name, combination_number in combination_number_dict.items():
        print(primer_name)
        print(combination_number)
