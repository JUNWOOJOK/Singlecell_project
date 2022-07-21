import sys
import numpy as np
input_file=sys.argv[1]
my_dict={}
n=0
with open(input_file) as file:
	for i in file:
		umi=i[16:26].strip()
		cellbarcode=i[0:16].strip()
		if cellbarcode not in my_dict.keys():
			my_dict[cellbarcode]=[umi]
		if cellbarcode in my_dict.keys():
			my_dict[cellbarcode]=my_dict[cellbarcode]+[umi]

my_dict2={}
for i,l in my_dict.items():
	my_dict2[i]=len(np.unique(l)) 

for i,l in my_dict2.items():
	print(f'{i}-1'.strip(),'\t',str(l).strip(),sep='')
