import pandas as pd
import numpy as np
import sys
input_sam=sys.argv[1]
my_dict={}
my_dict2={}
my_dict3={}
my_dict4={}
msp_gene=[]
aaa=pd.read_csv('final',sep='\t')
igc_species_name=aaa.columns

for i in igc_species_name:
	my_dict[i.strip()]=np.unique(aaa[i].astype('Int32').dropna().tolist()).tolist()

#print(my_dict[list(my_dict.keys())[1]])
with open(input_sam) as file:
	for i in file:
		if i.startswith('@'):
			continue
		else:
			match_id=i.split('\t')[2].strip()
			if match_id=="*":
				continue
			if match_id in my_dict2:
				my_dict2[match_id]+=1
			elif match_id not in my_dict2:
				my_dict2[match_id]=1
#print(my_dict2)
#print(len(my_dict2))
#print('complete')
geneset=my_dict2.keys()
aaa=set(geneset)
#print(aaa)

for i in my_dict.keys():
	my_dict3[i]=0
def species_count(name,allgene):
#	print(allgene)
	for i in  list(set(allgene) & set([int(i) for i in aaa])):
		my_dict3[name]+=int(my_dict2[str(i)])
		my_dict4[name]=len(list(set(allgene) & set([int(i) for i in aaa])))
		msp_gene.append(i)


for i,l in sorted(my_dict.items()):
	species_count(i,l)


#print(msp_gene)
other_gene=list(aaa-set(msp_gene))
num=0
for i in other_gene:
	num+=my_dict2[i.strip()]




for i,l in sorted(my_dict3.items(),key=lambda item: item[1]):
	if int(l)==0:   
		continue
	print(i,'\t',l)
print('rest','\t',num)

#for i,l in my_dict4.items():
#	if int(l)==0:
#		continue
#	print(i,'\t',l)
