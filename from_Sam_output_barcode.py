import pandas as pd
import numpy as np
import sys
import sys
import gzip
import time
from datetime import timedelta

input_sam=sys.argv[1]
fastqfile=sys.argv[2]
species_name=sys.argv[3]

#my_dict : MSPname:[gene1,gene2,gene34,gene11234]
my_dict={}
#my_dict2 : gene1 : 123 gene2 : 441 gene4 :1123
my_dict2={}
#my_dict3={} MSPname : 1234 Mspname2 :2234(count)
my_dict3={}
#my_dict4 =MSPname: 123(how many gene count) or [gene1,gene2,gene3]
my_dict4={}
#msp_gene that is assinged set()
msp_gene=[]
aaa=pd.read_csv('final',sep='\t')
igc_species_name=aaa.columns

for i in igc_species_name:
	my_dict[i.strip()]=np.unique(aaa[i].astype('Int32').dropna().tolist()).tolist()


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

#aaa = mapping 된 모든 gene 
geneset=my_dict2.keys()
aaa=set(geneset)


for i in my_dict.keys():
	my_dict3[i]=0
def species_count(name,allgene):
	for i in  list(set(allgene) & set([int(i) for i in aaa])):
		my_dict3[name]+=int(my_dict2[str(i)])
		my_dict4[name]=len(list(set(allgene) & set([int(i) for i in aaa])))
		msp_gene.append(i)


for i,l in sorted(my_dict.items()):
	species_count(i,l)


#other gene = samfile 에서 맵핑된 모든 gene들중에 MSP로 assign 되지 않은 gene들
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

start=time.process_time()
with gzip.open(fastqfile) as file:
	zxc=file.readlines()
end = time.process_time()
print("Iter 완성: ", timedelta(seconds=end-start))


gene_pool=my_dict[species_name]
gene_pool=[str(i) for i in gene_pool]
print(gene_pool)

with open(input_sam) as file, open(species_name,'w') as file2:
	n=0
	for i in file:
		n+=1
		match=i.split("\t")[2].strip()
		if match=="*":
			continue
		try:
			if match in gene_pool:
				file2.writelines([zxc[1+4*(n-1)].decode('utf-8').strip(),'\n'])
				print(n)
		except:
			pass

print('complete')
