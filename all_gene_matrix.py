import pandas as pd
import numpy as np
import sys
import gzip
import time
from datetime import timedelta
import gc
import os
input_sam=sys.argv[1]
fastqfile=sys.argv[2]
KUL=sys.argv[3]
BorN=sys.argv[4]

if BorN.strip()=='border':
    tissue_class='B'
if BorN.strip()=='tumor':
    tissue_class='T'
if BorN.strip()=='normal':
    tissue_class='N'
#species_name=sys.argv[3]
list_dir=os.listdir()
#my_dict : MSPname:[gene1,gene2,gene34,gene11234]
my_dict={}
#my_dict2 : gene1 : 123 gene2 : 441 gene4 :1123
my_dict2={}
#my_dict3={} MSPname : 1234 Mspname2 :2234(count)
my_dict3={}
#my_dict4 =MSPname: 123(how many gene count) or [gene1,gene2,gene3]
my_dict4={}
#my_dict5= gene : ['AAAAAATTTA'] ['AAAAAATTAA']
my_dict5={}

#msp_gene that is assinged set()
msp_gene=[]
aaa=pd.read_csv('final',sep='\t')
igc_species_name=aaa.columns

for i in igc_species_name:
    my_dict[i.strip()]=np.unique(aaa[i].astype('Int32').dropna().astype('str').tolist()).tolist()


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


for i in aaa:
    my_dict5[i]=[]



for i in my_dict.keys():
    my_dict3[i]=0
def species_count(name,allgene):
    for i in  list(set(allgene) & set([str(i) for i in aaa])):
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


my_dict3['rest']=num

msp_list=[]
if "/" in  input_sam:
    samfile_name=input_sam.split("/")[1].strip('.sam')+"_result"
    matrix_name=input_sam.split("/")[1].strip('.sam')+"_matrix"
else:
    samfile_name=str(input_sam).strip('.sam')+"_result"
    matrix_name=str(input_sam).strip('.sam')+"_matrix"

with open(samfile_name,'w') as file:
    for i,l in sorted(my_dict3.items(),key=lambda item: item[1]):
        if int(l)==0:
            continue
        msp_list.append(i)
        file.writelines([i,'\t',str(l),'\n'])

msp_list.remove('rest')


#print(msp_list)

#for i,l in my_dict4.items():
#    if int(l)==0:
#        continue
#    print(i,'\t',l)


start=time.process_time()
with gzip.open(fastqfile) as file:
    zxc=file.readlines()
end = time.process_time()



with open(input_sam) as file:
    n=0
    for i in file:
        n+=1
        match=i.split("\t")[2].strip()
        if match=="*":
            continue
        my_dict5[match]=my_dict5[match]+[str(zxc[1+4*(n-1)].decode('utf-8').strip())]
my_dict6={}

all_cell_barcode=[]
def umi_count(gene_name,barcode):
    dict_a={}
    dict_b={}
    global all_cell_barcode
    for i in barcode:
        cell_barcode=i[0:16]
        umi=i[16:26]
        all_cell_barcode+=[cell_barcode]
        if cell_barcode not in dict_a:
            dict_a[cell_barcode]=[umi]
        if cell_barcode in dict_a:
            dict_a[cell_barcode]=dict_a[cell_barcode]+[umi]
    for i,l in dict_a.items():
        dict_b[i]=len(set(l))
    my_dict6[gene_name]=list(dict_b.items())

#    for i,l in dict_b.items():
#        print(i,'\t',l)

del zxc
gc.collect()


for i,l in  my_dict5.items():
    umi_count(i,l)

all_cell_barcode=list(set(all_cell_barcode))



matrix=np.zeros((len(all_cell_barcode),len(aaa)),int)
matrix=pd.DataFrame(matrix,index=all_cell_barcode,columns=np.array(list(aaa)))
row=matrix.index
column=matrix.columns
#dict6 genename: [(cellid,umicount),(cellid,umicount),...]
lll=len(my_dict6)
n=0
for i,l in my_dict6.items():
    n+=1
    print(n/lll)
    for k,z in l:
        m=row.get_loc(k)
        tt=column.get_loc(i)
        matrix.iloc[[m],[tt]]=z
aaa=matrix.index.to_list()
aaa=[f'{KUL}_{tissue_class}_{i}-1' for i in aaa]
matrix.index=aaa
matrix.to_csv(f'{matrix_name}')


end = time.process_time()
print("Iter 완성: ", timedelta(seconds=end-start))
print(np.sum(matrix.values))
