import os
my_dict={}
file_list=os.listdir()
file_list=[i for i in file_list if 'result' in i and 'process' not in i]
for i in file_list:
	with open(i) as file:
		for line in file:
			id=line.split('\t')[0].strip()
			count=int(line.split('\t')[1].strip())
			if id in my_dict:
				my_dict[id]+=count
			elif id not in my_dict:
				my_dict[id]=count


for i,l in my_dict.items():
	print(i,'\t',l)

