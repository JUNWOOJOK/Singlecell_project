import os
dir_list=os.listdir()
dir_list=[i for i in dir_list if 'sam5' in i]
my_dict={}
for i in dir_list:
	with open(i) as file:
		for line in file:
			barcode=line.split("\t")[0]
			count=line.split("\t")[1]
			if barcode.strip() not in my_dict:
				my_dict[barcode.strip()]=int(count)
			else:
				my_dict[barcode.strip()]=int(my_dict[barcode.strip()])+int(count)

for i,l in my_dict.items():
	print(i,"\t",l,sep="")

