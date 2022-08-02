import sys
import gzip 
R1_file=sys.argv[1]
R2_file=sys.argv[2]
R2_name=str(R2_file)

with gzip.open(R1_file) as file:
        zxc=file.readlines()
n=0
with  gzip.open(R2_file) as file2,open(f'manifold_{R2_name}','w') as file3:
	for i in file2:
		n+=1
		if n%4==1:
			file3.writelines([i.decode('utf-8').strip(),'\t',zxc[n].decode('utf-8').strip(),'\n'])
		else:
			file3.writelines([i.decode('utf-8').strip(),'\n'])
