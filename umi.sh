
for i in $(ls|grep fastq|grep R1);do python ~/final_singlecell_process_for_all.py $(echo $i|cut -d "_" -f 1-4).sam $i ${i}_all_barcode;done

for i in $(ls|grep barcode);do python ~/umicout.py $i > $(echo $i|sed s/.fastq.gz_all_barcode//g).sam5;done
