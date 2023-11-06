####################################################
###############Host mapping:
##Salmon v1.8.0


cd ../Clean_data
index=/mnt/data/project/genome_reference/Human_GRCh38_v40/salmon_index

##For paired end:
ls *gz|cut -d "_" -f 1 | sort -u | while read id;do (nohup salmon quant -i $index  -l A  --gcBias -1 ${id}*.sra_1_val_1.fq.gz   -2 ${id}*.sra_2_val_2.fq.gz  -o ../Salmon_outputs/${id}_output &);done

##For single end:
ls *gz| sort -u | while read id;do (nohup salmon quant -i $index  -l A  --gcBias -r ${id} -o ../Salmon_outputs/${id}_output  &);done
