###########################################################
#########Download

##Skip this step if the raw sequencing data is a fastq file
##If the raw sequencing data is an SRR file, run the following code


###Download SRR file:
nohup wget -c -i download.txt >/dev/null 2>&1 &


###SRR to fastq

#fastq-dump v3.0.2 

cd /mnt/data/CD/Bulk_RNA_seq/data1/raw_fastq
ls * |while read id;do ( fasterq-dump -e 16  ./${id});done
