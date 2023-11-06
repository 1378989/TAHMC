########################################################################
############### Microbiome mapping
####Kraken2 v2.1.2

###For paired end:
cat ./cat|while read id
do  
kraken2 --db /mnt/data/data/kraken2_database --threads 40 --report /mnt/data/CD/Bulk_RNA_seq/data/kraken2/${id}.kreport2  --use-mpa-style --gzip-compressed --paired /mnt/data/CD/Bulk_RNA_seq/data/Clean_data/${id}_1_val_1.fq.gz /mnt/data/CD/Bulk_RNA_seq/data/Clean_data/${id}_2_val_2.fq.gz > /mnt/data/CD/Bulk_RNA_seq/data/kraken2/${id}.kraken2 ;
done 



####For single end:

cat ./down |while read id
do  
kraken2 --db /mnt/data/data/kraken2_database --threads 40 --report /mnt/data/CD/Bulk_RNA_seq/data/kraken2/${id}.kreport2  --use-mpa-style --gzip-compressed /mnt/data/CD/Bulk_RNA_seq/data/Clean_data/${id}_trimmed.fq.gz > /mnt/data/CD/Bulk_RNA_seq/data/kraken2/${id}.kraken2 ;
done 
