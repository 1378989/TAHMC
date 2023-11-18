A new pipeline, named transcriptome analysis of host-microbiome crosstalk (TAHMC), was proposed for restoring both host gene expression and microbial quantification from bulk RNA-seq concurrently.







The procedures of TAHMC

1)Bulk RNA-seq data processing and mapping：

2)Microbial data decontamination and normalization：



We applied TAHMC to a total of 1855 samples (1554 CD samples and 301 control samples) of RNA-Seq data。





The proposed TAHMC accurately restored both host gene expression and microbial quantification from CD RNA-seq data, thereby revealing the potential causal associations between changes in microbial composition, diversity within CD mucosal tissue, and host gene expression disorders. Furthermore, TAHMC could be generalized and applied to other organs and diseases for investigating the association between tissue-resident microbes and host gene expression with general interesting.



#########Download

##Skip this step if the raw sequencing data is a fastq file
##If the raw sequencing data is an SRR file, run the following code


###Download SRR file:
nohup wget -c -i download.txt >/dev/null 2>&1 &


###SRR to fastq

#fastq-dump v3.0.2 

cd /mnt/data/CD/Bulk_RNA_seq/data1/raw_fastq
ls * |while read id;do ( fasterq-dump -e 16  ./${id});done
