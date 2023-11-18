TAHMC: Transcriptome Analysis of Host-Microbiome Crosstalk
================
Huijun Chang
6/11/2023

## Introduction
We propose a novel pipeline, dubbed Transcriptome Analysis of Host-Microbiome Crosstalk (TAHMC), designed to concurrently restore both host gene expression and microbial quantification from bulk RNA-seq data.

The procedures of TAHMC:
-  1) Bulk RNA-seq data processing and mapping
-  2) Microbial data decontamination and normalization
-  (1) Batch filtering:
-  (2) Correlation filtering:
-  (3) Prevalence and read count filtering:
-  (4) Blacklist filtering:
-  (5) Manual literature inspection filtering:

Please see the related manuscript for more information.


### Using the codes:
-  01.Download
-  02.Quality checks (QC)
-  03.Host mapping
-  04.Microbiome mapping
-  05.Mapping results compilation
-  06.decontamination
-  07.normalization

##  01.Download
``` bash
##Skip this step if the raw sequencing data is a fastq file. If the raw sequencing data is an SRR file, run the following code
###Download SRR file:
nohup wget -c -i download.txt >/dev/null 2>&1 &
###SRR to fastq
#fastq-dump v3.0.2 
cd /mnt/data/CD/Bulk_RNA_seq/data1/raw_fastq
ls * |while read id;do ( fasterq-dump -e 16  ./${id});done
```


##  02.QC
``` bash
#quality checks and adaptor trimming
########Trim_galore v0.6.10
#####For paired end:
ls *gz|cut -d"_" -f 1 | sort -u |  \
	while read id;do \
  nohup trim_galore -q 25 --phred33 --length 36 \
  --stringency 3 --paired -o ../Clean_data *${id}*.gz & \
  done

#####For single end:
ls *fastq.gz > config
paste config > 1
cat config
vim qc.sh

bin_trim_galore=trim_galore
dir='../Clean_data'
cat $1 |while read id
do
	    arr=(${id})
	    fq1=${arr[0]}
	    nohup $bin_trim_galore -q 25 --phred33 --length 36 --stringency 3 -o $dir $fq1 &
done
bash qc.sh config

########FastQC v0.12.1 
mkdir ../fastqc
fastqc *fastq.gz -o ../fastqc -t 16
multiqc ../fastqc -o ../fastqc/multiqc
```



##  03.Host mapping
``` bash
###############Host mapping:
##Salmon v1.8.0
cd ../Clean_data
index=/mnt/data/project/genome_reference/Human_GRCh38_v40/salmon_index
##For paired end:
ls *gz|cut -d "_" -f 1 | sort -u | while read id;do (nohup salmon quant -i $index  -l A  --gcBias -1 ${id}*.sra_1_val_1.fq.gz   -2 ${id}*.sra_2_val_2.fq.gz  -o ../Salmon_outputs/${id}_output &);done

##For single end:
ls *gz| sort -u | while read id;do (nohup salmon quant -i $index  -l A  --gcBias -r ${id} -o ../Salmon_outputs/${id}_output  &);done
```





##  04.Microbiome mapping
``` bash
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
```

The following codes are R scripts:

##  05.Mapping results compilation

##  06.decontamination

##  07.normalization





## Application

We applied TAHMC to a total of 1855 samples (1554 CD samples and 301 control samples) of RNA-Seq data.

The proposed TAHMC accurately restored both host gene expression and microbial quantification from CD RNA-seq data, thereby revealing the potential causal associations between changes in microbial composition, diversity within CD mucosal tissue, and host gene expression disorders. Furthermore, TAHMC could be generalized and applied to other organs and diseases for investigating the association between tissue-resident microbes and host gene expression with general interesting.




