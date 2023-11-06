##################################
############
#############quality checks and adaptor trimming

##Trim_galore v0.6.10

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

#——————————————————————————————————————————————————————————————
##########################################################################################
########FastQC v0.12.1 

mkdir ../fastqc

fastqc *fastq.gz -o ../fastqc -t 16

multiqc ../fastqc -o ../fastqc/multiqc
