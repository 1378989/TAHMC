######################################################################
######################################################################
######################For host mapping results:
######
library(GenomicFeatures)
suppressWarnings(txdb <- makeTxDbFromGFF(file = '/mnt/data/CD/RNA_seq/ref/gencode.v40.annotation.gtf.gz',format = 'gtf',organism = 'Homo sapiens'))
keytypes(txdb)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
library(tximport)

#--------------------------------------------------------------------------------------------------------------------------------
#PRJNA248469
###################
samples <- list.files(path = "/mnt/data/CD/RNA_seq/data/Input/datasets/PRJNA248469/Salmon_outputs", full.names = T, pattern="output$")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "/mnt/data/CD/RNA_seq/data/Input/datasets/PRJNA248469/Salmon_outputs/", "") %>%
  str_replace("_output", "")
all(file.exists(files))


txi <- tximport(files, type = "salmon", tx2gene = tx2gene,ignoreAfterBar = T,countsFromAbundance="lengthScaledTPM")


data <- txi$abundance%>% 
  data.frame()

library(rtracklayer)
gtf <- rtracklayer::import('/mnt/data/CD/RNA_seq/ref/gencode.v40.annotation.gtf.gz')
gtf_df <- as.data.frame(gtf) 
dim(gtf_df)



##############
library(tidyr)
#gene_name, gene_id, gene_biotype
geneid_df <- dplyr::select(gtf_df,c(gene_name,gene_id,gene_type))
##
sort(table( geneid_df$gene_type ) )
##
length(unique( geneid_df$gene_id ) )
#
index <- duplicated(geneid_df$gene_id)
table(index)
##
geneid_df <- geneid_df[!index, ]
dim(geneid_df )
##
#save(geneid_df,file = "gene name-gene ID-gene biotype.Rdata")

##########################################
##
library(dplyr)
library(tidyr)

##
#load(file = "expr_df.Rdata")
head(data)
data<-data%>%rownames_to_column(var="gene_id")
## 
data$gene_id<-str_split(data$gene_id,"[.]",simplify = T)[,1]
expr_df<-data
expr_df <- expr_df[!(duplicated(expr_df$gene_id)),]
##
#load(file = "gene name-gene ID-gene biotype.Rdata")
head(geneid_df)
geneid_df<-geneid_df[!(is.na(geneid_df$gene_name)),]
geneid_df$gene_id<-str_split(geneid_df$gene_id,"[.]",simplify = T)[,1]
geneid_df<-geneid_df[!(grepl("ENSG",geneid_df$gene_name,)),]


##mRNA:
mRNA_exprSet <-geneid_df %>% 
  dplyr::filter(gene_type=="protein_coding") %>%
  dplyr::select(c(gene_name,gene_id,gene_type)) %>% 
  dplyr::inner_join(expr_df,by ="gene_id") %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_type,sep = " | ")
dim(mRNA_exprSet)
#save(mRNA_exprSet,file = "./Clear/mRNA_exprSet.Rdata")

##lncRNA
#sort(table(geneid_df$gene_type))
LncRNA_exprSet <- geneid_df %>% 
  dplyr::filter(gene_type=="lncRNA") %>%
  dplyr::select(c(gene_name,gene_id,gene_type)) %>% 
  dplyr::inner_join(expr_df,by ="gene_id") %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_type,sep = " | ")
#save(LncRNA_exprSet,file = "./Clear/LncRNA_exprSet.Rdata")

##########For mRNA matrix data:

gene_symbol<-(str_split(mRNA_exprSet$gene_id,"[ | ]",simplify = T))[,1]
head(gene_symbol)
mRNA_exprSet$gene_id<-gene_symbol
table(duplicated(gene_symbol))
mRNA_exprSet<-mRNA_exprSet[!(duplicated(mRNA_exprSet$gene_id)),]
mRNA_exprSet[1:4,1:4]
rownames(mRNA_exprSet)<-mRNA_exprSet$gene_id
mRNA_exprSet<-mRNA_exprSet[,-1]
save(mRNA_exprSet,file = "./mRNA_exprSet.Rdata")

########For mRNA matrix data:
gene_symbol<-(str_split(LncRNA_exprSet$gene_id,"[ | ]",simplify = T))[,1]
head(gene_symbol)
LncRNA_exprSet$gene_id<-gene_symbol
table(duplicated(gene_symbol))
LncRNA_exprSet<-LncRNA_exprSet[!(duplicated(LncRNA_exprSet$gene_id)),]
LncRNA_exprSet[1:4,1:4]
rownames(LncRNA_exprSet)<-LncRNA_exprSet$gene_id
LncRNA_exprSet<-LncRNA_exprSet[,-1]

save(LncRNA_exprSet,file = "./LncRNA_exprSet.Rdata")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
######################################################################
######################################################################
######################For microbiome mapping results:
#####

#######For genus level:
#---------------------------------------------------------------------------------------------------------------------------------------------------

read_kraken_reports = function(files, sample_names = NULL, study_name = NULL){
  require(data.table)
  require(tibble)
  
  if(is.null(sample_names)){sample_names = files}
  if(is.null(study_name)){study_name = NA}
  if(length(study_name) == 1){study_name = rep(study_name, length(files))}
  
  df = list()
  n = 0
  for(i in 1:length(files)){
    if(round(i/length(files)*100, 2) > n){n = round(i/length(files)*100, 2); cat(paste0('\r',n,'% done   '))}
    x = read.delim(files[i], header = F)
    df[[i]] = data.frame(study = study_name[i], sample = sample_names[i],
                         name = x$V1, reads = x$V2)
  }
  df = rbindlist(df) %>% tibble()
  cat('\n')
  df
}

#---------------------------------------------------------------------------------------------------------------------------------------------------
##################PRJNA248469################
reports = list()
f = list.files('/mnt/data/CD/RNA_seq/data/Input/datasets/PRJNA248469/kraken2', full.names = T) %>% str_subset('.kreport2')
s = f %>% str_remove('.*/mnt/data/CD/RNA_seq/data/Input/datasets/PRJNA248469/kraken2/') %>% str_remove('\\..*')
reports$PRJNA248469 = read_kraken_reports(f,s,'PRJNA248469')
data_kr = rbindlist(reports) %>% tibble()


data_kr<-data_kr[grepl("g__",data_kr$name),]
data_kr<-data_kr[!(grepl("s__",data_kr$name)),]
data_kr<-data_kr[!(grepl("g__Homo",data_kr$name)),]

taxonomy<-unique(data_kr$name)

data_kr$name<-str_extract(data_kr$name, "g__.+")

unique(data_kr$name)

####
data.new<-matrix(ncol = length(unique(data_kr$sample)),nrow = length(unique(data_kr$name)),data = 0)
rownames(data.new)<-unique(data_kr$name)
colnames(data.new)<-unique(data_kr$sample)
data.new<-data.new%>%as.data.frame%>%rownames_to_column("genus")

df<-data.new%>%reshape2::melt(id.vars = c("genus"), 
                              value.name = "value", variable.name = "sample")

df$id<-paste(df$sample,df$genus,sep = ",")
data_kr$id<-paste(data_kr$sample,data_kr$name,sep = ",")
df$value<-data_kr$reads[match(df$id,data_kr$id)]
df$value<-ifelse(is.na(df$value),0,df$value)
df$id<-NULL
df_new = pivot_wider(df, names_from = 'sample', values_from = 'value')
df_new<-df_new%>%column_to_rownames("genus")
save(df_new,file = "/mnt/data/CD/RNA_seq/data/Input/datasets/PRJNA248469/microbe_genus_count.Rdata")

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

###########
####For species level:


read_kraken_reports = function(files, sample_names = NULL, study_name = NULL){
  require(data.table)
  require(tibble)
  
  if(is.null(sample_names)){sample_names = files}
  if(is.null(study_name)){study_name = NA}
  if(length(study_name) == 1){study_name = rep(study_name, length(files))}
  
  df = list()
  n = 0
  for(i in 1:length(files)){
    if(round(i/length(files)*100, 2) > n){n = round(i/length(files)*100, 2); cat(paste0('\r',n,'% done   '))}
    x = read.delim(files[i], header = F)
    df[[i]] = data.frame(study = study_name[i], sample = sample_names[i],
                         name = x$V1, reads = x$V2)
  }
  df = rbindlist(df) %>% tibble()
  cat('\n')
  df
}

#---------------------------------------------------------------------------------------------------------------------------------------------------

######PRJNA248469########
reports = list()
f = list.files('/mnt/data/CD/RNA_seq/data/Input/datasets/PRJNA248469/kraken2', full.names = T) %>% str_subset('.kreport2')
s = f %>% str_remove('.*/mnt/data/CD/RNA_seq/data/Input/datasets/PRJNA248469/kraken2/') %>% str_remove('\\..*')
reports$PRJNA248469 = read_kraken_reports(f,s,'PRJNA248469')
data_kr = rbindlist(reports) %>% tibble()


data_kr<-data_kr[grepl("s__",data_kr$name),]
data_kr<-data_kr[!(grepl("s__Homo sapiens",data_kr$name)),]

taxonomy<-unique(data_kr$name)

data_kr$name<-str_extract(data_kr$name, "g__.+")

unique(data_kr$name)

####
data.new<-matrix(ncol = length(unique(data_kr$sample)),nrow = length(unique(data_kr$name)),data = 0)

micr=data.frame(micr=unique(data_kr$name))
micr$number=paste0("A",1:nrow(micr))
data_kr$number<-micr$number[match(data_kr$name,micr$micr)]

rownames(data.new)<-unique(data_kr$number)
colnames(data.new)<-unique(data_kr$sample)



data.new<-data.new%>%as.data.frame%>%rownames_to_column("species")

df<-data.new%>%reshape2::melt(id.vars = c("species"), 
                              value.name = "value", variable.name = "sample")

df$id<-paste(df$sample,df$species,sep = ",")
data_kr$id<-paste(data_kr$sample,data_kr$number,sep = ",")
df$value<-data_kr$reads[match(df$id,data_kr$id)]
df$value<-ifelse(is.na(df$value),0,df$value)
df$id<-NULL
df_new = pivot_wider(df, names_from = 'sample', values_from = 'value')
df_new$species<-data_kr$name[match(df_new$species,data_kr$number)]

save(df_new,file = "/mnt/data/CD/RNA_seq/data/Input/datasets/PRJNA248469/microbe_species_count.Rdata")

#--------------------------------------------------------------------------------------------------------------------------------------------------------------













