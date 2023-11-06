######################################################################################
#### Merging microbiome raw count matrix data from all datasets:
###########genus level:

########A five-step filtering pipeline consisting of batch filtering, correlation filtering, prevalence and read count filtering, blacklist filtering, and manual literature inspection filtering

####(1)Batch filtering:
#datasets with fewer than 50 samples were excluded from this analysis.
#-----------------------------------------------------------------------------------------------------------------------------
require(foreach)
require(tidyverse)
require(ggplot2)
require(data.table)
require(doParallel)
require(compositions)
require(VennDiagram)
registerDoParallel(cores=9)

load("./micro/merge/micro_count_raw_genus_all.Rdata")
data.new<-micro_count_raw_genus_all%>%as.data.frame%>%rownames_to_column("genus")
rm(micro_count_raw_genus_all)
df<-data.new%>%reshape2::melt(id.vars = c("genus"), 
                              value.name = "value", variable.name = "sample")
meta<-read.csv("./Supplement/metadata_datasets.csv")

df$value[df$value>0]=1
df$value[df$value<=0]=0
df$dataset<-meta$source[match(df$sample,meta$sample)]
df$group<-meta$group[match(df$sample,meta$sample)]

#nsample>50:
df<-df%>%group_by(dataset)%>%
dplyr::mutate(nsample = n())%>%
subset(nsample > 50*length(unique(df$genus)))



d = data.frame()
for (i in 1:nrow(data.new)) {
  micro<-data.new$genus[i]
  df_pre<-df[df$genus%in%micro,]
  prev_stats<-df_pre%>%group_by(dataset) %>%summarise(prevalence = sum(value)/ n())
  max_prev <- max(prev_stats$prevalence)
  min_prev <- min(prev_stats$prevalence)
  fold<-max_prev/min_prev
  max_level <- pull(prev_stats, dataset)[which(prev_stats$prevalence == max(prev_stats$prevalence))][1]
  min_level <- pull(prev_stats, dataset)[which(prev_stats$prevalence == min(prev_stats$prevalence))][1]
  crumb <- tibble(taxa = micro, 
                  max_level = max_level,  max_prev = max_prev,
                   min_level = min_level,min_prev = min_prev,
                  fold= fold)
  
  
  d = rbind(d, data.frame(crumb))
}





##
#decontams<-data.frame()
#for (i in c(seq(0.15,0.50,0.05))) {
 # for (j in c(seq(1,2,0.1))) {
 #   ngenus<-sum(d$max_prev>i&d$fold>j)
  #  #percent<-sum(d$taxa[c(d$max_prev>i&d$fold>j)]%in%contamination)/length(d$taxa[c(d$max_prev>i&d$fold>j)])
#    decontam<-tibble(max_prev = i, fold_diff = j, n_contam =ngenus,
 #                    #percent=percent
  #                   )
 #   gusto_decontams= rbind(gusto_decontams, data.frame(decontam))
 # }
#}


#
contam<-d$taxa[d$max_prev>0.15&d$fold>1.5]
decontam<-d$taxa[!(d$max_prev>0.15&d$fold>1.5)]



####(2)Correlation filtering:
require(foreach)
require(tidyverse)
require(ggplot2)
require(data.table)
require(doParallel)
require(compositions)
require(VennDiagram)
registerDoParallel(cores=9)


##########transform to log(cpm+1)
#####Total reads
load("./micro/merge/micro_count_raw_genus_all.Rdata")
Total_reads<-read.csv("./Total_reads/Total_reads.csv")


micro<-micro_count_raw_genus_all
rm(micro_count_raw_genus_all)

for (i in 1:ncol(micro)) {
  micro[,i]<-c(micro[,i]/Total_reads[i,2])*10^6
}
micro<-log2(micro+1)

#

dave = data.frame()

#
suppressWarnings(for (contaminant_taxon in contam) {
  for (non_contaminant_taxon in decontam) {
  sample=meta[meta$source%in%d[d$taxa==contaminant_taxon,]$max_level,]$sample
  spearman_test <- cor.test(as.numeric(micro[contaminant_taxon, sample]), as.numeric(micro[non_contaminant_taxon, sample]))
  rho <- spearman_test$estimate
  morsel <- tibble( level = d[d$taxa==contaminant_taxon,]$max_level,
                                        non_contaminant_taxon = non_contaminant_taxon, 
                                        contaminant_taxon = contaminant_taxon,
                                        rho = rho)

  dave  = rbind(dave , data.frame(morsel))
  }
})



cor_contam<-data.frame()
for (i in seq(0.5,0.8,0.1)){
  for (j in seq(2,20,2)) {
    
  }
  dave.new<-na.omit(dave)
  dave.new<-dave.new[dave.new$rho>i,]
  length(unique(dave.new$non_contaminant_taxon))
  cor_data<-tibble(cor=i,
                   #frequence=sum(unique(dave.new$non_contaminant_taxon)%in%contamination)/length(unique(dave.new$non_contaminant_taxon))
                   length=length(unique(dave.new$non_contaminant_taxon)))
  cor_contam  = rbind(cor_contam , data.frame(cor_data))
                      
}


#
dave.new<-na.omit(dave)
dave.new<-dave.new[dave.new$rho>0.8,]
length(unique(dave.new$non_contaminant_taxon))
cor_contam<-as.data.frame(table(dave.new$non_contaminant_taxon))
cor_contam<-as.character(cor_contam[cor_contam$Freq>0,]$Var1)



####(3)Prevalence and read count filtering:
####
#read count<100:

load("./micro/merge/micro_count_raw_genus_all.Rdata")
micro_reads<-rowSums(micro_count_raw_genus_all)%>%as.data.frame()
colnames(micro_reads)<-"reads"
micro_reads$genus<-rownames(micro_reads)
min_reads<-micro_reads[micro_reads$reads<=100,]$genus
####
#prevalence <5%:

load("./micro/merge/micro_count_raw_genus_all.Rdata")
data.new<-micro_count_raw_genus_all%>%as.data.frame%>%rownames_to_column("genus")
rm(micro_count_raw_genus_all)
df<-data.new%>%reshape2::melt(id.vars = c("genus"), 
                              value.name = "value", variable.name = "sample")
meta<-read.csv("./Supplement/metadata.csv")

df$value[df$value>0]=1
df$value[df$value<=0]=0
df$dataset<-meta$source[match(df$sample,meta$sample)]
df$group<-meta$group[match(df$sample,meta$sample)]

#nsample>50:
df<-df%>%group_by(dataset)%>%
dplyr::mutate(nsample = n())%>%
subset(nsample > 50*length(unique(df$genus)))

d = data.frame()
for (i in 1:nrow(data.new)) {
  micro<-data.new$genus[i]
  df_pre<-df[df$genus%in%micro,]
  prev_stats<-df_pre%>%group_by(dataset) %>%summarise(prevalence = sum(value)/ n())
  max_prev <- max(prev_stats$prevalence)
  min_prev <- min(prev_stats$prevalence)
  fold<-max_prev/min_prev
  max_level <- pull(prev_stats, dataset)[which(prev_stats$prevalence == max(prev_stats$prevalence))][1]
  min_level <- pull(prev_stats, dataset)[which(prev_stats$prevalence == min(prev_stats$prevalence))][1]
  crumb <- tibble(taxa = micro, 
                  max_level = max_level,  max_prev = max_prev,
                   min_level = min_level,min_prev = min_prev,
                  fold= fold)
  
  
  d = rbind(d, data.frame(crumb))
}

low_contam<-d[d$max_prev<0.05,]$taxa



####(4)Blacklist filtering:

contam<-read.csv("./Supplement/contaminant.csv")
#contam$genera_studies<-paste0("g__",contam$genera_studies)
contam$genera_Salter<-paste0("g__",contam$genera_Salter)

####(5)Manual literature inspection filtering:

####To further fully remove potential contaminants from the data, following the four steps listed above we then applied a manual literature search to remove other microbial genera that were unlikely to be present in the samples.
