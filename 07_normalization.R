###########################normalization:
#---------------------------------------------------------------------------------------------------------------------------------
###################For microbiome data:
##1. count to log2(CPM+1)
##2. batch remove


##count to log2(CPM+1)
load("./micro/merge/micro_count_raw_genus_all.Rdata")
Total_reads<-read.csv("./Total_reads/Total_reads.csv")

for (i in 1:ncol(micro)) {
  micro[,i]<-c(micro[,i]/Total_reads[i,2])*10^6
}
micro<-log2(micro+1)
#-------------------------------------------------------------------------------
####batch remove

metadata<-read.csv("./Supplement/metadata_datasets.csv")
load("./micro/merge/micro_filter_norm_genus.Rdata")
micro<-micro_norm_filter_genus_Bac
identical(colnames(micro),metadata$sample)
meta<-data.frame(sample=colnames(micro))
meta$dataset<-metadata$source[match(meta$sample,metadata$sample)]
meta$group<-metadata$group[match(meta$sample,metadata$sample)]
rownames(meta)<-meta$sample
identical(rownames(meta),colnames(micro))

library(tinyarray)
#####before：
draw_pca(exp = micro, group_list = factor(meta$group))
draw_pca(exp = micro, group_list = factor(meta$dataset))


######after：
library(sva)
mod <- model.matrix(~factor(meta$group))
expr_combat <- ComBat(dat = micro, batch = meta$dataset,mod = mod)
draw_pca(exp = expr_combat, group_list = factor(meta$group))
draw_pca(exp = expr_combat, group_list = factor(meta$dataset))

#################################
save(expr_combat,file = "./micro/merge/micro_rebatch_filter_genus_Bac.Rdata")





#---------------------------------------------------------------------------------------------------------------------------------
###################For host gene:

metadata<-read.csv("./Supplement/metadata_datasets.csv")
####：
load("./gene/merge/norm_mRNA.Rdata")
identical(colnames(exprset),metadata$sample)
meta<-data.frame(sample=colnames(exprset))
meta$dataset<-metadata$source[match(meta$sample,metadata$sample)]
meta$group<-metadata$group[match(meta$sample,metadata$sample)]
rownames(meta)<-meta$sample
identical(rownames(meta),colnames(exprset))

library(tinyarray)
##before：
draw_pca(exp = exprset, group_list = factor(meta$group))
draw_pca(exp = exprset, group_list = factor(meta$dataset))


##after：
library(sva)
mod <- model.matrix(~factor(meta$group))
expr_combat <- ComBat(dat = exprset, batch = meta$dataset,mod = mod)
draw_pca(exp = expr_combat, group_list = factor(meta$group))
draw_pca(exp = expr_combat, group_list = factor(meta$dataset))
#################################
save(expr_combat,file = "./gene/merge/rebatch_mRNA")
