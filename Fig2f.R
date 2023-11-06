rm(list = ls())
library(tidyverse)
library(microbiomeMarker)
load("./micro/merge/micro_count_filter_genus_Bac.Rdata")
metadata<-read.csv("./Supplement/metadata_new.csv")

micro<-micro_count_filter_genus_Bac
rm(micro_count_filter_genus_Bac)
identical(colnames(micro),metadata$sample)
meta<-data.frame(sample=colnames(micro))
meta$dataset<-metadata$source[match(meta$sample,metadata$sample)]
meta$group<-metadata$group[match(meta$sample,metadata$sample)]
rownames(meta)<-meta$sample
identical(rownames(meta),colnames(micro))
meta$group<-factor(meta$group,levels = c("Control","CD"))
meta<-meta%>%remove_rownames()%>%column_to_rownames("sample")
meta

library(tinyarray)
draw_pca(exp = micro, group_list = factor(meta$group))
draw_pca(exp =micro, group_list = factor(meta$dataset))

######batch remove：
micro<-as.matrix(micro)
combat_exp <- ComBat_seq(micro, batch=meta$dataset, group=meta$group)
combat_exp <- as.data.frame(combat_exp)
head(combat_exp)
micro<-combat_exp


#
load("./micro/taxonomy.Rdata")
taxonomy$genus<-str_extract(taxonomy$taxonomy, "g__.+")
taxonomy_new<-taxonomy$taxonomy[match(rownames(micro),taxonomy$genus)]%>%as.data.frame()
colnames(taxonomy_new)<-"tax"

table<-matrix(nrow = 653,ncol = 6)
colnames(table)<-c("Kingdom","Phylum","Class","Order","Family","Genus")
table<-as.data.frame(table)
table$Kingdom<-"k__Bacteria"
table$Phylum<-str_extract(taxonomy_new$tax, "p__.+")
table$Class<-str_extract(taxonomy_new$tax, "c__.+")
table$Order<-str_extract(taxonomy_new$tax, "o__.+")
table$Family<-str_extract(taxonomy_new$tax, "f__.+")
table$Genus<-str_extract(taxonomy_new$tax, "g__.+")

table$Phylum<-str_split(table$Phylum,"[|]",simplify = T)[,1]
table$Class<-str_split(table$Class,"[|]",simplify = T)[,1]
table$Order<-str_split(table$Order,"[|]",simplify = T)[,1]
table$Family<-str_split(table$Family,"[|]",simplify = T)[,1]
table$Genus<-str_split(table$Genus,"[|]",simplify = T)[,1]



identical(table$Genus,rownames(micro))
rownames(micro)<-paste0("OTU",1:653)
rownames(table)<-paste0("OTU",1:653)
table<-as.matrix(table)


###############---------------------------------

library(phyloseq)
physeq <- phyloseq(
  otu_table(micro,taxa_are_rows = TRUE), 
  tax_table(table), 
  sample_data(meta)
)
physeq



library(microbiomeMarker)
set.seed(12345)
data <- normalize(physeq, method = "rarefy")
data



adj <- tax_table(data) %>% apply(.,2,function(x) length(unique(x))) %>% sum

lefse <- run_lefse(
  data,  
  norm = "CPM",
  group = "group",
  multigrp_strat = TRUE,
  wilcoxon_cutoff = 0.05,
  kw_cutoff = 0.05, 
  lda_cutoff = 4
) 


lefse %>% marker_table() -> res.diff

dim(res.diff) 
res.diff %>% head()



cols <- RColorBrewer::brewer.pal(8, "Dark2")

res.diff$feature

res.diff$bac<-c("p__Firmicutes","o__Eubacteriales","c__Clostridia","f__Lachnospiraceae",
                "f__Bacteroidaceae","f__Oscillospiraceae","g__Phocaeicola",
                "g__Mycolicibacterium","g__Faecalibacterium","c__Bacteroidia","p__Proteobacteria",
                "c__Gammaproteobacteria","p__Actinobacteria","c__Actinomycetia"
                )


res.diff$ef_lda<-ifelse(res.diff$enrich_group%in%"Control",c(-res.diff$ef_lda),res.diff$ef_lda)
res.diff$type<-str_extract(res.diff$feature, "p__.+")
res.diff$type<-str_split(res.diff$type,"[|]",simplify = T)[,1]
data<-res.diff


p=ggplot(data)+
  geom_col(aes(reorder(bac, ef_lda), ef_lda, fill = type))+
  scale_fill_manual(values = c("#F57F20","#2078B4","#34A048","#E21F26","#8481BA","#E380AC"))+
  geom_segment(aes(y = 0, yend = 0,x = 0, xend = 12.4))+
  theme_classic()+
  ylim(-6,6)+
  coord_flip()+
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
  )+
  ylab("ef_lda")+
  geom_text(data=data[which(data$ef_lda > 0), ],aes(x = bac, y = 0, label = bac), 
            hjust = 1.2, size = 3)+
  geom_text(data=data[which(data$ef_lda < 0), ],aes(x = bac, y = 0, label = bac), 
            hjust = -0.1, size = 3)+
  ggtitle("Gene Ontology")+
  scale_x_discrete(expand=expansion(add=c(0,2)))+

  geom_segment(aes(y = -0.5, yend = -5,x = 18.3, xend = 18.3),
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), 
               size = 0.5)+
  geom_segment(aes(y = 0, yend = 5,x = 18.3, xend = 18.3),
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), 
               size = 0.5)+
  annotate("text", x = 5, y = -5, label = "Control")+
  annotate("text", x = 15, y =5, label = "CD")+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  theme(panel.grid = element_line(color = 'black', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank())+
  theme(axis.title = element_text(size = 18, colour="black"),
        axis.text = element_text(size = 16, colour = "black"),
        legend.text = element_text(size = 16))



pdf('/mnt/data/IBD/RNA_seq/data/Output/Others/08/Lefse1.pdf',width = 12,height = 12)
p
dev.off()
