#-----------------------------------------------------------------

load("./Input/matrix/gene/merge/tpm_mRNA.Rdata")
########################Transcriptomic data preprocess:
########define a function



filter_genes <- function(genes,qt){
  genes.sd <- transform(as.data.frame(genes), SD=apply(as.data.frame(genes),1,sd, na.rm = TRUE))
  ## select top genes with high SD (~ variability) across samples
  SD_quantile <- quantile(genes.sd$SD) ## identical to summary(genes.sd$SD)
  SD_cutoff <- SD_quantile[qt] ## 2nd quantile -- 25th quantile.
  genes.sd <- genes.sd[order(genes.sd$SD, decreasing = T),]
  top.variable.genes <- rownames(genes.sd[genes.sd$SD > SD_cutoff,])
  ## subset these genes from gene table
  select <- which(rownames(genes) %in% top.variable.genes)
  genes <- genes[select,]
}


FUN =  function(xxx) {
  (sum(xxx==0)/length(xxx))<0.8 }



########### All  data:
load("/mnt/data/CD/RNA_seq/data/Input/matrix/gene/merge/tpm_LncRNA.Rdata")
dim(exp)
exp_tpm <- exp%>%.[apply(., 1, FUN),]%>%filter_genes(2)
dim(exp_tpm)
exp_tpm<-log2(exp_tpm+1)
rm(exp)

##load metadata:
metadata<-read.csv("/mnt/data/CD/RNA_seq/data/Supplement/metadata_new.csv")
identical(colnames(exp_tpm),metadata$sample)
meta<-data.frame(sample=colnames(exp_tpm))
meta$dataset<-metadata$source[match(meta$sample,metadata$sample)]
meta$group<-metadata$group[match(meta$sample,metadata$sample)]
rownames(meta)<-meta$sample
identical(rownames(meta),colnames(exp_tpm))

library(tinyarray)
##before：
draw_pca(exp =exp_tpm, group_list = factor(meta$group))
draw_pca(exp = exp_tpm, group_list = factor(meta$dataset))


##after：
library(sva)
mod <- model.matrix(~factor(meta$group))
expr_combat <- ComBat(dat = exp_tpm, batch = meta$dataset,mod = mod)
draw_pca(exp = expr_combat, group_list = factor(meta$group))
draw_pca(exp = expr_combat, group_list = factor(meta$dataset))
rm(exp_tpm)

CD_id<-meta[meta$group%in%"CD",]$sample
expr_combat<-expr_combat[,colnames(expr_combat)%in%CD_id]
dim(expr_combat)

identical(colnames(expr_combat),rownames(df_norm_all))
identical(colnames(expr_combat),rownames(df_norm_Bac))

rm(meta,metadata,mod,CD_id,data_type)



microbiome_all<-df_norm_Bac
expr_combat<-t(expr_combat)





###########################
########

m2 <- c()
pvalue <- c()
pvalue.man <- c()
microbiome_list <- list(microbiome_all)
data_type <- c( "all_combined")
procrustes_result <- data.frame(X1=NA,X2=NA,MDS1=NA,MDS2=NA,type=NA)
for (i in 1:1) {
  if(T){
    genes <-expr_combat
    genes[genes<0] <- 0
    microbes <- microbiome_list[[i]] %>% .[rownames(genes),]
    microbes[microbes<0] <- 0
    microbe.dist <- vegan::vegdist(microbes,method = 'bray')
    gene.dist <- vegan::vegdist(genes,method = 'bray')
    mds.m <- monoMDS(microbe.dist)
    mds.g <- monoMDS(gene.dist)
    
    pro.m.g <- procrustes(mds.g,mds.m,symmetric = TRUE)
    summary(pro.m.g)
    
    # M2 calculating
    set.seed(123)
    pro.m.g_t <- protest(mds.g,mds.m,permutations = 9999)
    # pro.m.g_t
    m2 <- c(m2,pro.m.g_t$ss)
    pvalue <- c(pvalue,round(pro.m.g_t$signif,5))
    
    Pro_Y <- cbind(data.frame(pro.m.g$Yrot),data.frame(pro.m.g$X))
    Y <- Pro_Y %>% dplyr::mutate(type=rep(data_type[i]))
    Pro_X <- data.frame(pro.m.g$rotation)
    procrustes_result <- rbind(procrustes_result,Y)
    
    # mantel analysis
    micro.dist <- vegan::vegdist(microbes,method = 'bray')
    gene.dist <- vegan::vegdist(genes,method = 'bray')
    set.seed(123)
    man.m.g <- mantel(micro.dist,gene.dist,method = 'spearman', permutations = 9999, na.rm = TRUE)
    pvalue.man <- c(pvalue.man,round(man.m.g$signif,5))
  }
}





result <- procrustes_result %>% .[-1,] %>% dplyr::mutate(type=as.factor(.$type)) %>% dplyr::mutate(type=factor(.$type,levels = c("all_combined")))
all<-data.frame(m2=m2,pvalue=pvalue,pvalue.man=pvalue.man)
write.table(result,file = "/mnt/data/CD/RNA_seq/data/Output/Others/07/Procrustes_mantel_result_all_CD_lncRNA.txt",sep = "\t")
write.table(all,file = "/mnt/data/CD/RNA_seq/data/Output/Others/07/Procrustes_mantel_all_all_CD_lncRNA.txt",sep = "\t")

# plot
result<-read.table(file = "/mnt/data/CD/RNA_seq/data/Output/Others/07/Procrustes_mantel_result_all_CD_lncRNA.txt",sep = "\t")
all<-read.table(file = "/mnt/data/CD/RNA_seq/data/Output/Others/07/Procrustes_mantel_all_all_CD_lncRNA.txt",sep = "\t")
p=ggplot(result) +
  geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + 
                                                                   MDS2)/2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#2F75B3", size = 1) +
  geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend 
                   = MDS2), 
               arrow = arrow(length = unit(0, 'cm')),
               color = "#EF856A", size = 1) +
  geom_point(aes(X1, X2,shape='23'), fill ="#2F75B3" , size = 2.8) +
  geom_point(aes(MDS1, MDS2,shape='21'), fill = "#EF856A", size = 0.4) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'),
        axis.ticks.length = unit(0.3,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_text(colour='black', size=14),
        axis.title.y=element_text(colour='black', size=14),
        axis.text=element_text(colour='black',size=12)) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  scale_shape_manual(name = "",
                     values = c('23' = 23,'21' = 21),
                     labels = c('Microbiota abundance(Bray-Curtis)', "Host gene expression(Aitchison's)"))+
  labs(title="Correlation between Host gene expression and microbiome abundance 
in LUSC") +
  annotate('text', label = 'Procrustes analysis:\n M2 = 0.946, p-value = 0.004',size = 4.2,hjust = 0) +
  theme(plot.title = element_text(size=14,colour = "black",hjust = 0,face = "bold"))+
  labs(x='Dimension 1',y='Dimension 2')+
  theme(legend.position = 'bottom')




pdf('/mnt/data/CD/RNA_seq/data/Output/Others/07/procrusters_analysis_lncRNA.pdf',width = 15,height =6)
p
dev.off()

m2