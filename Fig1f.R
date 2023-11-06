# load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(psych)
library(ggpubr)
library(vegan)
library(reshape2)



# load metadata 
df_metadata <- read.csv("/mnt/data/project/Input/Poore/Kraken/Metadata-TCGA-Kraken-17625-Samples.csv")
colnames(df_metadata)[1]<-"id"
# only included primary tumor
df_metadata_PT <- df_metadata %>% filter(sample_type=="Primary Tumor")



# load microbial data from Poore's study

df_normalization_filter_likely <- read.csv("/mnt/data/project/Input/Poore/Kraken/Kraken-TCGA-Voom-SNM-Likely-Contaminants-Removed-Data.csv")

# include primary tumor

df_normalization_filter_likely_PT <- df_normalization_filter_likely %>% filter(X %in% df_metadata_PT$id) %>% remove_rownames %>% column_to_rownames(var="X")


# fish out patients containing both WGS and RNA-Seq
df_patients <- unique(df_metadata_PT[, c("case_uuid", "experimental_strategy")])
df_patients$type <- NA
for (i in 1:nrow(df_patients)){
  df_sub <- df_patients %>% filter(case_uuid==df_patients$case_uuid[i])
  number_type <- length(unique(df_sub$experimental_strategy))
  df_patients$type[i] <- number_type
}
df_patients_two <- df_patients %>% filter(type==2)
pateints_two <- unique(df_patients_two$case_uuid)

# fish out patients containing more than one WGS
df_patients_WGS <-df_metadata_PT %>% filter(experimental_strategy=="WGS")
df_patients_WGS_no <- table(df_patients_WGS$case_uuid) %>% as.data.frame()

df_patients_WGS_two <- df_patients_WGS_no %>% filter(Freq>1)
pateints_WGS_two <- unique(df_patients_WGS_two$Var1)



# compare microbota between WGS and RNA-Seq or between WGS itself based on Bray-Curtis distance
################################################################################################################
dist2list <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(inDist))
}

# WGS vs. RNA-Seq for individuals
function_bray <- function(in_matrix){
  v_cor <- vector()
  for (patient in pateints_two){
    meta_sub <- df_metadata_PT %>% filter(case_uuid==patient)
    df_data <- in_matrix[meta_sub$id, ]
    df_data[df_data<0] <- 0
    distance <- vegdist(df_data, method = 'bray')
    df_dist <- dist2list(distance)
    df_dist_m1 <- merge(df_dist, meta_sub[,c("id","experimental_strategy")], by.x="row", by.y = "id", all.x = TRUE)
    df_dist_m2 <- merge(df_dist_m1, meta_sub[,c("id","experimental_strategy")], by.x="col", by.y = "id", all.x = TRUE)
    df_dist_m2$filter <- ifelse(df_dist_m2$experimental_strategy.x!=df_dist_m2$experimental_strategy.y, "Y", "X")
    df_dist_m2 <- df_dist_m2 %>% filter(filter=="Y")
    v_cor <- append(v_cor, df_dist_m2$value)
  }
  return(v_cor)
}


v_vor_nor_f_likely_bray <- function_bray(df_normalization_filter_likely_PT)


# WGS for individuals
function_bray_WGS <- function(in_matrix){
  v_cor <- vector()
  for (patient in pateints_WGS_two){
    meta_sub <- df_metadata_PT %>% filter(case_uuid==patient)
    df_data <- in_matrix[meta_sub$id, ]
    df_data[df_data<0] <- 0
    distance <- vegdist(df_data, method = 'bray')
    df_dist <- dist2list(distance)
    v_cor <- append(v_cor, df_dist$value)
  }
  return(v_cor)
}


v_vor_nor_f_likely_bray_WGS <- function_bray_WGS(df_normalization_filter_likely_PT)

# combine results

df_bray_nor_f_likely <- data.frame(bray=c(v_vor_nor_f_likely_bray,v_vor_nor_f_likely_bray_WGS),
                                   class=c(rep("WGS_RNASeq", length(v_vor_nor_f_likely_bray)), rep("WGS", length(v_vor_nor_f_likely_bray_WGS))))



# visualization
function_plot_b <- function(in_matrix){
  pd <- ggplot(in_matrix, aes(x=bray, fill=class)) +
    geom_density(alpha=0.6) +
    theme_classic()+
    scale_fill_manual(values = c("#BEB8DC","#96C37D"))+
    ylab("Density") +
    xlab("Bray-Curtis dissimilarity")+
    theme(axis.title = element_text(size = 18, colour="black"),
          axis.text = element_text(size = 16, colour = "black"),
          legend.text = element_text(size = 16),
          legend.position = "bottom")
  pb <- ggplot(in_matrix, aes(x=class, y=bray, fill=class)) +
    geom_boxplot(outlier.size=0.1) +
    #theme_base()+
    scale_fill_manual(values = c("#BEB8DC","#96C37D"))+
    coord_flip()+
    theme(panel.background = element_blank())+
    clean_theme()
  p_m <- ggarrange(pb, pd,
                   ncol = 1, nrow = 2,  align = "hv", 
                   widths = c(1, 1), heights = c(1, 2.5),
                   common.legend = TRUE)
  return(p_m)
}


p_bray_f_likely <- function_plot_b(df_bray_nor_f_likely)


pdf(file = "./Comparison_Bray–Curtis_dissimilarity_WGS_RNASeq.pdf", width = 12, height = 8)
library(patchwork)
p_bray_f_likely
dev.off()




#------------------------------------------------------------------------------------------------------------------------------------------------------# compare the normalized overall microbiota between WGS their self
################################################################################################################
function_r_WGS <- function(in_matrix, method){
  v_cor <- vector()
  for (patient in pateints_WGS_two){
    meta_sub <- df_metadata_PT %>% filter(case_uuid==patient)
    df_sub <- in_matrix[meta_sub$id, ] %>% t() %>% as.data.frame()
    m <- corr.test(df_sub, method=method)
    co_eff <- as.vector(m$r)
    co_eff <- co_eff[co_eff!=1]
    v_cor <- append(v_cor, co_eff)
  }
  return(v_cor)
}


v_vor_nor_f_likely_pearson_WGS <- function_r_WGS(df_normalization_filter_likely_PT, "pearson")



df_cor_nor_f_likely <- data.frame(pearson=c(v_vor_nor_f_likely_pearson,v_vor_nor_f_likely_pearson_WGS),
                                  class=c(rep("WGS_RNASeq", length(v_vor_nor_f_likely_pearson)), rep("WGS", length(v_vor_nor_f_likely_pearson_WGS))))

# visualization
function_plot_p <- function(in_matrix){
  pd <- ggplot(in_matrix, aes(x=pearson, fill=class)) +
    geom_density(alpha=0.6) +
    theme_classic()+
    scale_fill_manual(values = c("#db6968","#606f8a"))+
    ylab("Density") +
    xlab("Pearson correlation coefficient")+
    theme(axis.title = element_text(size = 18, colour="black"),
          axis.text = element_text(size = 16, colour = "black"),
          legend.text = element_text(size = 16),
          legend.position = "bottom")
  pb <- ggplot(in_matrix, aes(x=class, y=pearson, fill=class)) +
    geom_boxplot(outlier.size=0.1) +
    #theme_base()+
    scale_fill_manual(values = c("#db6968","#606f8a"))+
    coord_flip()+
    theme(panel.background = element_blank())+
    clean_theme()
  p_m <- ggarrange(pb, pd,
                   ncol = 1, nrow = 2,  align = "hv", 
                   widths = c(1, 1), heights = c(1, 2.5),
                   common.legend = TRUE)
  return(p_m)
}


p_nor_f_likely <- function_plot_p(df_cor_nor_f_likely)


pdf(file = "./Comparison_correlation_WGS_RNASeq.pdf", width = 12, height = 8)
p_nor_f_likely
dev.off()
























