

#-------load data:

############## Input datasets ##############

load(file = "/mnt/data/CD/RNA_seq/data/Output/Others/18/IBD_case_mRNA.Rdata")
case_mRNA<-IBD.lasso.stabsel[,c(1,2,14)]
load(file = "/mnt/data/IBD/RNA_seq/data/Output/Others/18/IBD_control_mRNA.Rdata")
control_mRNA<-IBD.lasso.stabsel[,c(1,2,14)]

load(file = "/mnt/data/CD/RNA_seq/data/Output/Others/18/IBD_case_LncRNA.Rdata")
case_LncRNA<-IBD.lasso.stabsel[,c(1,2,14)]
load(file = "/mnt/data/IBD/RNA_seq/data/Output/Others/18/IBD_control_LncRNA.Rdata")
control_LncRNA<-IBD.lasso.stabsel[,c(1,2,14)]

rm(IBD.lasso.stabsel)



unique(case_LncRNA$gene)
unique(case_LncRNA$taxa)

unique(case_mRNA$gene)
unique(case_mRNA$taxa)

########################################mRNA##################################################
#------------------mRNA__CD--------------

node<-data.frame(Part=c(rep("Genes",length(unique(case_mRNA$gene))),rep("Microbes",length(unique(case_mRNA$taxa)))),
                    node=c(unique(case_mRNA$gene),unique(case_mRNA$taxa)))

node$node<-gsub("g__","",node$node)


#large:

diameter.1 <- 8
nodes.1<-node%>%dplyr::filter(Part != "Microbes") %>% 
  sample_n(n()) %>%
  mutate(angle = seq(pi / 2, by = -2 * pi / nrow(.), length.out = nrow(.))) %>% 
  mutate(x = (diameter.1/2) * cos(angle), 
         y = (diameter.1/2) * sin(angle) 
  )





#small:

diameter.2 <- 3.5   
nodes.2 <- node%>%dplyr::filter(Part == "Microbes") %>%
  sample_n(n()) %>%
  mutate(angle = seq(pi / 2, by = -2 * pi / nrow(.), length.out = nrow(.))) %>%  
  mutate(x = (diameter.2/2) * cos(angle) + 7, 
         y = (diameter.2/2) * sin(angle) + 2.25    
  )



data.node.plot <- rbind(nodes.1, nodes.2)[-1]

edges <-case_mRNA
colnames(edges)<-c("source","target","value")
edges$target<-gsub("g__","",edges$target)



edges <- edges %>%
  mutate(From.x = data.node.plot$x[match(edges$source, data.node.plot$node)],
         From.y = data.node.plot$y[match(edges$source, data.node.plot$node)],
         To.x = data.node.plot$x[match(edges$target, data.node.plot$node)],
         To.y = data.node.plot$y[match(edges$target, data.node.plot$node)])

edges$correlation<-ifelse(edges$value>0,"Positive","Negative")
edges$correlation<-factor(edges$correlation,levels = c("Positive","Negative"))


load("/mnt/data/IBD/RNA_seq/data/Input/matrix/micro/taxonomy.Rdata")
taxonomy$genus<-str_extract(taxonomy$taxonomy, "g__.+")
taxonomy$phylum<- str_extract(taxonomy$taxonomy, "p__.+")
taxonomy$phylum<-str_split(taxonomy$phylum,"[|]",simplify = T)[,1]
taxonomy$genus<-gsub("g__","",taxonomy$genus)
data.node.plot$phylum<-taxonomy$phylum[match(data.node.plot$node,taxonomy$genus)]


DEGs<-read.csv("/mnt/data/CD/RNA_seq/data/Output/Others/13/DEG_mRNA_Deseq2_count.csv")
data.node.plot$phylum<-ifelse(data.node.plot$node%in%DEGs$X,DEGs$group,data.node.plot$phylum)


p1=ggplot(data = data.node.plot) +
  geom_segment(data = edges,
               aes(xend= To.x, yend = To.y, x = From.x, y = From.y, color = factor(correlation), linewidth = factor(correlation)),
               show.legend = F)+
  scale_linewidth_discrete(range = c(0.07, 0.25))+
  geom_point(aes(x = x, y = y, fill = phylum),
             show.legend = F, size = 7, shape = 21, stroke = 0)+
  coord_fixed()+
  scale_color_manual(values = c("Negative" = "#828282", "Positive" = "#ff7373")) +
  scale_fill_manual(values = c("Up"="#0096DB",
                               "Down"="#0096DB",
                               "None"="#0096DB",
                                "p__Acidobacteria" = "#FF8907",
                               "p__Actinobacteria" = "#FF8907", 
                               "p__Bacteroidota" = "#FF8907", 
                               "p__Chlorobi" = "#FF8907", 
                               "p__Spirochaetes" = "#FF8907", 
                               "p__Firmicutes" = "#FF8907",
                               "p__Fusobacteria" = "#FF8907",
                               "p__Proteobacteria"="#FF8907",
                               "p__Synergistetes"="#FF8907",
                               "p__Thermotogae"="#FF8907",
                               "p__Verrucomicrobia"="#FF8907"))+  
  geom_text(data = data.frame(part = c("", ""),
                              x = c(0,7),
                              y = c(0,2.25)),
            aes(label = part, x = x, y = y), size = 4.7)+ 
  geom_text(data = data.frame(
    label = c(
      paste("Nodes:", nrow(data.node.plot)),
      paste("Edges:", nrow(edges)),
      paste("Positive edges:", sum(edges$correlation == "Positive"))
    ),
    x = 5,
    y = c(-3,-3.5,-4)),
    aes(x = x, y = y, label = label),hjust = 0, size = 4.1)+
  theme_void()











