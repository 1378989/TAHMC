####################################Step1:
#--------------------------------
#############
#---------------microbe:-------------------

df_core_features<-read_tsv(file = "/mnt/data/CD/RNA_seq/data/Output/Others/12/Maaslin2/all_results.tsv")
df_core_features<-df_core_features[df_core_features$qval<0.05,]
core_features<-df_core_features$feature


library(tidyverse)
library(glmnet)
load("/mnt/data/CD/RNA_seq/data/Input/matrix/micro/merge/micro_norm_filter_genus_Bac.Rdata")
metadata<-read.csv("/mnt/data/CD/RNA_seq/data/Supplement/metadata_new.csv")
micro<-micro_norm_filter_genus_Bac
rm(micro_norm_filter_genus_Bac)
identical(colnames(micro),metadata$sample)
meta<-data.frame(sample=colnames(micro))
meta$dataset<-metadata$source[match(meta$sample,metadata$sample)]
meta$group<-metadata$group[match(meta$sample,metadata$sample)]
meta$gender<-metadata$gender[match(meta$sample,metadata$sample)]
rownames(meta)<-meta$sample
identical(rownames(meta),colnames(micro))
meta$group<-factor(meta$group,levels = c("Control","CD"))
micro<-micro[core_features,]

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

data<-expr_combat
data<-t(data)%>%as.data.frame%>%rownames_to_column("sample")
data$gender<-meta$gender[match(data$sample,meta$sample)]
data<-data[,c(1,ncol(data),2:c(ncol(data)-1))]
data$gender<-factor(data$gender,levels = c("female","male"))
#data$group<-ifelse(data$group%in%"Control",0,1)
#data<-data[data$sample%in%train_id,]
data<-data%>%remove_rownames%>%column_to_rownames("sample")
table(data$gender)
microbes_all<-data
microbes_all$gender<-ifelse(microbes_all$gender%in%"female",1,0)
save(microbes_all,file = "/mnt/data/IBD/RNA_seq/data/Output/Others/14/microbes_all.Rdata")


#---------------mRNA:-------------------
#--
library(tidyverse)
library(glmnet)
load("/mnt/data/IBD/RNA_seq/data/Input/matrix/gene/merge/norm_mRNA.Rdata")
metadata<-read.csv("/mnt/data/IBD/RNA_seq/data/Supplement/metadata_new.csv")

identical(colnames(exprset),metadata$sample)
meta<-data.frame(sample=colnames(exprset))
meta$dataset<-metadata$source[match(meta$sample,metadata$sample)]
meta$group<-metadata$group[match(meta$sample,metadata$sample)]
rownames(meta)<-meta$sample
identical(rownames(meta),colnames(exprset))
meta$group<-factor(meta$group,levels = c("Control","CD"))
#
DEG<-read.csv("/mnt/data/IBD/RNA_seq/data/Output/Others/09/DEG_mRNA_Deseq2_count.csv")
DEG_genes<-DEG[DEG$group%in%c("Up","Down"),]$X
exprset<-exprset[DEG_genes,]

library(tinyarray)
#####before：
draw_pca(exp = exprset, group_list = factor(meta$group))
draw_pca(exp =exprset, group_list = factor(meta$dataset))


######after：
library(sva)
mod <- model.matrix(~factor(meta$group))
expr_combat <- ComBat(dat = exprset, batch = meta$dataset,mod = mod)
draw_pca(exp = expr_combat, group_list = factor(meta$group))
draw_pca(exp = expr_combat, group_list = factor(meta$dataset))

data<-expr_combat
data<-t(data)%>%as.data.frame%>%rownames_to_column("sample")
#data<-data[data$sample%in%train_id,]
data<-data%>%remove_rownames%>%column_to_rownames("sample")
table(data$group)
mRNA_all<-data
save(mRNA_all,file = "/mnt/data/IBD/RNA_seq/data/Output/Others/14/mRNA_all.Rdata")




#---------------lncRNA:-------------------
#--
library(tidyverse)
library(glmnet)
load("/mnt/data/IBD/RNA_seq/data/Input/matrix/gene/merge/norm_LncRNA.Rdata")
metadata<-read.csv("/mnt/data/IBD/RNA_seq/data/Supplement/metadata_new.csv")

identical(colnames(exprset),metadata$sample)
meta<-data.frame(sample=colnames(exprset))
meta$dataset<-metadata$source[match(meta$sample,metadata$sample)]
meta$group<-metadata$group[match(meta$sample,metadata$sample)]
rownames(meta)<-meta$sample
identical(rownames(meta),colnames(exprset))
meta$group<-factor(meta$group,levels = c("Control","CD"))
#
DEG<-read.csv("/mnt/data/IBD/RNA_seq/data/Output/Others/09/DEG_lncRNA_Deseq2_count.csv")
DEG_genes<-DEG[DEG$group%in%c("Up","Down"),]$X
exprset<-exprset[DEG_genes,]

library(tinyarray)
#####before：
draw_pca(exp = exprset, group_list = factor(meta$group))
draw_pca(exp =exprset, group_list = factor(meta$dataset))


######after：
library(sva)
mod <- model.matrix(~factor(meta$group))
expr_combat <- ComBat(dat = exprset, batch = meta$dataset,mod = mod)
draw_pca(exp = expr_combat, group_list = factor(meta$group))
draw_pca(exp = expr_combat, group_list = factor(meta$dataset))

data<-expr_combat
data<-t(data)%>%as.data.frame%>%rownames_to_column("sample")
#data<-data[data$sample%in%train_id,]
data<-data%>%remove_rownames%>%column_to_rownames("sample")
table(data$group)
LncRNA_all<-data
save(LncRNA_all,file = "/mnt/data/IBD/RNA_seq/data/Output/Others/14/LncRNA_all.Rdata")


##################################################################################
#

Contr<-meta[meta$group%in%"Control",]$sample
CD<-meta[meta$group%in%"CD",]$sample



microbes_Contr<-microbes_all[rownames(microbes_all)%in%Contr,]
microbes_CD<-microbes_all[rownames(microbes_all)%in%CD,]


mRNA_Contr<-mRNA_all[rownames(mRNA_all)%in%Contr,]
mRNA_CD<-mRNA_all[rownames(mRNA_all)%in%CD,]


LncRNA_Contr<-LncRNA_all[rownames(LncRNA_all)%in%Contr,]
LncRNA_CD<-LncRNA_all[rownames(LncRNA_all)%in%CD,]




save(microbes_Contr,file = "/mnt/data/IBD/RNA_seq/data/Output/Others/14/microbes_Contr.Rdata")
save(microbes_CD,file = "/mnt/data/IBD/RNA_seq/data/Output/Others/14/microbes_CD.Rdata")
save(mRNA_Contr,file = "/mnt/data/IBD/RNA_seq/data/Output/Others/14/mRNA_Contr.Rdata")
save(mRNA_CD,file = "/mnt/data/IBD/RNA_seq/data/Output/Others/14/mRNA_CD.Rdata")
save(LncRNA_Contr,file = "/mnt/data/IBD/RNA_seq/data/Output/Others/14/LncRNA_Contr.Rdata")
save(LncRNA_CD,file = "/mnt/data/IBD/RNA_seq/data/Output/Others/14/LncRNA_CD.Rdata")



#---------------------------------------
rm(list = ls())
load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/LncRNA_all.Rdata")
load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/mRNA_all.Rdata")


gene_list1<-data.frame(gene=colnames(mRNA_all))
gene_list2<-data.frame(gene=colnames(LncRNA_all))

#write.csv(gene_list1,file = "/mnt/data/IBD/RNA_seq/data/Output/Others/14/gene_list1.csv",row.names = F)
#write.csv(gene_list2,file = "/mnt/data/IBD/RNA_seq/data/Output/Others/14/gene_list2.csv",row.names = F)



#############################################################Step2:


######
#---------------------------------------hdi_lasso_proj_loocv_MSI---------------------------------
################################################################################################

################################################################################################


##################################mRNA————————————CD:######################

rm(list=ls()) ## ## check to make sure no precomputed dataframe before clearing

## Setup for parallel processing on MSI
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("doParallel") ## package methods is not loaded by default by RScript. 
check.packages(packages)

## now initiate parallel processors.
cl <- makeCluster(60) # 8 workers per core -- itasca
registerDoParallel(cl)
print(paste0("Number of workers: ", getDoParWorkers()))



################# Input genes and taxa matrix ###########
input.genes.list="/mnt/data/IBD/RNA_seq/data/Output/Others/14/gene_list1.txt"
output.dir="/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/lasso_hdi_loocv"

load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/mRNA_CD.Rdata")
load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/microbes_CD.Rdata")


genes <-mRNA_CD
microbes <- microbes_CD

## convert both dataframes to matrix. cv.glmnet expects a matrix of predictors, not a data frame
y <- as.matrix(genes) #response
x <- as.matrix(microbes) #predictors

x.uniq <- unique(x, MARGIN = 2)
dim(x.uniq)
x <- x.uniq

print(paste0("# genes = ",dim(y)[2]))

print(paste0("# microbes = ",dim(x)[2]))


############### Functions ############

estimate.sigma <- function(x, y_i, bestlambda, tol) {
  
  lasso.fit = glmnet(x,y_i,alpha = 1)
  # beta <- coef(lasso.fit, s = bestlambda)[-1] ## This gives coefficients of fitted model, not predicted coeff. 
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
    sigma.flag = 0
  } else{
    
    sigmas <- numeric()
    ## Repeat 5x and take median of finite values else assign sigma to ss_res and raise the flag 
    for( j in 1:5){
      bestlambda = get.lambda(x, y_i, 10, 10)
      lasso.fit = glmnet(x,y_i,alpha = 1)
      y = as.vector(y_i)
      yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
      beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
      df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
      n = length(y_i)
      ss_res = sum((y - yhat)^2)
      if((n-df-1) >= 1) {
        sigmas[j] = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
      }
    }
    if(length(sigmas != 0 )){
      
      sigma = median(sigmas, na.rm = T)
      sigma.flag = 1
    } else{
      sigma = 1 ## conservative option
      # sigma = ss_res ## lenient option
      sigma.flag = 2
    }
    
  }
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}

## Inspired from the details for estimateSigma() function in selectiveInference package: https://cran.r-project.org/web/packages/selectiveInference/selectiveInference.pdf
estimate.sigma.loocv <- function(x, y_i, bestlambda, tol) {
  
  lasso.fit = glmnet(x,y_i,alpha = 1)
  # beta <- coef(lasso.fit, s = bestlambda)[-1] ## This gives coefficients of fitted model, not predicted coeff. 
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
    sigma.flag = 0
  } else{
    sigma = 1 ## conservative option
    # sigma = ss_res ## lenient option
    sigma.flag = 2
  }
  
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}


fit.cv.lasso <- function(x, y_i, kfold, repeats){
  
  lambdas = NULL
  r.sqr.final <- numeric()
  r.sqr.final.adj <- numeric()
  r.sqr.CV.test <- numeric()
  lasso.cv.list <- list()
  
  for (i in 1:repeats){
    
    ## glmnet CV
    fit <- cv.glmnet(x, y_i, alpha=1, nfolds=kfold, type.measure = "mse", keep =TRUE, grouped=FALSE)  
    errors = data.frame(fit$lambda,fit$cvm)
    lambdas <- rbind(lambdas,errors)
    
    ## Get R^2 of final model
    r.sqr.final[i] <- r_squared(as.vector(y_i), 
                                as.vector(predict(fit$glmnet.fit, 
                                                  newx = x, s = fit$lambda.min)))
    ## Get adjusted R^2
    r.sqr.final.adj[i] <- adj_r_squared(r.sqr.final[i], n = nrow(x), 
                                        p = sum(as.vector(coef(fit$glmnet.fit, 
                                                               s = fit$lambda.min)) > 0))
    
  }
  
  # take mean cvm for each lambda
  lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)
  # dim(lambdas)
  # select the best one
  bestindex = which(lambdas[2]==min(lambdas[2]))
  bestlambda = lambdas[bestindex,1]
  
  return(list(bestlambda = bestlambda, r.sqr = median(r.sqr.final), 
              r.sqr.adj = median(r.sqr.final.adj)
  ))
}

## functions to compute R2
## Adapted from https://rpubs.com/kaz_yos/alasso
r_squared <- function(y, yhat) {
  ybar <- mean(y)
  ## Total SS
  ss_tot <- sum((y - ybar)^2)
  ## Residual SS
  ss_res <- sum((y - yhat)^2)
  ## R^2 = 1 - ss_res/ ss_tot
  1 - (ss_res / ss_tot)
}
## Function for Adjusted R^2
## n sample size, p number of prameters
adj_r_squared <- function(r_squared, n, p) {
  1 - (1 - r_squared) * (n - 1) / (n - p - 1)
}

############### Estimate sigma (standard deviation of the error term or noise) ###############

estimate.sigma.fit.hdi <- function(x, y, gene_name){
  ## Import all the libraries for the current node/core on MSI 
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("glmnet","hdi","methods","doParallel") ## package methods is not loaded by default by RScript on MSI 
  check.packages(packages)
  
  print(paste0("Processing gene:", gene_name));flush.console()
  
  ## Extract the expression for this gene (response variable)
  y_i <- y[,grep(paste0("^",gene_name,"$"),colnames(y))]
  
  ## Make sure y_i is numeric before model fitting 
  stopifnot(class(y_i) == "numeric")
  
  ## Fit lasso CV model
  fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i), repeats = 1)
  bestlambda <- fit.model$bestlambda
  r.sqr <- fit.model$r.sqr
  r.sqr.adj <- fit.model$r.sqr.adj
  
  ## Estimate sigma using the estimated lambda param
  sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
  sigma <- sigma.myfun$sigmahat
  beta <- as.vector(sigma.myfun$betahat)[-1]
  sigma.flag <- sigma.myfun$sigmaflag
  
  ## Inference 
  lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
  lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))
  lasso.FDR.df <- data.frame(gene = rep(gene_name, length(lasso.proj.fit$pval)), 
                             taxa = names(lasso.proj.fit$pval.corr), 
                             r.sqr = r.sqr, r.sqr.adj = r.sqr.adj,
                             pval = lasso.proj.fit$pval, padj = lasso.proj.fit$pval.corr, 
                             ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
                             sigma = sigma, sigma.flag = sigma.flag,
                             row.names=NULL)
  
  
  return(lasso.FDR.df)
  
}


## Parallel run on MSI
gene_list <- scan(input.genes.list, what="", sep="\n")



## invoke parallel computation for fitting model for each gene
suppressWarnings(parallel_time <- system.time({
  parallel_res <- foreach(i=1:length(gene_list), .combine = rbind) %dopar% estimate.sigma.fit.hdi(x,y,gene_list[i])
}))

stopCluster(cl)
registerDoSEQ()

## time taken for parallel computation.
print(parallel_time)

## Print result for this node's genes list
filename <- strsplit(input.genes.list,"/")
filename <- strsplit(filename[[1]][length(filename[[1]])],"\\.")[[1]][1]
write.table(parallel_res, file=paste0(output.dir,"/",filename,"_lasso_hdi.txt"), quote=F, sep="\t",col.names = NA, row.names = TRUE)


#-------------------------------------Control ----mRNA：------------------------------

##################################mRNA————————————Control:######################

rm(list=ls()) ## ## check to make sure no precomputed dataframe before clearing

## Setup for parallel processing on MSI
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("doParallel") ## package methods is not loaded by default by RScript. 
check.packages(packages)

## now initiate parallel processors.
cl <- makeCluster(60) # 8 workers per core -- itasca
registerDoParallel(cl)
print(paste0("Number of workers: ", getDoParWorkers()))



################# Input genes and taxa matrix ###########
input.genes.list="/mnt/data/IBD/RNA_seq/data/Output/Others/14/gene_list1.txt"
output.dir="/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/lasso_hdi_loocv"

load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/mRNA_Contr.Rdata")
load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/microbes_Contr.Rdata")


genes <-mRNA_Contr
microbes <- microbes_Contr

## convert both dataframes to matrix. cv.glmnet expects a matrix of predictors, not a data frame
y <- as.matrix(genes) #response
x <- as.matrix(microbes) #predictors

x.uniq <- unique(x, MARGIN = 2)
dim(x.uniq)
x <- x.uniq

print(paste0("# genes = ",dim(y)[2]))

print(paste0("# microbes = ",dim(x)[2]))


#colSums(x!=0)%>%sort()
#colSums(x!=0)>0.05*nrow(x)

x=x[,colSums(x!=0)>0.05*nrow(x)]
dim(x)


############### Functions ############

estimate.sigma <- function(x, y_i, bestlambda, tol) {
  
  lasso.fit = glmnet(x,y_i,alpha = 1)
  # beta <- coef(lasso.fit, s = bestlambda)[-1] ## This gives coefficients of fitted model, not predicted coeff. 
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
    sigma.flag = 0
  } else{
    
    sigmas <- numeric()
    ## Repeat 5x and take median of finite values else assign sigma to ss_res and raise the flag 
    for( j in 1:5){
      bestlambda = get.lambda(x, y_i, 10, 10)
      lasso.fit = glmnet(x,y_i,alpha = 1)
      y = as.vector(y_i)
      yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
      beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
      df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
      n = length(y_i)
      ss_res = sum((y - yhat)^2)
      if((n-df-1) >= 1) {
        sigmas[j] = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
      }
    }
    if(length(sigmas != 0 )){
      
      sigma = median(sigmas, na.rm = T)
      sigma.flag = 1
    } else{
      sigma = 1 ## conservative option
      # sigma = ss_res ## lenient option
      sigma.flag = 2
    }
    
  }
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}

## Inspired from the details for estimateSigma() function in selectiveInference package: https://cran.r-project.org/web/packages/selectiveInference/selectiveInference.pdf
estimate.sigma.loocv <- function(x, y_i, bestlambda, tol) {
  
  lasso.fit = glmnet(x,y_i,alpha = 1)
  # beta <- coef(lasso.fit, s = bestlambda)[-1] ## This gives coefficients of fitted model, not predicted coeff. 
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
    sigma.flag = 0
  } else{
    sigma = 1 ## conservative option
    # sigma = ss_res ## lenient option
    sigma.flag = 2
  }
  
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}


fit.cv.lasso <- function(x, y_i, kfold, repeats){
  
  lambdas = NULL
  r.sqr.final <- numeric()
  r.sqr.final.adj <- numeric()
  r.sqr.CV.test <- numeric()
  lasso.cv.list <- list()
  
  for (i in 1:repeats){
    
    ## glmnet CV
    fit <- cv.glmnet(x, y_i, alpha=1, nfolds=kfold, type.measure = "mse", keep =TRUE, grouped=FALSE)  
    errors = data.frame(fit$lambda,fit$cvm)
    lambdas <- rbind(lambdas,errors)
    
    ## Get R^2 of final model
    r.sqr.final[i] <- r_squared(as.vector(y_i), 
                                as.vector(predict(fit$glmnet.fit, 
                                                  newx = x, s = fit$lambda.min)))
    ## Get adjusted R^2
    r.sqr.final.adj[i] <- adj_r_squared(r.sqr.final[i], n = nrow(x), 
                                        p = sum(as.vector(coef(fit$glmnet.fit, 
                                                               s = fit$lambda.min)) > 0))
    
  }
  
  # take mean cvm for each lambda
  lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)
  # dim(lambdas)
  # select the best one
  bestindex = which(lambdas[2]==min(lambdas[2]))
  bestlambda = lambdas[bestindex,1]
  
  return(list(bestlambda = bestlambda, r.sqr = median(r.sqr.final), 
              r.sqr.adj = median(r.sqr.final.adj)
  ))
}

## functions to compute R2
## Adapted from https://rpubs.com/kaz_yos/alasso
r_squared <- function(y, yhat) {
  ybar <- mean(y)
  ## Total SS
  ss_tot <- sum((y - ybar)^2)
  ## Residual SS
  ss_res <- sum((y - yhat)^2)
  ## R^2 = 1 - ss_res/ ss_tot
  1 - (ss_res / ss_tot)
}
## Function for Adjusted R^2
## n sample size, p number of prameters
adj_r_squared <- function(r_squared, n, p) {
  1 - (1 - r_squared) * (n - 1) / (n - p - 1)
}

############### Estimate sigma (standard deviation of the error term or noise) ###############

estimate.sigma.fit.hdi <- function(x, y, gene_name){
  ## Import all the libraries for the current node/core on MSI 
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("glmnet","hdi","methods","doParallel") ## package methods is not loaded by default by RScript on MSI 
  check.packages(packages)
  
  print(paste0("Processing gene:", gene_name));flush.console()
  
  ## Extract the expression for this gene (response variable)
  y_i <- y[,grep(paste0("^",gene_name,"$"),colnames(y))]
  
  ## Make sure y_i is numeric before model fitting 
  stopifnot(class(y_i) == "numeric")
  
  ## Fit lasso CV model
  fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i), repeats = 1)
  bestlambda <- fit.model$bestlambda
  r.sqr <- fit.model$r.sqr
  r.sqr.adj <- fit.model$r.sqr.adj
  
  ## Estimate sigma using the estimated lambda param
  sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
  sigma <- sigma.myfun$sigmahat
  beta <- as.vector(sigma.myfun$betahat)[-1]
  sigma.flag <- sigma.myfun$sigmaflag
  
  ## Inference 
  lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
  lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))
  lasso.FDR.df <- data.frame(gene = rep(gene_name, length(lasso.proj.fit$pval)), 
                             taxa = names(lasso.proj.fit$pval.corr), 
                             r.sqr = r.sqr, r.sqr.adj = r.sqr.adj,
                             pval = lasso.proj.fit$pval, padj = lasso.proj.fit$pval.corr, 
                             ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
                             sigma = sigma, sigma.flag = sigma.flag,
                             row.names=NULL)
  
  
  return(lasso.FDR.df)
  
}


## Parallel run on MSI
gene_list <- scan(input.genes.list, what="", sep="\n")



## invoke parallel computation for fitting model for each gene
suppressWarnings(parallel_time <- system.time({
  parallel_res <- foreach(i=1:length(gene_list), .combine = rbind) %dopar% estimate.sigma.fit.hdi(x,y,gene_list[i])
}))


stopCluster(cl)
registerDoSEQ()

## time taken for parallel computation.
print(parallel_time)

## Print result for this node's genes list
filename <- strsplit(input.genes.list,"/")
filename <- strsplit(filename[[1]][length(filename[[1]])],"\\.")[[1]][1]
write.table(parallel_res, file=paste0(output.dir,"/",filename,"_lasso_hdi.txt"), quote=F, sep="\t",col.names = NA, row.names = TRUE)




#############################################################################################################################################
#----------------------------------------------------------------------------------------------
#-------------------------------------CD ----LncRNA：------------------------------

##################################LncRNA————————————CD:######################

rm(list=ls()) ## ## check to make sure no precomputed dataframe before clearing

## Setup for parallel processing on MSI
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("doParallel") ## package methods is not loaded by default by RScript. 
check.packages(packages)

## now initiate parallel processors.
cl <- makeCluster(60) # 8 workers per core -- itasca




registerDoParallel(cl)
print(paste0("Number of workers: ", getDoParWorkers()))



################# Input genes and taxa matrix ###########
input.genes.list="/mnt/data/IBD/RNA_seq/data/Output/Others/14/gene_list2.txt"
output.dir="/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/lasso_hdi_loocv"

load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/LncRNA_CD.Rdata")
load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/microbes_CD.Rdata")


genes <-LncRNA_CD
microbes <- microbes_CD

## convert both dataframes to matrix. cv.glmnet expects a matrix of predictors, not a data frame
y <- as.matrix(genes) #response
x <- as.matrix(microbes) #predictors

x.uniq <- unique(x, MARGIN = 2)
dim(x.uniq)
x <- x.uniq

print(paste0("# genes = ",dim(y)[2]))

print(paste0("# microbes = ",dim(x)[2]))


############### Functions ############

estimate.sigma <- function(x, y_i, bestlambda, tol) {
  
  lasso.fit = glmnet(x,y_i,alpha = 1)
  # beta <- coef(lasso.fit, s = bestlambda)[-1] ## This gives coefficients of fitted model, not predicted coeff. 
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
    sigma.flag = 0
  } else{
    
    sigmas <- numeric()
    ## Repeat 5x and take median of finite values else assign sigma to ss_res and raise the flag 
    for( j in 1:5){
      bestlambda = get.lambda(x, y_i, 10, 10)
      lasso.fit = glmnet(x,y_i,alpha = 1)
      y = as.vector(y_i)
      yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
      beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
      df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
      n = length(y_i)
      ss_res = sum((y - yhat)^2)
      if((n-df-1) >= 1) {
        sigmas[j] = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
      }
    }
    if(length(sigmas != 0 )){
      
      sigma = median(sigmas, na.rm = T)
      sigma.flag = 1
    } else{
      sigma = 1 ## conservative option
      # sigma = ss_res ## lenient option
      sigma.flag = 2
    }
    
  }
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}

## Inspired from the details for estimateSigma() function in selectiveInference package: https://cran.r-project.org/web/packages/selectiveInference/selectiveInference.pdf
estimate.sigma.loocv <- function(x, y_i, bestlambda, tol) {
  
  lasso.fit = glmnet(x,y_i,alpha = 1)
  # beta <- coef(lasso.fit, s = bestlambda)[-1] ## This gives coefficients of fitted model, not predicted coeff. 
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
    sigma.flag = 0
  } else{
    sigma = 1 ## conservative option
    # sigma = ss_res ## lenient option
    sigma.flag = 2
  }
  
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}


fit.cv.lasso <- function(x, y_i, kfold, repeats){
  
  lambdas = NULL
  r.sqr.final <- numeric()
  r.sqr.final.adj <- numeric()
  r.sqr.CV.test <- numeric()
  lasso.cv.list <- list()
  
  for (i in 1:repeats){
    
    ## glmnet CV
    fit <- cv.glmnet(x, y_i, alpha=1, nfolds=kfold, type.measure = "mse", keep =TRUE, grouped=FALSE)  
    errors = data.frame(fit$lambda,fit$cvm)
    lambdas <- rbind(lambdas,errors)
    
    ## Get R^2 of final model
    r.sqr.final[i] <- r_squared(as.vector(y_i), 
                                as.vector(predict(fit$glmnet.fit, 
                                                  newx = x, s = fit$lambda.min)))
    ## Get adjusted R^2
    r.sqr.final.adj[i] <- adj_r_squared(r.sqr.final[i], n = nrow(x), 
                                        p = sum(as.vector(coef(fit$glmnet.fit, 
                                                               s = fit$lambda.min)) > 0))
    
  }
  
  # take mean cvm for each lambda
  lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)
  # dim(lambdas)
  # select the best one
  bestindex = which(lambdas[2]==min(lambdas[2]))
  bestlambda = lambdas[bestindex,1]
  
  return(list(bestlambda = bestlambda, r.sqr = median(r.sqr.final), 
              r.sqr.adj = median(r.sqr.final.adj)
  ))
}

## functions to compute R2
## Adapted from https://rpubs.com/kaz_yos/alasso
r_squared <- function(y, yhat) {
  ybar <- mean(y)
  ## Total SS
  ss_tot <- sum((y - ybar)^2)
  ## Residual SS
  ss_res <- sum((y - yhat)^2)
  ## R^2 = 1 - ss_res/ ss_tot
  1 - (ss_res / ss_tot)
}
## Function for Adjusted R^2
## n sample size, p number of prameters
adj_r_squared <- function(r_squared, n, p) {
  1 - (1 - r_squared) * (n - 1) / (n - p - 1)
}

############### Estimate sigma (standard deviation of the error term or noise) ###############

estimate.sigma.fit.hdi <- function(x, y, gene_name){
  ## Import all the libraries for the current node/core on MSI 
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("glmnet","hdi","methods","doParallel") ## package methods is not loaded by default by RScript on MSI 
  check.packages(packages)
  
  print(paste0("Processing gene:", gene_name));flush.console()
  
  ## Extract the expression for this gene (response variable)
  y_i <- y[,grep(paste0("^",gene_name,"$"),colnames(y))]
  
  ## Make sure y_i is numeric before model fitting 
  stopifnot(class(y_i) == "numeric")
  
  ## Fit lasso CV model
  fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i), repeats = 1)
  bestlambda <- fit.model$bestlambda
  r.sqr <- fit.model$r.sqr
  r.sqr.adj <- fit.model$r.sqr.adj
  
  ## Estimate sigma using the estimated lambda param
  sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
  sigma <- sigma.myfun$sigmahat
  beta <- as.vector(sigma.myfun$betahat)[-1]
  sigma.flag <- sigma.myfun$sigmaflag
  
  ## Inference 
  lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
  lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))
  lasso.FDR.df <- data.frame(gene = rep(gene_name, length(lasso.proj.fit$pval)), 
                             taxa = names(lasso.proj.fit$pval.corr), 
                             r.sqr = r.sqr, r.sqr.adj = r.sqr.adj,
                             pval = lasso.proj.fit$pval, padj = lasso.proj.fit$pval.corr, 
                             ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
                             sigma = sigma, sigma.flag = sigma.flag,
                             row.names=NULL)
  
  
  return(lasso.FDR.df)
  
}


## Parallel run on MSI
gene_list <- scan(input.genes.list, what="", sep="\n")



## invoke parallel computation for fitting model for each gene
suppressWarnings(parallel_time <- system.time({
  parallel_res <- foreach(i=1:length(gene_list), .combine = rbind) %dopar% estimate.sigma.fit.hdi(x,y,gene_list[i])
}))


stopCluster(cl)
registerDoSEQ()

## time taken for parallel computation.
print(parallel_time)

## Print result for this node's genes list
filename <- strsplit(input.genes.list,"/")
filename <- strsplit(filename[[1]][length(filename[[1]])],"\\.")[[1]][1]
write.table(parallel_res, file=paste0(output.dir,"/",filename,"_lasso_hdi.txt"), quote=F, sep="\t",col.names = NA, row.names = TRUE)








#############################################################################################################################################
#----------------------------------------------------------------------------------------------
#-------------------------------------Control----LncRNA：------------------------------

##################################LncRNA————————————Control:######################

rm(list=ls()) ## ## check to make sure no precomputed dataframe before clearing

## Setup for parallel processing on MSI
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("doParallel") ## package methods is not loaded by default by RScript. 
check.packages(packages)

## now initiate parallel processors.
cl <- makeCluster(60) # 8 workers per core -- itasca
registerDoParallel(cl)
print(paste0("Number of workers: ", getDoParWorkers()))



################# Input genes and taxa matrix ###########
input.genes.list="/mnt/data/IBD/RNA_seq/data/Output/Others/14/gene_list2.txt"
output.dir="/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/lasso_hdi_loocv"

load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/LncRNA_Contr.Rdata")
load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/microbes_Contr.Rdata")


genes <-LncRNA_Contr
microbes <- microbes_Contr



## convert both dataframes to matrix. cv.glmnet expects a matrix of predictors, not a data frame
y <- as.matrix(genes) #response
x <- as.matrix(microbes) #predictors

x.uniq <- unique(x, MARGIN = 2)
dim(x.uniq)
x <- x.uniq

colSums(x!=0)%>%sort()

print(paste0("# genes = ",dim(y)[2]))

print(paste0("# microbes = ",dim(x)[2]))



#colSums(x!=0)>0.05*nrow(x)

#x=x[,colSums(x!=0)>0.05*nrow(x)]
#dim(x)
############### Functions ############

estimate.sigma <- function(x, y_i, bestlambda, tol) {
  
  lasso.fit = glmnet(x,y_i,alpha = 1)
  # beta <- coef(lasso.fit, s = bestlambda)[-1] ## This gives coefficients of fitted model, not predicted coeff. 
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
    sigma.flag = 0
  } else{
    
    sigmas <- numeric()
    ## Repeat 5x and take median of finite values else assign sigma to ss_res and raise the flag 
    for( j in 1:5){
      bestlambda = get.lambda(x, y_i, 10, 10)
      lasso.fit = glmnet(x,y_i,alpha = 1)
      y = as.vector(y_i)
      yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
      beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
      df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
      n = length(y_i)
      ss_res = sum((y - yhat)^2)
      if((n-df-1) >= 1) {
        sigmas[j] = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
      }
    }
    if(length(sigmas != 0 )){
      
      sigma = median(sigmas, na.rm = T)
      sigma.flag = 1
    } else{
      sigma = 1 ## conservative option
      # sigma = ss_res ## lenient option
      sigma.flag = 2
    }
    
  }
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}

## Inspired from the details for estimateSigma() function in selectiveInference package: https://cran.r-project.org/web/packages/selectiveInference/selectiveInference.pdf
estimate.sigma.loocv <- function(x, y_i, bestlambda, tol) {
  
  lasso.fit = glmnet(x,y_i,alpha = 1)
  # beta <- coef(lasso.fit, s = bestlambda)[-1] ## This gives coefficients of fitted model, not predicted coeff. 
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  beta = predict(lasso.fit,s=bestlambda, type="coef") ## predicted coefficients
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1)) ## note, if we get rid of intercept when extracting beta, so should be drop -1?
    sigma.flag = 0
  } else{
    sigma = 1 ## conservative option
    # sigma = ss_res ## lenient option
    sigma.flag = 2
  }
  
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}


fit.cv.lasso <- function(x, y_i, kfold, repeats){
  
  lambdas = NULL
  r.sqr.final <- numeric()
  r.sqr.final.adj <- numeric()
  r.sqr.CV.test <- numeric()
  lasso.cv.list <- list()
  
  for (i in 1:repeats){
    
    ## glmnet CV
    fit <- cv.glmnet(x, y_i, alpha=1, nfolds=kfold, type.measure = "mse", keep =TRUE, grouped=FALSE)  
    errors = data.frame(fit$lambda,fit$cvm)
    lambdas <- rbind(lambdas,errors)
    
    ## Get R^2 of final model
    r.sqr.final[i] <- r_squared(as.vector(y_i), 
                                as.vector(predict(fit$glmnet.fit, 
                                                  newx = x, s = fit$lambda.min)))
    ## Get adjusted R^2
    r.sqr.final.adj[i] <- adj_r_squared(r.sqr.final[i], n = nrow(x), 
                                        p = sum(as.vector(coef(fit$glmnet.fit, 
                                                               s = fit$lambda.min)) > 0))
    
  }
  
  # take mean cvm for each lambda
  lambdas <- aggregate(lambdas[, 2], list(lambdas$fit.lambda), mean)
  # dim(lambdas)
  # select the best one
  bestindex = which(lambdas[2]==min(lambdas[2]))
  bestlambda = lambdas[bestindex,1]
  
  return(list(bestlambda = bestlambda, r.sqr = median(r.sqr.final), 
              r.sqr.adj = median(r.sqr.final.adj)
  ))
}

## functions to compute R2
## Adapted from https://rpubs.com/kaz_yos/alasso
r_squared <- function(y, yhat) {
  ybar <- mean(y)
  ## Total SS
  ss_tot <- sum((y - ybar)^2)
  ## Residual SS
  ss_res <- sum((y - yhat)^2)
  ## R^2 = 1 - ss_res/ ss_tot
  1 - (ss_res / ss_tot)
}
## Function for Adjusted R^2
## n sample size, p number of prameters
adj_r_squared <- function(r_squared, n, p) {
  1 - (1 - r_squared) * (n - 1) / (n - p - 1)
}

############### Estimate sigma (standard deviation of the error term or noise) ###############

estimate.sigma.fit.hdi <- function(x, y, gene_name){
  ## Import all the libraries for the current node/core on MSI 
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("glmnet","hdi","methods","doParallel") ## package methods is not loaded by default by RScript on MSI 
  check.packages(packages)
  
  print(paste0("Processing gene:", gene_name));flush.console()
  
  ## Extract the expression for this gene (response variable)
  y_i <- y[,grep(paste0("^",gene_name,"$"),colnames(y))]
  
  ## Make sure y_i is numeric before model fitting 
  stopifnot(class(y_i) == "numeric")
  
  ## Fit lasso CV model
  fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i), repeats = 1)
  bestlambda <- fit.model$bestlambda
  r.sqr <- fit.model$r.sqr
  r.sqr.adj <- fit.model$r.sqr.adj
  
  ## Estimate sigma using the estimated lambda param
  sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
  sigma <- sigma.myfun$sigmahat
  beta <- as.vector(sigma.myfun$betahat)[-1]
  sigma.flag <- sigma.myfun$sigmaflag
  
  ## Inference 
  lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
  lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))
  lasso.FDR.df <- data.frame(gene = rep(gene_name, length(lasso.proj.fit$pval)), 
                             taxa = names(lasso.proj.fit$pval.corr), 
                             r.sqr = r.sqr, r.sqr.adj = r.sqr.adj,
                             pval = lasso.proj.fit$pval, padj = lasso.proj.fit$pval.corr, 
                             ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
                             sigma = sigma, sigma.flag = sigma.flag,
                             row.names=NULL)
  
  
  return(lasso.FDR.df)
  
}


## Parallel run on MSI
gene_list <- scan(input.genes.list, what="", sep="\n")



## invoke parallel computation for fitting model for each gene
suppressWarnings(parallel_time <- system.time({
  parallel_res <- foreach(i=1:length(gene_list), .combine = rbind) %dopar% estimate.sigma.fit.hdi(x,y,gene_list[i])
}))


stopCluster(cl)
registerDoSEQ()

## time taken for parallel computation.
print(parallel_time)

## Print result for this node's genes list
filename <- strsplit(input.genes.list,"/")
filename <- strsplit(filename[[1]][length(filename[[1]])],"\\.")[[1]][1]
write.table(parallel_res, file=paste0(output.dir,"/",filename,"_lasso_hdi.txt"), quote=F, sep="\t",col.names = NA, row.names = TRUE)




### This script performs stability selection for lasso model to identify
## robustly associated microbes with host genes.
## Here we use the "stabs" R package's implementation of 
## stability selection.
## Stabs paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0575-3
## Stabs manual: https://cran.r-project.org/web/packages/stabs/stabs.pdf 

## This code is run in parallel on MSI supercomputing nodes.


#-----------------------stabs_stability_selection-----------------------------

#--------------------mRNA----------CD --------------------


######################### Initial setup ######################
## before initiating parallel 
rm(list=ls()) ## clear workspace

## Setup for parallel processing on MSI
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("doParallel", "data.table") ## package methods is not loaded by default by RScript on MSI
check.packages(packages)

## now initiate parallel processors.
cl <- makeCluster(60) # 8 workers per core -- itasca
registerDoParallel(cl)
print(paste0("Number of workers: ", getDoParWorkers()))



input.genes.list="/mnt/data/IBD/RNA_seq/data/Output/Others/14/gene_list1.txt"
output.dir="/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/stabsel_path"

load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/mRNA_CD.Rdata")
load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/microbes_CD.Rdata")


genes <-mRNA_CD
microbes <- microbes_CD

print(paste0("# genes = ",dim(genes)[2]))

print(paste0("# microbes = ",dim(microbes)[2]))

## Ensure same sampleIDs in both genes and microbes matrices
stopifnot(all(rownames(genes) == rownames(microbes)))

## convert both dataframes to matrix. cv.glmnet expects a matrix of predictors, not a data frame
y <- as.matrix(genes) #response
x <- as.matrix(microbes) #predictors


##################### Stability selection ####################

stabs_stability <- function( x, y, gene_name){
  
  ## Import all the libraries for the current node/core on MSI 
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("stabs","glmnet","methods","doParallel") ## package methods is not loaded by default by RScript. 
  check.packages(packages)
  
  ## Extract the expression for this gene (response variable)
  y_i <- y[,grep(paste0("^",gene_name,"$"),colnames(y))]
  
  ## Make sure y_i is numeric before model fitting 
  stopifnot(class(y_i) == "numeric")
  
  ## perform stability selection using glmnet lasso
  stab.glmnet <- stabsel(x = x, y = y_i,
                         fitfun = glmnet.lasso, cutoff = 0.6,
                         PFER = 1)
  
  taxa.selected <- names(stab.glmnet$selected)
  if(length(taxa.selected) == 0) taxa.selected <-"None"
  
  taxa.selected.df <- as.data.frame(cbind(gene_name,taxa.selected))
  
  return(taxa.selected.df)
  
}

## Parallel run on MSI
gene_list <- scan(input.genes.list, what="", sep="\n")

## invoke parallel computation for fitting model for each gene
parallel_time <- system.time({
  # sink(paste0(output.dir,"/log.txt"), append=TRUE)
  parallel_res <- foreach(i=1:length(gene_list), .combine = rbind) %dopar% {
    # cat(paste("Processing gene:",gene_list[i],"\n"))
    # sink() #end diversion of output
    stabs_stability(x,y,gene_list[i])
  }
})
stopCluster(cl)
registerDoSEQ()

## time taken for parallel computation.
print(parallel_time)

filename <- strsplit(input.genes.list,"/")
filename <- strsplit(filename[[1]][length(filename[[1]])],"\\.")[[1]][1]
write.table(parallel_res, file=paste0(output.dir,"/",filename,"_stabs_stabsel_output.txt"), quote=F, sep="\t",col.names = NA, row.names = TRUE)



#--------------------mRNA----------Control --------------------


######################### Initial setup ######################
## before initiating parallel 
rm(list=ls()) ## clear workspace

## Setup for parallel processing on MSI
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("doParallel", "data.table") ## package methods is not loaded by default by RScript on MSI
check.packages(packages)

## now initiate parallel processors.
cl <- makeCluster(60) # 8 workers per core -- itasca
registerDoParallel(cl)
print(paste0("Number of workers: ", getDoParWorkers()))



################# Input genes and taxa matrix ###########
input.genes.list="/mnt/data/IBD/RNA_seq/data/Output/Others/14/gene_list1.txt"
output.dir="/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/stabsel_path"

load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/mRNA_Contr.Rdata")
load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/microbes_Contr.Rdata")


genes <-mRNA_Contr
microbes <- microbes_Contr

print(paste0("# genes = ",dim(genes)[2]))

print(paste0("# microbes = ",dim(microbes)[2]))

## Ensure same sampleIDs in both genes and microbes matrices
stopifnot(all(rownames(genes) == rownames(microbes)))

## convert both dataframes to matrix. cv.glmnet expects a matrix of predictors, not a data frame
y <- as.matrix(genes) #response
x <- as.matrix(microbes) #predictors


#colSums(x!=0)%>%sort()
#colSums(x!=0)>0.05*nrow(x)

x=x[,colSums(x!=0)>0.05*nrow(x)]
dim(x)

##################### Stability selection ####################

stabs_stability <- function( x, y, gene_name){
  
  ## Import all the libraries for the current node/core on MSI 
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("stabs","glmnet","methods","doParallel") ## package methods is not loaded by default by RScript. 
  check.packages(packages)
  
  ## Extract the expression for this gene (response variable)
  y_i <- y[,grep(paste0("^",gene_name,"$"),colnames(y))]
  
  ## Make sure y_i is numeric before model fitting 
  stopifnot(class(y_i) == "numeric")
  
  ## perform stability selection using glmnet lasso
  stab.glmnet <- stabsel(x = x, y = y_i,
                         fitfun = glmnet.lasso, cutoff = 0.6,
                         PFER = 1)
  
  taxa.selected <- names(stab.glmnet$selected)
  if(length(taxa.selected) == 0) taxa.selected <-"None"
  
  taxa.selected.df <- as.data.frame(cbind(gene_name,taxa.selected))
  
  return(taxa.selected.df)
  
}

## Parallel run on MSI
gene_list <- scan(input.genes.list, what="", sep="\n")

## invoke parallel computation for fitting model for each gene
parallel_time <- system.time({
  # sink(paste0(output.dir,"/log.txt"), append=TRUE)
  parallel_res <- foreach(i=1:length(gene_list), .combine = rbind) %dopar% {
    # cat(paste("Processing gene:",gene_list[i],"\n"))
    # sink() #end diversion of output
    stabs_stability(x,y,gene_list[i])
  }
})
stopCluster(cl)
registerDoSEQ()

## time taken for parallel computation.
print(parallel_time)

filename <- strsplit(input.genes.list,"/")
filename <- strsplit(filename[[1]][length(filename[[1]])],"\\.")[[1]][1]
write.table(parallel_res, file=paste0(output.dir,"/",filename,"_stabs_stabsel_output.txt"), quote=F, sep="\t",col.names = NA, row.names = TRUE)




#-----------------------stabs_stability_selection-----------------------------

#--------------------lncRNA----------CD --------------------


######################### Initial setup ######################
## before initiating parallel 
rm(list=ls()) ## clear workspace

## Setup for parallel processing on MSI
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("doParallel", "data.table") ## package methods is not loaded by default by RScript on MSI
check.packages(packages)

## now initiate parallel processors.
cl <- makeCluster(60) # 8 workers per core -- itasca
registerDoParallel(cl)
print(paste0("Number of workers: ", getDoParWorkers()))



################# Input genes and taxa matrix ###########
input.genes.list="/mnt/data/IBD/RNA_seq/data/Output/Others/14/gene_list2.txt"
output.dir="/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/stabsel_path"

load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/LncRNA_CD.Rdata")
load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/microbes_CD.Rdata")


genes <-LncRNA_CD
microbes <- microbes_CD

print(paste0("# genes = ",dim(genes)[2]))

print(paste0("# microbes = ",dim(microbes)[2]))

## Ensure same sampleIDs in both genes and microbes matrices
stopifnot(all(rownames(genes) == rownames(microbes)))

## convert both dataframes to matrix. cv.glmnet expects a matrix of predictors, not a data frame
y <- as.matrix(genes) #response
x <- as.matrix(microbes) #predictors


##################### Stability selection ####################

stabs_stability <- function( x, y, gene_name){
  
  ## Import all the libraries for the current node/core on MSI 
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("stabs","glmnet","methods","doParallel") ## package methods is not loaded by default by RScript. 
  check.packages(packages)
  
  ## Extract the expression for this gene (response variable)
  y_i <- y[,grep(paste0("^",gene_name,"$"),colnames(y))]
  
  ## Make sure y_i is numeric before model fitting 
  stopifnot(class(y_i) == "numeric")
  
  ## perform stability selection using glmnet lasso
  stab.glmnet <- stabsel(x = x, y = y_i,
                         fitfun = glmnet.lasso, cutoff = 0.6,
                         PFER = 1)
  
  taxa.selected <- names(stab.glmnet$selected)
  if(length(taxa.selected) == 0) taxa.selected <-"None"
  
  taxa.selected.df <- as.data.frame(cbind(gene_name,taxa.selected))
  
  return(taxa.selected.df)
  
}

## Parallel run on MSI
gene_list <- scan(input.genes.list, what="", sep="\n")

## invoke parallel computation for fitting model for each gene
parallel_time <- system.time({
  # sink(paste0(output.dir,"/log.txt"), append=TRUE)
  parallel_res <- foreach(i=1:length(gene_list), .combine = rbind) %dopar% {
    # cat(paste("Processing gene:",gene_list[i],"\n"))
    # sink() #end diversion of output
    stabs_stability(x,y,gene_list[i])
  }
})
stopCluster(cl)
registerDoSEQ()

## time taken for parallel computation.
print(parallel_time)

filename <- strsplit(input.genes.list,"/")
filename <- strsplit(filename[[1]][length(filename[[1]])],"\\.")[[1]][1]
write.table(parallel_res, file=paste0(output.dir,"/",filename,"_stabs_stabsel_output.txt"), quote=F, sep="\t",col.names = NA, row.names = TRUE)



#--------------------lncRNA----------Control --------------------


######################### Initial setup ######################
## before initiating parallel 
rm(list=ls()) ## clear workspace

## Setup for parallel processing on MSI
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("doParallel", "data.table") ## package methods is not loaded by default by RScript on MSI
check.packages(packages)

## now initiate parallel processors.
cl <- makeCluster(60) # 8 workers per core -- itasca
registerDoParallel(cl)
print(paste0("Number of workers: ", getDoParWorkers()))



################# Input genes and taxa matrix ###########
input.genes.list="/mnt/data/IBD/RNA_seq/data/Output/Others/14/gene_list2.txt"
output.dir="/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/stabsel_path"

load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/LncRNA_Contr.Rdata")
load("/mnt/data/IBD/RNA_seq/data/Output/Others/14/microbes_Contr.Rdata")


genes <-LncRNA_Contr
microbes <- microbes_Contr

print(paste0("# genes = ",dim(genes)[2]))

print(paste0("# microbes = ",dim(microbes)[2]))

## Ensure same sampleIDs in both genes and microbes matrices
stopifnot(all(rownames(genes) == rownames(microbes)))

## convert both dataframes to matrix. cv.glmnet expects a matrix of predictors, not a data frame
y <- as.matrix(genes) #response
x <- as.matrix(microbes) #predictors


#colSums(x!=0)%>%sort()
#colSums(x!=0)>0.05*nrow(x)

x=x[,colSums(x!=0)>0.05*nrow(x)]
dim(x)

##################### Stability selection ####################

stabs_stability <- function( x, y, gene_name){
  
  ## Import all the libraries for the current node/core on MSI 
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repos = "http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  
  packages <- c("stabs","glmnet","methods","doParallel") ## package methods is not loaded by default by RScript. 
  check.packages(packages)
  
  ## Extract the expression for this gene (response variable)
  y_i <- y[,grep(paste0("^",gene_name,"$"),colnames(y))]
  
  ## Make sure y_i is numeric before model fitting 
  stopifnot(class(y_i) == "numeric")
  
  ## perform stability selection using glmnet lasso
  stab.glmnet <- stabsel(x = x, y = y_i,
                         fitfun = glmnet.lasso, cutoff = 0.6,
                         PFER = 1)
  
  taxa.selected <- names(stab.glmnet$selected)
  if(length(taxa.selected) == 0) taxa.selected <-"None"
  
  taxa.selected.df <- as.data.frame(cbind(gene_name,taxa.selected))
  
  return(taxa.selected.df)
  
}

## Parallel run on MSI
gene_list <- scan(input.genes.list, what="", sep="\n")

## invoke parallel computation for fitting model for each gene
parallel_time <- system.time({
  # sink(paste0(output.dir,"/log.txt"), append=TRUE)
  parallel_res <- foreach(i=1:length(gene_list), .combine = rbind) %dopar% {
    # cat(paste("Processing gene:",gene_list[i],"\n"))
    # sink() #end diversion of output
    stabs_stability(x,y,gene_list[i])
  }
})
stopCluster(cl)
registerDoSEQ()

## time taken for parallel computation.
print(parallel_time)

filename <- strsplit(input.genes.list,"/")
filename <- strsplit(filename[[1]][length(filename[[1]])],"\\.")[[1]][1]
write.table(parallel_res, file=paste0(output.dir,"/",filename,"_stabs_stabsel_output.txt"), quote=F, sep="\t",col.names = NA, row.names = TRUE)
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------mRNA------CD--------------

## Initial setup
rm(list=ls()) ##check workspace to make sure no precomputed data before clearing
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(RColorBrewer)
library(gridExtra)
library(corrplot)
library(data.table)


################# Functions ###################
## massage taxa names
massage_taxaname <- function(taxa){
  # microbe_name <- colnames(gene_x_microbe)[4] ## for debugging
  microbe_name <- taxa
  ## remove k__ or p__ etc. from the string
  microbe_name <- gsub("[A-z]__","", microbe_name)
  ## remove trailing ;
  microbe_name <- gsub("\\;+$","", microbe_name)
  return(microbe_name)
}


## Function to combine all the gene outputs together
get.combined.output <- function(filepath){
  filenames <- list.files(filepath, pattern="*.txt")
  summary_list <- list()
  count <- 1 
  for( file in filenames){
    # df <- read.table(file = paste0(filepath,"/",file), sep="\t",head=T, row.names = 1, check.names = F)
    df <- data.frame(fread(paste0(filepath,"/",file), sep="\t", head=T), row.names = 1, check.names = F)
    summary_list[[count]] <- df
    count <- count + 1
  }
  all_genes_result_df <- do.call(rbind,summary_list)
}

## correct for multiple hypothesis testing
MHT.correction <- function(df){
  ## Bonferonni
  assoc_Bonferonni <- p.adjust(df$pval, method = "bonferroni")
  df$Bonferonni_global <- assoc_Bonferonni
  ## BY
  assoc_BY <- p.adjust(df$pval, method="BY")
  df$BY_global <- assoc_BY
  ## BH 
  assoc_BH <- p.adjust(df$pval, method="BH")
  df$BH_global <- assoc_BH
  
  return(df)
}

create_scatterplot_IBD <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  ## add IBS_subytpe info before plotting
  IBS_subtype <- taxa.table$IBS_subtype + 1 ## since 0 is not accessible in R
  df <- data.frame(microbe,gene,IBS_subtype,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    ## Note, colors are set by IBS_subtype
    geom_point(fill= col.fill[IBS_subtype], colour=border.col[IBS_subtype], shape=21, size=3, stroke = 1)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  
  return(scatterplot)
}

create_scatterplot_Healthy <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  df <- data.frame(microbe,gene,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill, colour=border.col, shape=21, size=3, stroke = 1)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}


## scatterplot function
create.scatterplot <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  df <- data.frame(microbe,gene,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill, colour=border.col, shape=21, size=3, stroke = 2)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}

create.scatterplot.2 <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex)
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  condition <- taxa.table$condition + 1 ## since 0 is not accessible in R
  df <- data.frame(microbe,gene,condition,pval,padj)
  
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill[condition], colour=border.col, shape=21, size=3, stroke = 2)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}


############### Combine MSI output files (SKIP IF ALREADY DONE) ################

## Input dataset
dataset <- "IBD/case/"

combined.output <- get.combined.output("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/lasso_hdi_loocv/")
dim(combined.output)
write.table(combined.output, file ="/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/combined_output/combined_lasso_output.txt", sep="\t", row.names = T, col.names = NA )


combined.output.filepath <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/combined_output/combined_lasso_output.txt"
## read.table() takes too long! Replace with fread(), see below.
## fread() to read large file. 
combined.output <- data.frame(fread(combined.output.filepath, sep="\t", head=T), row.names = 1)
dim(combined.output)


########### Perform MHT correction on p-values for all assoc. ############

combined.output <- MHT.correction(combined.output)
## Order output by BH_global
combined.output <- combined.output[order(combined.output$BH_global),]


## write to file
filename <- "combined_lasso_output_MHT_corrected.txt"
write.table(combined.output, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/combined_output/",filename), sep="\t", row.names = T, col.names = NA )

################# SKIP IF DONE -- Combine output from MSI stability selection run ##############
## Combine all the gene outputs together
stabsel_path <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/stabsel_path/"
stability.combined.output <- get.combined.output(stabsel_path)
dim(stability.combined.output)
length(unique(stability.combined.output$gene_name)) 



## Number of associations with a stability selected covariate (taxa or gender or condition)
dim(stability.combined.output[stability.combined.output$taxa.selected != "None",])

## Number of genes with stability selected covariate (taxa or gender or condition)
length(unique(stability.combined.output[stability.combined.output$taxa.selected != "None",]$gene_name))

## Number of genes with stability selected taxa
length(unique(stability.combined.output[stability.combined.output$taxa.selected != "None" & stability.combined.output$taxa.selected != "gender",]$gene_name)) 


## write to file
stabsel.combined.filename <- "combined_stabsel.txt"
write.table(stability.combined.output, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/combined_output/", stabsel.combined.filename), sep="\t", row.names = T, col.names = NA )



################# Test overlap with stability selected output ###############

## Input dataset and stabsel output filename
stabsel_path <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/combined_output/"
stabsel_file <- "combined_stabsel.txt"

stabsel_output_filepath <- paste0(stabsel_path,stabsel_file)
stabsel.output <- read.table(stabsel_output_filepath,sep="\t",head=T, row.names = 1, check.names = F)
dim(stabsel.output)




## massage column names of stability selected output for downstream processing
colnames(stabsel.output) <- c("gene","taxa")

## Identify stability selected associations
stabsel.output.filt <- stabsel.output[stabsel.output$taxa != "None",] 
dim(stabsel.output.filt)

## Identify unique stability selected genes
genes.stabsel <- unique(stabsel.output.filt$gene)
length(genes.stabsel)

## Identify unique stability selected predictors
taxa.stabsel <- unique(stabsel.output.filt$taxa)
length(taxa.stabsel)


## Overlap between combined.output.filt and stability selected assoc
overlap.assoc.stabsel <- merge(combined.output,stabsel.output.filt, by = c("gene","taxa"))
overlap.assoc.stabsel <- overlap.assoc.stabsel[order(overlap.assoc.stabsel$BH_global),]
dim(overlap.assoc.stabsel)

length(unique(overlap.assoc.stabsel$gene))

length(unique(overlap.assoc.stabsel$taxa))


## write to file before filtering by FDR cutoffs
filename <- "combined_output_MHT_stabsel.txt"
write.table(overlap.assoc.stabsel, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/combined_output/",filename), sep="\t", row.names = T, col.names = NA )

## How many stability selected assoc. also significant at FDR < 0.1
overlap.assoc.stabsel.BH.0.1 <- overlap.assoc.stabsel[overlap.assoc.stabsel$BH_global < 0.1,]
dim(overlap.assoc.stabsel.BH.0.1)

length(unique(overlap.assoc.stabsel.BH.0.1$gene))

length(unique(overlap.assoc.stabsel.BH.0.1$taxa)) 


############# Filter out non-taxa associations (e.g. gene-gender assoc.)######################

combined.output <- overlap.assoc.stabsel.BH.0.1

## Identify gender associations
select <- which(combined.output$taxa =="gender")
## Spit out gender associated output to a file for later inspection
gender.assoc <- combined.output[select,]; dim(gender.assoc)

length(unique(gender.assoc$gene))

## write to file
write.table(gender.assoc, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/combined_output/gene_gender_FDR_0.1.txt"), sep="\t", row.names = T, col.names = NA )

## filter out gender associations.
combined.output.filt <- combined.output[-select,]
dim(combined.output.filt)

## Additional gene-covariate filtering for IBD or IBS
if(dataset == "IBD/case/" || dataset == "IBS/case/"){
  ## set subtype label
  subtype_covariate <- "IBD_subtype" #IBD_subtype
  select <- which(combined.output.filt$taxa == subtype_covariate)
  condition.assoc <- combined.output.filt[select,]; dim(condition.assoc)
  
  length(unique(condition.assoc$gene)) 
  
  # write gene-condition associations to file
  # write.table(condition.assoc, file = paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V2/combined_output/gene_IBS_subtype_FDR_0.1.txt"), sep="\t", row.names = T, col.names = NA )
  
  # #keep non-condition associations.
  if(length(select) != 0)
    combined.output.filt <- combined.output.filt[-select,]; dim(combined.output.filt)
}




## no. of unique genes
length(unique(combined.output.filt$gene))

## no. of unique taxa
length(unique(combined.output.filt$taxa))

## set datatset and filename and write to file
dataset <- "IBD/control/"
filename <- "IBD_case_gene_taxa_FDR_0.1.txt"
write.table(combined.output.filt, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/combined_output/",filename), sep="\t", row.names = T, col.names = NA )



#---------------------------------mRNA------Control--------------

## Initial setup
rm(list=ls()) ##check workspace to make sure no precomputed data before clearing
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(RColorBrewer)
library(gridExtra)
library(corrplot)
library(data.table)


################# Functions ###################
## massage taxa names
massage_taxaname <- function(taxa){
  # microbe_name <- colnames(gene_x_microbe)[4] ## for debugging
  microbe_name <- taxa
  ## remove k__ or p__ etc. from the string
  microbe_name <- gsub("[A-z]__","", microbe_name)
  ## remove trailing ;
  microbe_name <- gsub("\\;+$","", microbe_name)
  return(microbe_name)
}


## Function to combine all the gene outputs together
get.combined.output <- function(filepath){
  filenames <- list.files(filepath, pattern="*.txt")
  summary_list <- list()
  count <- 1 
  for( file in filenames){
    # df <- read.table(file = paste0(filepath,"/",file), sep="\t",head=T, row.names = 1, check.names = F)
    df <- data.frame(fread(paste0(filepath,"/",file), sep="\t", head=T), row.names = 1, check.names = F)
    summary_list[[count]] <- df
    count <- count + 1
  }
  all_genes_result_df <- do.call(rbind,summary_list)
}

## correct for multiple hypothesis testing
MHT.correction <- function(df){
  ## Bonferonni
  assoc_Bonferonni <- p.adjust(df$pval, method = "bonferroni")
  df$Bonferonni_global <- assoc_Bonferonni
  ## BY
  assoc_BY <- p.adjust(df$pval, method="BY")
  df$BY_global <- assoc_BY
  ## BH 
  assoc_BH <- p.adjust(df$pval, method="BH")
  df$BH_global <- assoc_BH
  
  return(df)
}

create_scatterplot_IBD <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  ## add IBS_subytpe info before plotting
  IBS_subtype <- taxa.table$IBS_subtype + 1 ## since 0 is not accessible in R
  df <- data.frame(microbe,gene,IBS_subtype,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    ## Note, colors are set by IBS_subtype
    geom_point(fill= col.fill[IBS_subtype], colour=border.col[IBS_subtype], shape=21, size=3, stroke = 1)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  
  return(scatterplot)
}

create_scatterplot_Healthy <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  df <- data.frame(microbe,gene,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill, colour=border.col, shape=21, size=3, stroke = 1)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}


## scatterplot function
create.scatterplot <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  df <- data.frame(microbe,gene,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill, colour=border.col, shape=21, size=3, stroke = 2)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}

create.scatterplot.2 <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex)
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  condition <- taxa.table$condition + 1 ## since 0 is not accessible in R
  df <- data.frame(microbe,gene,condition,pval,padj)
  
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill[condition], colour=border.col, shape=21, size=3, stroke = 2)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}


############### Combine MSI output files (SKIP IF ALREADY DONE) ################

## Input dataset
dataset <- "IBD/control/"

combined.output <- get.combined.output("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/lasso_hdi_loocv/")
dim(combined.output)
write.table(combined.output, file ="/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/combined_output/combined_lasso_output.txt", sep="\t", row.names = T, col.names = NA )


combined.output.filepath <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/combined_output/combined_lasso_output.txt"
## read.table() takes too long! Replace with fread(), see below.
## fread() to read large file. 
combined.output <- data.frame(fread(combined.output.filepath, sep="\t", head=T), row.names = 1)
dim(combined.output)


########### Perform MHT correction on p-values for all assoc. ############

combined.output <- MHT.correction(combined.output)
## Order output by BH_global
combined.output <- combined.output[order(combined.output$BH_global),]


## write to file
filename <- "combined_lasso_output_MHT_corrected.txt"
write.table(combined.output, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/combined_output/",filename), sep="\t", row.names = T, col.names = NA )

################# SKIP IF DONE -- Combine output from MSI stability selection run ##############
## Combine all the gene outputs together
stabsel_path <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/stabsel_path/"
stability.combined.output <- get.combined.output(stabsel_path)
dim(stability.combined.output)
length(unique(stability.combined.output$gene_name)) 



## Number of associations with a stability selected covariate (taxa or gender or condition)
dim(stability.combined.output[stability.combined.output$taxa.selected != "None",])

## Number of genes with stability selected covariate (taxa or gender or condition)
length(unique(stability.combined.output[stability.combined.output$taxa.selected != "None",]$gene_name))

## Number of genes with stability selected taxa
length(unique(stability.combined.output[stability.combined.output$taxa.selected != "None" & stability.combined.output$taxa.selected != "gender",]$gene_name)) 


## write to file
stabsel.combined.filename <- "combined_stabsel.txt"
write.table(stability.combined.output, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/combined_output/", stabsel.combined.filename), sep="\t", row.names = T, col.names = NA )



################# Test overlap with stability selected output ###############

## Input dataset and stabsel output filename
stabsel_path <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/combined_output/"
stabsel_file <- "combined_stabsel.txt"

stabsel_output_filepath <- paste0(stabsel_path,stabsel_file)
stabsel.output <- read.table(stabsel_output_filepath,sep="\t",head=T, row.names = 1, check.names = F)
dim(stabsel.output)




## massage column names of stability selected output for downstream processing
colnames(stabsel.output) <- c("gene","taxa")

## Identify stability selected associations
stabsel.output.filt <- stabsel.output[stabsel.output$taxa != "None",] 
dim(stabsel.output.filt)

## Identify unique stability selected genes
genes.stabsel <- unique(stabsel.output.filt$gene)
length(genes.stabsel)

## Identify unique stability selected predictors
taxa.stabsel <- unique(stabsel.output.filt$taxa)
length(taxa.stabsel)


## Overlap between combined.output.filt and stability selected assoc
overlap.assoc.stabsel <- merge(combined.output,stabsel.output.filt, by = c("gene","taxa"))
overlap.assoc.stabsel <- overlap.assoc.stabsel[order(overlap.assoc.stabsel$BH_global),]
dim(overlap.assoc.stabsel)

length(unique(overlap.assoc.stabsel$gene))

length(unique(overlap.assoc.stabsel$taxa))


## write to file before filtering by FDR cutoffs
filename <- "combined_output_MHT_stabsel.txt"
write.table(overlap.assoc.stabsel, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/combined_output/",filename), sep="\t", row.names = T, col.names = NA )

## How many stability selected assoc. also significant at FDR < 0.1
overlap.assoc.stabsel.BH.0.1 <- overlap.assoc.stabsel[overlap.assoc.stabsel$BH_global < 0.1,]
dim(overlap.assoc.stabsel.BH.0.1)

length(unique(overlap.assoc.stabsel.BH.0.1$gene))

length(unique(overlap.assoc.stabsel.BH.0.1$taxa)) 


############# Filter out non-taxa associations (e.g. gene-gender assoc.)######################

combined.output <- overlap.assoc.stabsel.BH.0.1

## Identify gender associations
select <- which(combined.output$taxa =="gender")
## Spit out gender associated output to a file for later inspection
gender.assoc <- combined.output[select,]; dim(gender.assoc)

length(unique(gender.assoc$gene))

## write to file
write.table(gender.assoc, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/combined_output/gene_gender_FDR_0.1.txt"), sep="\t", row.names = T, col.names = NA )

## filter out gender associations.
combined.output.filt <- combined.output[-select,]
dim(combined.output.filt)

## Additional gene-covariate filtering for IBD or IBS
if(dataset == "IBD/control/" || dataset == "IBS/control/"){
  ## set subtype label
  subtype_covariate <- "IBD_subtype" #IBD_subtype
  select <- which(combined.output.filt$taxa == subtype_covariate)
  condition.assoc <- combined.output.filt[select,]; dim(condition.assoc)
  
  length(unique(condition.assoc$gene)) 
  
  # write gene-condition associations to file
  # write.table(condition.assoc, file = paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V2/combined_output/gene_IBS_subtype_FDR_0.1.txt"), sep="\t", row.names = T, col.names = NA )
  
  # #keep non-condition associations.
  if(length(select) != 0)
    combined.output.filt <- combined.output.filt[-select,]; dim(combined.output.filt)
}




## no. of unique genes
length(unique(combined.output.filt$gene))

## no. of unique taxa
length(unique(combined.output.filt$taxa))

## set datatset and filename and write to file
dataset <- "IBD/control/"
filename <- "IBD_control_gene_taxa_FDR_0.1.txt"
write.table(combined.output.filt, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/combined_output/",filename), sep="\t", row.names = T, col.names = NA )





#---------------------------------lncRNA------CD--------------

## Initial setup
rm(list=ls()) ##check workspace to make sure no precomputed data before clearing
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(RColorBrewer)
library(gridExtra)
library(corrplot)
library(data.table)


################# Functions ###################
## massage taxa names
massage_taxaname <- function(taxa){
  # microbe_name <- colnames(gene_x_microbe)[4] ## for debugging
  microbe_name <- taxa
  ## remove k__ or p__ etc. from the string
  microbe_name <- gsub("[A-z]__","", microbe_name)
  ## remove trailing ;
  microbe_name <- gsub("\\;+$","", microbe_name)
  return(microbe_name)
}


## Function to combine all the gene outputs together
get.combined.output <- function(filepath){
  filenames <- list.files(filepath, pattern="*.txt")
  summary_list <- list()
  count <- 1 
  for( file in filenames){
    # df <- read.table(file = paste0(filepath,"/",file), sep="\t",head=T, row.names = 1, check.names = F)
    df <- data.frame(fread(paste0(filepath,"/",file), sep="\t", head=T), row.names = 1, check.names = F)
    summary_list[[count]] <- df
    count <- count + 1
  }
  all_genes_result_df <- do.call(rbind,summary_list)
}

## correct for multiple hypothesis testing
MHT.correction <- function(df){
  ## Bonferonni
  assoc_Bonferonni <- p.adjust(df$pval, method = "bonferroni")
  df$Bonferonni_global <- assoc_Bonferonni
  ## BY
  assoc_BY <- p.adjust(df$pval, method="BY")
  df$BY_global <- assoc_BY
  ## BH 
  assoc_BH <- p.adjust(df$pval, method="BH")
  df$BH_global <- assoc_BH
  
  return(df)
}

create_scatterplot_IBD <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  ## add IBS_subytpe info before plotting
  IBS_subtype <- taxa.table$IBS_subtype + 1 ## since 0 is not accessible in R
  df <- data.frame(microbe,gene,IBS_subtype,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    ## Note, colors are set by IBS_subtype
    geom_point(fill= col.fill[IBS_subtype], colour=border.col[IBS_subtype], shape=21, size=3, stroke = 1)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  
  return(scatterplot)
}

create_scatterplot_Healthy <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  df <- data.frame(microbe,gene,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill, colour=border.col, shape=21, size=3, stroke = 1)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}


## scatterplot function
create.scatterplot <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  df <- data.frame(microbe,gene,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill, colour=border.col, shape=21, size=3, stroke = 2)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}

create.scatterplot.2 <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex)
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  condition <- taxa.table$condition + 1 ## since 0 is not accessible in R
  df <- data.frame(microbe,gene,condition,pval,padj)
  
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill[condition], colour=border.col, shape=21, size=3, stroke = 2)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}


############### Combine MSI output files (SKIP IF ALREADY DONE) ################

## Input dataset
dataset <- "IBD/case/"

combined.output <- get.combined.output("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/lasso_hdi_loocv/")
dim(combined.output)
write.table(combined.output, file ="/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/combined_output/combined_lasso_output.txt", sep="\t", row.names = T, col.names = NA )


combined.output.filepath <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/combined_output/combined_lasso_output.txt"
## read.table() takes too long! Replace with fread(), see below.
## fread() to read large file. 
combined.output <- data.frame(fread(combined.output.filepath, sep="\t", head=T), row.names = 1)
dim(combined.output)


########### Perform MHT correction on p-values for all assoc. ############

combined.output <- MHT.correction(combined.output)
## Order output by BH_global
combined.output <- combined.output[order(combined.output$BH_global),]


## write to file
filename <- "combined_lasso_output_MHT_corrected.txt"
write.table(combined.output, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/combined_output/",filename), sep="\t", row.names = T, col.names = NA )

################# SKIP IF DONE -- Combine output from MSI stability selection run ##############
## Combine all the gene outputs together
stabsel_path <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/stabsel_path/"
stability.combined.output <- get.combined.output(stabsel_path)
dim(stability.combined.output)
length(unique(stability.combined.output$gene_name)) 



## Number of associations with a stability selected covariate (taxa or gender or condition)
dim(stability.combined.output[stability.combined.output$taxa.selected != "None",])

## Number of genes with stability selected covariate (taxa or gender or condition)
length(unique(stability.combined.output[stability.combined.output$taxa.selected != "None",]$gene_name))

## Number of genes with stability selected taxa
length(unique(stability.combined.output[stability.combined.output$taxa.selected != "None" & stability.combined.output$taxa.selected != "gender",]$gene_name)) 


## write to file
stabsel.combined.filename <- "combined_stabsel.txt"
write.table(stability.combined.output, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/combined_output/", stabsel.combined.filename), sep="\t", row.names = T, col.names = NA )



################# Test overlap with stability selected output ###############

## Input dataset and stabsel output filename
stabsel_path <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/combined_output/"
stabsel_file <- "combined_stabsel.txt"

stabsel_output_filepath <- paste0(stabsel_path,stabsel_file)
stabsel.output <- read.table(stabsel_output_filepath,sep="\t",head=T, row.names = 1, check.names = F)
dim(stabsel.output)




## massage column names of stability selected output for downstream processing
colnames(stabsel.output) <- c("gene","taxa")

## Identify stability selected associations
stabsel.output.filt <- stabsel.output[stabsel.output$taxa != "None",] 
dim(stabsel.output.filt)

## Identify unique stability selected genes
genes.stabsel <- unique(stabsel.output.filt$gene)
length(genes.stabsel)

## Identify unique stability selected predictors
taxa.stabsel <- unique(stabsel.output.filt$taxa)
length(taxa.stabsel)


## Overlap between combined.output.filt and stability selected assoc
overlap.assoc.stabsel <- merge(combined.output,stabsel.output.filt, by = c("gene","taxa"))
overlap.assoc.stabsel <- overlap.assoc.stabsel[order(overlap.assoc.stabsel$BH_global),]
dim(overlap.assoc.stabsel)

length(unique(overlap.assoc.stabsel$gene))

length(unique(overlap.assoc.stabsel$taxa))


## write to file before filtering by FDR cutoffs
filename <- "combined_output_MHT_stabsel.txt"
write.table(overlap.assoc.stabsel, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/combined_output/",filename), sep="\t", row.names = T, col.names = NA )

## How many stability selected assoc. also significant at FDR < 0.1
overlap.assoc.stabsel.BH.0.1 <- overlap.assoc.stabsel[overlap.assoc.stabsel$BH_global < 0.1,]
dim(overlap.assoc.stabsel.BH.0.1)

length(unique(overlap.assoc.stabsel.BH.0.1$gene))

length(unique(overlap.assoc.stabsel.BH.0.1$taxa)) 


############# Filter out non-taxa associations (e.g. gene-gender assoc.)######################

combined.output <- overlap.assoc.stabsel.BH.0.1

## Identify gender associations
select <- which(combined.output$taxa =="gender")
## Spit out gender associated output to a file for later inspection
gender.assoc <- combined.output[select,]; dim(gender.assoc)

length(unique(gender.assoc$gene))

## write to file
write.table(gender.assoc, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/combined_output/gene_gender_FDR_0.1.txt"), sep="\t", row.names = T, col.names = NA )

## filter out gender associations.
combined.output.filt <- combined.output[-select,]
dim(combined.output.filt)

## Additional gene-covariate filtering for IBD or IBS
if(dataset == "IBD/case/" || dataset == "IBS/case/"){
  ## set subtype label
  subtype_covariate <- "IBD_subtype" #IBD_subtype
  select <- which(combined.output.filt$taxa == subtype_covariate)
  condition.assoc <- combined.output.filt[select,]; dim(condition.assoc)
  
  length(unique(condition.assoc$gene)) 
  
  # write gene-condition associations to file
  # write.table(condition.assoc, file = paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V2/combined_output/gene_IBS_subtype_FDR_0.1.txt"), sep="\t", row.names = T, col.names = NA )
  
  # #keep non-condition associations.
  if(length(select) != 0)
    combined.output.filt <- combined.output.filt[-select,]; dim(combined.output.filt)
}




## no. of unique genes
length(unique(combined.output.filt$gene))

## no. of unique taxa
length(unique(combined.output.filt$taxa))

## set datatset and filename and write to file
dataset <- "IBD/control/"
filename <- "IBD_case_gene_taxa_FDR_0.1.txt"
write.table(combined.output.filt, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/CD/combined_output/",filename), sep="\t", row.names = T, col.names = NA )


#---------------------------------lncRNA------Control--------------

## Initial setup
rm(list=ls()) ##check workspace to make sure no precomputed data before clearing
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(RColorBrewer)
library(gridExtra)
library(corrplot)
library(data.table)


################# Functions ###################
## massage taxa names
massage_taxaname <- function(taxa){
  # microbe_name <- colnames(gene_x_microbe)[4] ## for debugging
  microbe_name <- taxa
  ## remove k__ or p__ etc. from the string
  microbe_name <- gsub("[A-z]__","", microbe_name)
  ## remove trailing ;
  microbe_name <- gsub("\\;+$","", microbe_name)
  return(microbe_name)
}


## Function to combine all the gene outputs together
get.combined.output <- function(filepath){
  filenames <- list.files(filepath, pattern="*.txt")
  summary_list <- list()
  count <- 1 
  for( file in filenames){
    # df <- read.table(file = paste0(filepath,"/",file), sep="\t",head=T, row.names = 1, check.names = F)
    df <- data.frame(fread(paste0(filepath,"/",file), sep="\t", head=T), row.names = 1, check.names = F)
    summary_list[[count]] <- df
    count <- count + 1
  }
  all_genes_result_df <- do.call(rbind,summary_list)
}

## correct for multiple hypothesis testing
MHT.correction <- function(df){
  ## Bonferonni
  assoc_Bonferonni <- p.adjust(df$pval, method = "bonferroni")
  df$Bonferonni_global <- assoc_Bonferonni
  ## BY
  assoc_BY <- p.adjust(df$pval, method="BY")
  df$BY_global <- assoc_BY
  ## BH 
  assoc_BH <- p.adjust(df$pval, method="BH")
  df$BH_global <- assoc_BH
  
  return(df)
}

create_scatterplot_IBD <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  ## add IBS_subytpe info before plotting
  IBS_subtype <- taxa.table$IBS_subtype + 1 ## since 0 is not accessible in R
  df <- data.frame(microbe,gene,IBS_subtype,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    ## Note, colors are set by IBS_subtype
    geom_point(fill= col.fill[IBS_subtype], colour=border.col[IBS_subtype], shape=21, size=3, stroke = 1)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  
  return(scatterplot)
}

create_scatterplot_Healthy <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  df <- data.frame(microbe,gene,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill, colour=border.col, shape=21, size=3, stroke = 1)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}


## scatterplot function
create.scatterplot <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe/metabolite names e.g s__Anaerotruncus_sp._G3(2012), melatonin (pg/mg)
  microbe_name_regex <- gsub("\\(","\\\\(",microbe_name, perl = T)
  microbe_name_regex <- gsub("\\)","\\\\)",microbe_name_regex, perl = T)
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name_regex, perl = T)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex, perl = T)
  
  
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  df <- data.frame(microbe,gene,pval,padj)
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill, colour=border.col, shape=21, size=3, stroke = 2)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}

create.scatterplot.2 <- function(genes.table,taxa.table,microbe_name,gene_name,pval,padj,col.fill, border.col){
  
  
  gene <- as.numeric(genes.table[,grep(paste0("^",gene_name,"$"),colnames(genes.table))]) ## sample x genes
  ## deal with special characters in microbe names e.g [ ]
  microbe_name_regex <- gsub("\\[","\\\\[",microbe_name)
  microbe_name_regex <- gsub("\\]","\\\\]",microbe_name_regex)
  microbe <- as.numeric(taxa.table[,grep(paste0("^",microbe_name_regex,"$"),colnames(taxa.table))]) ## samples x taxa
  condition <- taxa.table$condition + 1 ## since 0 is not accessible in R
  df <- data.frame(microbe,gene,condition,pval,padj)
  
  
  ## Massage taxaname before plotting
  microbe_name <- massage_taxaname(as.character(microbe_name))
  microbe_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
  
  scatterplot <- ggplot(df,aes(x=microbe,y=gene)) +
    ## For 3x3 output
    geom_point(fill= col.fill[condition], colour=border.col, shape=21, size=3, stroke = 2)+
    geom_smooth(method="lm",se=TRUE, color ="black")+
    labs(x=microbe_name,y=gene_name)+ #x-axis and y-axis labels
    theme_bw()+
    ## For 3x3 output
    ggtitle(paste0("pval=",pval,", padj=",padj))+
    theme(plot.title = element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10, face = "italic"))
  
  return(scatterplot)
}


############### Combine MSI output files (SKIP IF ALREADY DONE) ################

## Input dataset
dataset <- "IBD/control/"

combined.output <- get.combined.output("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/lasso_hdi_loocv/")
dim(combined.output)
write.table(combined.output, file ="/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/combined_output/combined_lasso_output.txt", sep="\t", row.names = T, col.names = NA )


combined.output.filepath <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/combined_output/combined_lasso_output.txt"
## read.table() takes too long! Replace with fread(), see below.
## fread() to read large file. 
combined.output <- data.frame(fread(combined.output.filepath, sep="\t", head=T), row.names = 1)
dim(combined.output)


########### Perform MHT correction on p-values for all assoc. ############

combined.output <- MHT.correction(combined.output)
## Order output by BH_global
combined.output <- combined.output[order(combined.output$BH_global),]


## write to file
filename <- "combined_lasso_output_MHT_corrected.txt"
write.table(combined.output, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/combined_output/",filename), sep="\t", row.names = T, col.names = NA )

################# SKIP IF DONE -- Combine output from MSI stability selection run ##############
## Combine all the gene outputs together
stabsel_path <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/stabsel_path/"
stability.combined.output <- get.combined.output(stabsel_path)
dim(stability.combined.output)
length(unique(stability.combined.output$gene_name)) 



## Number of associations with a stability selected covariate (taxa or gender or condition)
dim(stability.combined.output[stability.combined.output$taxa.selected != "None",])

## Number of genes with stability selected covariate (taxa or gender or condition)
length(unique(stability.combined.output[stability.combined.output$taxa.selected != "None",]$gene_name))

## Number of genes with stability selected taxa
length(unique(stability.combined.output[stability.combined.output$taxa.selected != "None" & stability.combined.output$taxa.selected != "gender",]$gene_name)) 


## write to file
stabsel.combined.filename <- "combined_stabsel.txt"
write.table(stability.combined.output, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/combined_output/", stabsel.combined.filename), sep="\t", row.names = T, col.names = NA )



################# Test overlap with stability selected output ###############

## Input dataset and stabsel output filename
stabsel_path <- "/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/combined_output/"
stabsel_file <- "combined_stabsel.txt"

stabsel_output_filepath <- paste0(stabsel_path,stabsel_file)
stabsel.output <- read.table(stabsel_output_filepath,sep="\t",head=T, row.names = 1, check.names = F)
dim(stabsel.output)




## massage column names of stability selected output for downstream processing
colnames(stabsel.output) <- c("gene","taxa")

## Identify stability selected associations
stabsel.output.filt <- stabsel.output[stabsel.output$taxa != "None",] 
dim(stabsel.output.filt)

## Identify unique stability selected genes
genes.stabsel <- unique(stabsel.output.filt$gene)
length(genes.stabsel)

## Identify unique stability selected predictors
taxa.stabsel <- unique(stabsel.output.filt$taxa)
length(taxa.stabsel)


## Overlap between combined.output.filt and stability selected assoc
overlap.assoc.stabsel <- merge(combined.output,stabsel.output.filt, by = c("gene","taxa"))
overlap.assoc.stabsel <- overlap.assoc.stabsel[order(overlap.assoc.stabsel$BH_global),]
dim(overlap.assoc.stabsel)

length(unique(overlap.assoc.stabsel$gene))

length(unique(overlap.assoc.stabsel$taxa))


## write to file before filtering by FDR cutoffs
filename <- "combined_output_MHT_stabsel.txt"
write.table(overlap.assoc.stabsel, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/combined_output/",filename), sep="\t", row.names = T, col.names = NA )

## How many stability selected assoc. also significant at FDR < 0.1
overlap.assoc.stabsel.BH.0.1 <- overlap.assoc.stabsel[overlap.assoc.stabsel$BH_global < 0.1,]
dim(overlap.assoc.stabsel.BH.0.1)

length(unique(overlap.assoc.stabsel.BH.0.1$gene))

length(unique(overlap.assoc.stabsel.BH.0.1$taxa)) 


############# Filter out non-taxa associations (e.g. gene-gender assoc.)######################

combined.output <- overlap.assoc.stabsel.BH.0.1

## Identify gender associations
select <- which(combined.output$taxa =="gender")
## Spit out gender associated output to a file for later inspection
gender.assoc <- combined.output[select,]; dim(gender.assoc)

length(unique(gender.assoc$gene))

## write to file
write.table(gender.assoc, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/combined_output/gene_gender_FDR_0.1.txt"), sep="\t", row.names = T, col.names = NA )

## filter out gender associations.
combined.output.filt <- combined.output[-select,]
dim(combined.output.filt)

## Additional gene-covariate filtering for IBD or IBS
if(dataset == "IBD/control/" || dataset == "IBS/control/"){
  ## set subtype label
  subtype_covariate <- "IBD_subtype" #IBD_subtype
  select <- which(combined.output.filt$taxa == subtype_covariate)
  condition.assoc <- combined.output.filt[select,]; dim(condition.assoc)
  
  length(unique(condition.assoc$gene)) 
  
  # write gene-condition associations to file
  # write.table(condition.assoc, file = paste0("./data/analysis/",dataset,"output_lasso_hdi_loocv_V2/combined_output/gene_IBS_subtype_FDR_0.1.txt"), sep="\t", row.names = T, col.names = NA )
  
  # #keep non-condition associations.
  if(length(select) != 0)
    combined.output.filt <- combined.output.filt[-select,]; dim(combined.output.filt)
}




## no. of unique genes
length(unique(combined.output.filt$gene))

## no. of unique taxa
length(unique(combined.output.filt$taxa))

## set datatset and filename and write to file
dataset <- "IBD/control/"
filename <- "IBD_control_gene_taxa_FDR_0.1.txt"
write.table(combined.output.filt, file = paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/lncRNA/Control/combined_output/",filename), sep="\t", row.names = T, col.names = NA )

################################################specific link
#----------------------------------------mRNA------------------------------------------------

############## Input datasets ##############

## IBD
## case
dataset <- "IBD/case"
filename <- "IBD_case_gene_taxa_FDR_0.1.txt"
filepath <- paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/CD/combined_output/", filename)
case <- data.frame(fread(filepath, sep="\t", head=T), row.names = 1)

dim(case)

## control 
dataset <- "IBD/control"
filename <- "IBD_control_gene_taxa_FDR_0.1.txt"
filepath <- paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/15/mRNA/Control/combined_output/", filename)

control <- data.frame(fread(filepath, sep="\t", head=T), row.names = 1)
dim(control)




############# Compute overlaps (FDR < 0.1) ###########

## overlapping associations
case_control_overlap_assoc <- merge(case, control, by = c("gene", "taxa"))
dim(case_control_overlap_assoc) 

## overlapping genes
case_control_overlap_genes <- intersect(case$gene, control$gene)
length(case_control_overlap_genes) 
case_control_overlap_genes

## overlapping taxa
case_control_overlap_taxa <- intersect(case$taxa, control$taxa)
length(case_control_overlap_taxa)
case_control_overlap_taxa

## write only case-specific genes to a file
case_specific_genes <- setdiff(case$gene, control$gene)
length(case_specific_genes)

dataset <- "disease/case"
filename <- "case_specific_genes_assoc_w_taxa.txt"
filepath <- paste0("/mnt/data/IBD/RNA_seq/data/Output/Others/16/mRNA/", filename)
write(case_specific_genes, file = filepath, sep="\n")
























