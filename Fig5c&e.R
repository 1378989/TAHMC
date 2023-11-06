#########################Machine learning & ROC plot：


CompareModel <- function(training, validation, method){
  
  training <-df_train
  validation  <- df_vad
  
#Grid search for parameter tuning
  Grid <- list( 
                adaboost = expand.grid(nIter = c(50,100,150,200,250) ,method= c('Adaboost.M1','Real adaboost')),
                LogitBoost = expand.grid(nIter = c(11,16,21,26,31,36,41,46,51,56,61,66,71,76,81,86,91,96,101,111,121,151) ),
                gbm=expand.grid(interaction.depth = c(1,6,9),n.trees = (1:10) * 30,shrinkage = 0.1, n.minobsinnode = 20)
             
  )
  
  TuneLength =  list( 
                      adaboost = nrow(Grid[['adaboost']]),
                      LogitBoost =  nrow(Grid[['LogitBoost']]),
                      gbm= nrow(Grid[['gbm']])
                    
  )
  
  
  ##model training with different algorithms
  ls_model <- lapply(method,function(m){

     f = 5  # f folds resampling
      r = 10 # r repeats
      n = f*r
      
      # sets random seeds for parallel running for each single resampling f-folds and r-repeats cross-validation
      seeds <- vector(mode = "list", length = n + 1)
      #the number of tuning parameter
      for(i in 1:n) seeds[[i]] <- sample.int(n=1000, TuneLength[[m]])
      
      #for the last model
      seeds[[n+1]]<-sample.int(1000, 1)
      
      
      ctrl <- trainControl(method="repeatedcv",
                           number = f, ## 5-folds cv
                           summaryFunction=twoClassSummary,   # Use AUC to pick the best model
                           classProbs=TRUE,
                           repeats = r, ## 10-repeats cv,
                           seeds = seeds
      )
      
      
      
      model.tune <- train(group ~ .,
                          data = training,
                          method = m,
                          metric="ROC",
                          trControl=ctrl,
                          tuneGrid = Grid[[m]]
      )
    }
    print(m)
    return(model.tune)
  }
  )
  
  ##model validation
  auc <- lapply(ls_model,function(model.tune){
    if(c((class(model.tune) == 'predictor')[1])){
      pData <- data.frame(class = validation$group, sample = rownames(validation),row.names = rownames(validation))
      phenoData <- new("AnnotatedDataFrame",data=pData)
      HMNSig <- t(validation[,-1])
     HMNSig.test <- ExpressionSet(assayData=as.matrix(HMNSig),phenoData=phenoData)
      prediction <- predict(model.tune, HMNSig.test,"CD", ngenes=nrow(HMNSig), dist = "cor")
      roc <- roc(response  = prediction@prediction[,'class_membership'],
                 predictor = as.numeric(prediction@prediction[,'z'])
      )
      roc_result <- coords(roc, "best")
      auc <- data.frame(ROC=roc$auc, Sens = roc_result$sensitivity, Spec = roc_result$specificity)
      
    }else {
      prob <- predict(model.tune,validation[,-1],type = "prob")
      pre <- predict(model.tune,validation[,-1])
      test_set <- data.frame(obs = validation$group, CD = prob[,'CD'], Control = prob[,'Control'], pred=pre)
      auc <- twoClassSummary(test_set, lev = levels(test_set$obs))
    }
    
    return(auc)
  }) %>% do.call(rbind,.)
  
  rownames(auc) <- method
  
  res <- list()
  
  res[['model']] <- ls_model
  res[['auc']] <- auc
  
  
  return(res)
  
}

#parallel processing
cl <- makePSOCKcluster(60)
registerDoParallel(cl)

res <- CompareModel(training = training,
                    validation = validation,
                    method = c("adaboost",'LogitBoost',"gbm"))



# stopCluster(cl) 
print(res[['auc']]) ## 'nb' achieves best performance

class(res$model[[which.max(res$auc[,1])]]$finalModel) # best model
res$model[[which.max(res$auc[,1])]]$finalModel$tuneValue ##parameters

# save(res,file='results/Machine Learning/res.Rdata')




#----------------------------------------------------------------------------------------------------------

training <-df_train
validation  <- df_vad



##===loading data===###
library(ROCit)
library(caret)
library(dplyr)
library(survminer)
library(survival)
library(reshape2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ROCit)
library(ggsci)
library(cancerclass)



models <-  c('adaboost','LogitBoost',"gbm")

roc <- lapply(1:3,function(i){
    prob <- predict(res[['model']][[i]],validation[,-1],type = "prob") # use 'nb' model
    pre <- predict(res[['model']][[i]],validation[,-1]) # use 'nb' model
    test_set <- data.frame(obs = validation$group, CD = prob[,'CD'], Control = prob[,'Control'], pred=pre)
    roc <- rocit(score = test_set$CD,
                 class = test_set$obs,
                 negref = 'Control')

})

suppressWarnings(colset <- pal_lancet("lanonc", alpha = 0.7)(15))
plot(roc[[1]], col = colset[1], 
     legend = FALSE, YIndex = F)

for(i in 2:3){
  lines(roc[[i]]$TPR~roc[[i]]$FPR, 
        col = colset[i], lwd = 2)
}

legend("bottomright", 
       col = colset[1:6],
       paste(models, 'AUC', round(res$auc[,1],3)), 
       lwd = 2)





###===Figure 5e===###

#ROC plot

rocplot <- function(data){
  prob <- predict(res[['model']][[which.max(res$auc[,1])]],data[,-1],type = "prob") # use 'nb' model
  pre <- predict(res[['model']][[which.max(res$auc[,1])]],data[,-1]) # use 'nb' model
  test_set <- data.frame(obs = data$group, CD= prob[,'CD'], Control= prob[,'Control'], pred=pre)
  roc <- ROCit::rocit(score = test_set$CD,
                      class = test_set$obs,
                      negref = 'Control')
  plot(roc,legend=F,YIndex = F)
  title(deparse(substitute(data)))
  text(x=0.8,y=0.6,labels = paste0("AUC: ",round(roc$AUC,2), " (95%CI: ",round(as.numeric(ciAUC(roc)[5]),2),'-',round(as.numeric(ciAUC(roc)[6]),2),")"))
  
}


rocplot(validation) # validation
#rocplot(df_test) # testing
rocplot(df_test3) 
stopCluster(cl)
registerDoSEQ()


pdf('/mnt/data/IBD/RNA_seq/data/Output/Others/22/model.pdf',width = 10,height = 10)
rocplot(df_test) 
dev.off()

