eval.model <- function(model, newData, newLabels, ml_dir, keep_all=FALSE, method="", feature="", seed) {
  
  # Evaluate predictions for model on new data
  
  feature<-rownames(model$glmnet.fit$beta)
  
  if (keep_all){
    # check that all features are used
    if (cv.glmfit$nzero[cv.glmfit$index["min",]] == max(cv.glmfit$nzero)){
      # if so, use best value of lambda
      lambda_use = cv.glmfit$lambda.min
    }else{
      # if not, use best model with all features
      if (cv.glmfit$name == "AUC") {
        idx_nzero_use = which.max(cv.glmfit$cvm[cv.glmfit$nzero == max(cv.glmfit$nzero)])
        lambda_use = cv.glmfit$lambda[cv.glmfit$nzero == max(cv.glmfit$nzero)][idx_nzero_use]
      }else{
        idx_nzero_use = which.min(cv.glmfit$cvm[cv.glmfit$nzero == max(cv.glmfit$nzero)])
        lambda_use = cv.glmfit$lambda[cv.glmfit$nzero == max(cv.glmfit$nzero)][idx_nzero_use]
      }
    }
  }else{
    lambda_use = cv.glmfit$lambda.min
  }
  
  testPredictions<-predict(model,newx=t(newData[feature,]),s=lambda_use,type="response")
  
  # trainError <- model$cvm[which(model$lambda==model$lambda.min)]
  
  # Compute performance statistics on the training data. 
  
  # confMat <- table(ifelse(trainPredictions < 0.5, "case", "control"), trainLabels$Sample_Group)
  # if (dim(confMat)[1] == 1) {
  #   if (rownames(confMat) == "control") {
  #     confMat <- rbind(c(0, 0), confMat)
  #   } else {
  #     confMat <- rbind(confMat, c(0, 0))
  #   }
  #   rownames(confMat) <- c("case", "control")
  # }
  # SN <- confMat[1]/(confMat[1] + confMat[2])
  # SP <- confMat[4]/(confMat[3] + confMat[4])
  # trainError <- (confMat[2] + confMat[3])/sum(confMat)
  
  # Make predictions on the test data.
  # Test data is denoted confMat, 'predictions'
  
  # testPredictions <- predict(model, newx = t(test[feature,]), s = model$lambda.min, type = "response")
  
  confMat <- table(ifelse(testPredictions < 0.5, "case", "control"), newLabels$Sample_Group)
  
  if (dim(confMat)[1] == 1) {
    if (rownames(confMat) == "control") {
      confMat <- rbind(c(0, 0), confMat)
    } else {
      confMat <- rbind(confMat, c(0, 0))
    }
    rownames(confMat) <- c("case", "control")
  }
  
  SN <- confMat[1]/(confMat[1] + confMat[2])
  SP <- confMat[4]/(confMat[3] + confMat[4])
  testError <- (confMat[2] + confMat[3])/sum(confMat)
  
  # Compute AUC curves and ROC data for test and lockbox using ROCR package
  
  predTest <- prediction(testPredictions,newLabels$Sample_Group)
  rocTest <- performance(predTest,"tpr","fpr")
  aucTest <- performance(predTest,"auc")
  precisionRecall<-performance(predTest,"prec","rec")
  
  # Computing AUC curves and ROC data using the pROC package
  
  outcomes<-data.frame(cbind(testPredictions,newLabels$Sample_Group))
  colnames(outcomes)[1]<-"predicted"
  colnames(outcomes)[2]<-"actual"
  if(class(try({pROC::roc(outcomes$actual,as.numeric(outcomes$predicted),smooth=T)})) == "try-error") {
    ROC2<-pROC::roc(outcomes$actual,as.numeric(outcomes$predicted),smooth=F)
  }else{
    ROC2<-pROC::roc(outcomes$actual,as.numeric(outcomes$predicted),smooth=T)
  }
  # ROC2<-pROC::roc(outcomes$actual,as.numeric(outcomes$predicted),smooth=T)
  AUCpROC<-pROC::auc(ROC2)
  
  save(ROC2,file=paste0(ml_dir, "/resultsfsvab/pROC_ROC_",seed,".R"))
  save(AUCpROC,file=paste0(ml_dir,"/resultsfsvab/pROC_AUC_",seed,".R"))
  save(outcomes,file=paste0(ml_dir,"/resultsfsvab/predictedOutcomes_",seed,".R"))
  
  print(ROC2)
  
  # Save out the glmnet model for future analysis
  
  save(model,file=paste0(ml_dir,"/resultsfsvab/CVglmfit_",seed,".R"))
  
  # Rerun the model selecting features incrementally
  
  coefs <- as.array(coef(model, s = lambda_use))
  featuresIncluded <- names(which(rowSums(coefs) != 0))[-1]
  
  
  #lengthFeaturesIncluded=nrow(featuresIncluded)
  #lengthFeaturesIncluded=lengthFeaturesIncluded-1
  
  #print(paste0("Total number of features this cycle: ",lengthFeaturesIncluded))
  
  #for (seed10 in 2:lengthFeaturesIncluded) {
  
  #incremental<-incremental.performance(train,test,trainLabels,newLabels,method="exact",feature=feature,seed10)
  
  #save(incremental, file = paste0("./resultsfsvab/Incremental_unsupervised","_",seed,"_",seed10,".R", sep = ""))
  
  #}
  
  # Return the model and performance statistics
  
  return(list(testError = testError, confMat = confMat, features = featuresIncluded, testPred = testPredictions, testLabels = newLabels[,1],aucTest = aucTest, rocTest = rocTest))
  
}