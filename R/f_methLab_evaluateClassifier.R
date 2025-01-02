methLab.evaluateClassifier<-function(data,pheno,prop=0.7,ml_dir,n_cores,nfolds,keep_all=FALSE,resultsLabel="",feature,identificationLabel,seeds){
  # use parallel backend
  message("Setting up socket cluster with ", n_cores, " cores")
  cluster = parallel::makeCluster(
    n_cores,
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = cluster)
  doRNG::registerDoRNG(seed = 12345, once = TRUE)
  
  # Create directories if they don't exist
  for (d in paste0(ml_dir, c("/indices", "/results"))) {
    dir.create(d, showWarnings = FALSE)
  }
  
  foreach(seed = 1:seeds, .packages = c("rio", "glmnet", "ROCR", "pROC")) %dorng% {
    # subset cases and controls, randomize each, split into dev and validation and recombine
    permutedIndicesCases <- which(pheno$Sample_Group=="case")
    permutedIndicesCases <- sample(permutedIndicesCases)
    
    permutedIndicesControls <- which(pheno$Sample_Group=="control")
    permutedIndicesControls <- sample(permutedIndicesControls)
    
    # below sets the split ratio, e.g. 0.7 for 70/30, 0.8 for 80/20, etc
    
    developmentCases1 <- permutedIndicesCases[1:round(length(permutedIndicesCases)*prop)]
    validationCases <- permutedIndicesCases[(length(developmentCases1)+1):length(permutedIndicesCases)]
    
    developmentControls1 <- permutedIndicesControls[1:round(length(permutedIndicesControls)*prop)]
    validationControls <- permutedIndicesControls[(length(developmentControls1)+1):length(permutedIndicesControls)]
    
    developmentIndices <- c(developmentCases1,developmentControls1)
    validationIndices <- c(validationCases,validationControls)
    
    export(data.frame(developmentIndices),file=paste0(ml_dir,"/indices/developmentIndices",identificationLabel,"_",seed,".csv"))
    export(data.frame(validationIndices),file=paste0(ml_dir, "/indices/validationIndices",identificationLabel,"_",seed,".csv"))
    
    # subset pheno and data for development and validation sets
    
    developmentPheno <- pheno[developmentIndices,]
    validationPheno <- pheno[validationIndices,]
    
    developmentData <- data[,developmentIndices]
    validationData <- data[,validationIndices]
    
    print(c(identificationLabel,seed))
    
    # pass data, pheno, and feature to classification function
    source("../R/f_classify_feature.R")
    devResults<-classify.feature(
      train=developmentData,
      test = validationData,
      trainLabels = developmentPheno,
      testLabels = validationPheno,
      ml_dir = ml_dir,
      nfolds = nfolds,
      keep_all = keep_all,
      method="exact",
      feature=feature,
      seed=seed
    )
    
    save(devResults, file = paste0(ml_dir, "/results/TrueLabels_",resultsLabel,"_",seed,".R"))
  }
  parallel::stopCluster(cluster)
}