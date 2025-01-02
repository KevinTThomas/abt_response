methLab.evaluateClassifier<-function(data,pheno,prop=0.7,ml_dir,n_cores,nfolds,keep_all=FALSE,resultsLabel="",feature,identificationLabel,seeds,force=FALSE){
  # Create directories if they don't exist
  for (d in paste0(ml_dir, c("/indices", "/results"))) {
    if (!dir.exists(d)) {dir.create(path = d, recursive = TRUE)}
  }
  
  # Check if "results" directory has been filled.
  # If `force` is FALSE, don't re-do the expensive calculation
  if (length(list.files(paste0(ml_dir, "/results")))>0 & force == FALSE) {
    message("Directory contains a populated 'results' folder. To force a recalculation, set `force` to `TRUE`.")
    return(invisible())
  }
  
  # use parallel backend
  message("Setting up socket cluster with ", n_cores, " cores.")
  cluster = parallel::makeCluster(
    n_cores,
    type = "PSOCK"
  )
  doParallel::registerDoParallel(cl = cluster)
  doRNG::registerDoRNG(seed = 12345, once = TRUE)
  
  foreach(seed = 1:seeds, .packages = c("rio", "glmnet", "ROCR", "pROC")) %dorng% {
    #attempt memory management
    gc(verbose = FALSE)
    # subset cases and controls, randomize each, split into development and validation
    permutedIndicesCases <- which(pheno$Sample_Group=="case")
    permutedIndicesCases <- sample(permutedIndicesCases)
    
    permutedIndicesControls <- which(pheno$Sample_Group=="control")
    permutedIndicesControls <- sample(permutedIndicesControls)
    
    developmentCases <- permutedIndicesCases[1:round(length(permutedIndicesCases)*prop)]
    validationCases <- permutedIndicesCases[(length(developmentCases)+1):length(permutedIndicesCases)]
    
    developmentControls <- permutedIndicesControls[1:round(length(permutedIndicesControls)*prop)]
    validationControls <- permutedIndicesControls[(length(developmentControls)+1):length(permutedIndicesControls)]
    
    developmentIndices <- c(developmentCases,developmentControls)
    validationIndices <- c(validationCases,validationControls)
    
    export(data.frame(developmentIndices),file=paste0(ml_dir,"/indices/developmentIndices",identificationLabel,"_",seed,".csv"))
    export(data.frame(validationIndices),file=paste0(ml_dir, "/indices/validationIndices",identificationLabel,"_",seed,".csv"))
    
    # subset phenotype data for development and validation sets
    
    developmentPheno <- pheno[developmentIndices,"Sample_Group",drop=FALSE]
    validationPheno <- pheno[validationIndices,"Sample_Group",drop=FALSE]
    
    developmentData <- data[,developmentIndices]
    validationData <- data[,validationIndices]
    
    # print(c(identificationLabel,seed))
    
    # pass data, pheno, and feature to classification function
    ## not
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
  message("Closing socket cluster.")
  parallel::stopCluster(cluster)
}