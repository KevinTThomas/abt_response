plot_auc_curves <- function (
  ml_dir,
  runIndexLabel,
  seed,
  trainTest = c("train", "test"),
  colIndiv = "lightblue",
  colSumm = "blue"
  # fileName = "ml_access_train",
  # fig.width = 7,
  # fig.height = 7
) {
  # Stop if trainTest isn't set correctly
  if (length(trainTest) != 1 || !(trainTest %in% c("train", "test"))) {
    stop("Variable `trainTest` must be set to either 'train' or 'test'.")
  }
  
  # Concatenate all ROC data from labels
  concatData <- data.frame()
  label = ifelse(is.null(runIndexLabel), "_results", paste0("_", runIndexLabel))
  devResults <- get(load(file = paste0(ml_dir, "/results/TrueLabels",label,"_1.R")))
  ## If test labels was given as a dataframe, convert to character vector and get ROC stats
  if ("data.frame" %in% class(devResults[[paste0(trainTest,"Labels")]])) {
    concatData <- roc(factor(devResults[[paste0(trainTest,"Labels")]][[1]], levels = c("control", "case")), devResults[[paste0(trainTest,"Pred")]][,1], quiet=TRUE) %>%
      (function(x) cbind(coords(x), aucROC = as.numeric(x$auc))) %>%
      mutate(seed = 1, idx = 1:n()) %>%
      rbind(concatData, .)
    if (seed >=2) {
      for (seed in 2:seed) {
        devResults <- get(load(file =paste0(ml_dir, "/results/TrueLabels",label,"_",seed,".R")))
        concatData <- roc(factor(devResults[[paste0(trainTest,"Labels")]][[1]], levels = c("control", "case")), devResults[[paste0(trainTest,"Pred")]][,1], quiet=TRUE) %>%
          (function(x) cbind(coords(x), aucROC = as.numeric(x$auc))) %>%
          mutate(seed = seed, idx = 1:n()) %>%
          rbind(concatData, .)
      }
    }
    ## Otherwise, get ROC stats directly from vector labels
  }else{
    concatData <- roc(factor(devResults[[paste0(trainTest,"Labels")]], levels = c("control", "case")), devResults[[paste0(trainTest,"Pred")]][,1], quiet=TRUE) %>%
      (function(x) cbind(coords(x), aucROC = as.numeric(x$auc))) %>%
      mutate(seed = 1, idx = 1:n()) %>%
      rbind(concatData, .)
    if (seed >= 2) {
      for (seed in 2:seed) {
        devResults <- get(load(file = paste0(ml_dir, "/results/TrueLabels",label,"_",seed,".R")))
        concatData <- roc(factor(devResults[[paste0(trainTest,"Labels")]], levels = c("control", "case")), devResults[[paste0(trainTest,"Pred")]][,1], quiet=TRUE) %>%
          (function(x) cbind(coords(x), aucROC = as.numeric(x$auc))) %>%
          mutate(seed = seed, idx = 1:n()) %>%
          rbind(concatData, .)
      }
    }
  }
  
  # Calculate mean and sd of AUC and store in dataframe
  AUC_summary <- concatData %>%
    dplyr::select(seed, aucROC) %>%
    unique() %>%
    summarize(aucROC_avg = mean(aucROC), aucROC_sd = sd(aucROC))
  
  # Calculate average sensitivity and specificity coordinates of all curves
  avgData<-concatData %>%
    mutate(
      bin = cut(threshold, breaks = c(-Inf,seq(0,1,by = 0.1),Inf)),
      bin = ifelse(is.na(bin), 1, bin)
    ) %>%
    group_by(bin) %>%
    summarize(sensitivity = mean(sensitivity), specificity = mean(specificity)) %>%
    arrange(-specificity, sensitivity)
  
  # Plot final data in graph
  gg = concatData %>%
    ggplot(aes(x = 1-specificity, y = sensitivity, group = as.factor(seed))) +
    geom_path(col = colIndiv, alpha = 0.3) +
    geom_line(data = avgData, col = colSumm, linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, col = "lightgrey") +
    coord_fixed() +
    labs(
      title = paste0(toupper(substr(trainTest, 1, 1)), substr(trainTest, 2, nchar(trainTest)), "ing"),
      subtitle = paste0("AUC = ", signif(AUC_summary$aucROC_avg,digits = 3), " (sd=", signif(AUC_summary$aucROC_sd,digits = 3), ")")
    ) +
    theme_bw()
  
  # Save an svg file, DOES NOT WORK
  # if (!dir.exists(paste0(ml_dir, "/figures"))) {dir.create(path = paste0(ml_dir, "/figures"))}
  # svglite(file = paste0(ml_dir, "/figures/", fileName, ".svg"), width = fig.width, height = fig.height)
  # cowplot::plot_grid(gg)
  # dev.off()
  
  return(gg)
}


