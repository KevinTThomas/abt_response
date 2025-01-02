ml_dir = "ML_results/ML_ACCESS"
runIndexLabel

num_features = c()
for (i in 1:40) {
  load(paste0(ml_dir,"/results/TrueLabels_", runIndexLabel, "_", i,".R"))
  num_features[i]=length(devResults$features)
}
# min(num_features)
# max(num_features)
# median(num_features)

# tibble::enframe(num_features, "seed", "num_features") %>%
#   ggplot(aes(x = num_features)) +
#   geom_histogram(fill = "black", binwidth = 2) +
#   theme_bw()

# svglite(file = "figures/ml_T4_num_feats.svg", width = 7, height = 7)
# read.csv(paste0(ml_dir, "/topFeatures.csv")) %>%
#   head(20) %>%
#   mutate(feature = factor(feature, levels = feature)) %>%
#   ggplot(aes(x = feature, y = counts)) +
#   geom_col(aes(fill = counts>35), show.legend = FALSE) +
#   scale_fill_manual(values = c("TRUE" = "blue", "FALSE", "grey50")) +
#   labs(x = "CpG probe", y = "Number of models used") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# dev.off()
# 
# png(filename = paste0(ml_dir, "/num_feats.png"), width = 720, height = 720)
# read.csv(paste0(ml_dir, "/topFeatures.csv")) %>%
#   head(20) %>%
#   mutate(feature = factor(feature, levels = feature)) %>%
#   ggplot(aes(x = feature, y = counts)) +
#   geom_col() +
#   labs(x = "CpG probe", y = "Number of models used") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# dev.off()
# 
# ml_dirs = c("results/ML_Tcell", "results/ML_Bcell", "results/ML_Mono")
# 
# for (i in seq_along(ml_dirs)){
#   assign(
#     x = paste0("gg_",i),
#     value = read.csv(paste0(ml_dirs[i], "/topFeatures.csv")) %>%
#       head(20) %>%
#       mutate(feature = factor(feature, levels = feature)) %>%
#       ggplot(aes(x = feature, y = counts)) +
#       geom_col(aes(fill = counts>35), show.legend = FALSE) +
#       scale_fill_manual(values = c("TRUE" = "blue", "FALSE", "grey50")) +
#       labs(x = "CpG probe", y = "Number of models used") +
#       theme_bw() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
#   )
# }
# 
# svglite(file = "figures/ml_celltypes_num_feats.svg", width = 7, height = 10)
# cowplot::plot_grid(
#   gg_1 + labs(title = "T Cells"),
#   gg_2 + labs(title = "B Cells"),
#   gg_3 + labs(title = "Monocytes"),
#   ncol = 1
# )
# dev.off()


load(file = "data/betaMatrix.R")

all_coef = matrix(nrow = nrow(dataClean))
for (i in 1:40) {
  load(paste(ml_dir,"/resultsfsvab/CVglmfit_",i,".R", sep = ""))
  idx = read.csv(paste0(ml_dir, "/indices/developmentIndicesUnsupervised_", i, ".csv"))$developmentIndices
  betas = coef(cv.glmfit, s = cv.glmfit$lambda.min)
  feature_sds = dataClean[,idx] %>%
    as.matrix() %>%
    matrixStats::rowSds()
  norm_coef = betas[-1,]*feature_sds
  all_coef = cbind(all_coef, norm_coef)
}
all_coef = all_coef[,-1]
colnames(all_coef) <- 1:40

topFeatCoef = all_coef[read.csv(paste0(ml_dir, "/topFeatures.csv"))$feature,]
mean_and_sd = function(x) {
  c(mean(abs(x)), sd(x))
}

# svglite(file = "figures/ml_all_feat_imp.svg", width = 7, height = 7)
# topFeatCoef %>%
#   apply(1, mean_and_sd) %>%
#   t() %>%
#   `colnames<-`(c("mean", "sd")) %>%
#   as.data.frame() %>%
#   rownames_to_column("feature") %>%
#   mutate(feature = factor(feature, levels = feature)) %>%
#   head(20) %>%
#   ggplot(aes(x = feature, y = mean)) +
#   geom_col() +
#   geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.3) +
#   labs(x = "CpG probe", y = "Standardized variable importance") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
# dev.off()

png(filename = paste0(ml_dir, "num_feats.png"), width = 7, height = 7, units = "in")
topFeatCoef %>%
  apply(1, mean_and_sd) %>%
  t() %>%
  `colnames<-`(c("mean", "sd")) %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>%
  mutate(feature = factor(feature, levels = feature)) %>%
  head(20) %>%
  ggplot(aes(x = feature, y = mean)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.3) +
  labs(x = "CpG probe", y = "Standardized variable importance") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
dev.off()