# Biomarker list
biomarkers <- c(
  "log_ab4240",
  "log_ptau217_42",
  "log_AMYLB42",
  "log_ORRESRAW",
  "log_TPP181",
  "log_TPP181_42"
)

# ROC curves for each marker
roc_df_all <- do.call(rbind, lapply(biomarkers, function(biomarker) {
  roc_result <- roc(all_data$group, all_data[[biomarker]])
  data.frame(
    specificity = roc_result$specificities,
    sensitivity = roc_result$sensitivities,
    auc = rep(auc(roc_result), length(roc_result$specificities)),
    feature = rep(biomarker, length(roc_result$specificities))
  )
}))

# ROC plot
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "black")
features_sorted <- biomarkers[order(biomarkers, decreasing = TRUE)]

roc_plot <- ggplot(roc_df_all, aes(x = 1 - specificity, y = sensitivity, color = feature)) +
  geom_step(size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = colors) +
  labs(title = "ROC curves", x = "1 - Specificity", y = "Sensitivity", color = "Feature") +
  ggprism::theme_prism(border = TRUE) +
  annotate(
    "text",
    x = 0.55,
    y = seq(0.05, 0.25, length.out = length(features_sorted)),
    label = paste0(features_sorted, " AUC: ", round(roc_df_all$auc[match(features_sorted, roc_df_all$feature)], 2)),
    size = 5, hjust = 0, color = rev(colors)
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1)

print(roc_plot)

# Youden cutoff
group_num <- ifelse(all_data$group %in% c("positive", "pos", "POS", "Positive", 1, "1"), 1, 0)

set.seed(123)
results_list <- list()

for (bm in biomarkers) {
  marker <- all_data[[bm]]
  idx <- !is.na(marker) & !is.na(group_num)
  marker_sub <- marker[idx]
  group_sub <- group_num[idx]

  roc_obj <- roc(response = group_sub, predictor = marker_sub)
  ci_auc <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000)

  youden_res <- coords(
    roc_obj, x = "best",
    ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"),
    transpose = FALSE
  )

  cutoff <- youden_res$threshold[1]
  cm <- coords(
    roc_obj, x = cutoff, input = "threshold",
    ret = c("tp", "tn", "fp", "fn"), transpose = FALSE
  )

  accuracy <- as.numeric((cm["tp"] + cm["tn"]) / (cm["tp"] + cm["tn"] + cm["fp"] + cm["fn"]))

  results_list[[bm]] <- data.frame(
    Biomarker = bm,
    AUC = as.numeric(roc_obj$auc),
    AUC_Lower = as.numeric(ci_auc[1]),
    AUC_Upper = as.numeric(ci_auc[3]),
    Cutoff = cutoff,
    Sensitivity = youden_res$sensitivity[1],
    Specificity = youden_res$specificity[1],
    PPV = youden_res$ppv[1],
    NPV = youden_res$npv[1],
    Accuracy = accuracy
  )
}

results_all <- do.call(rbind, results_list)
write_csv(results_all, "output/single_marker_roc_results.csv")
