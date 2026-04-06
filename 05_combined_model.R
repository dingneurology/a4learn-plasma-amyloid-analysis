# Combined model variables
base_vars <- c("AGEYR", "SEX", "APOE4_positive")
biomarker_vars <- c("log_ORRESRAW", "log_TPP181", "log_AMYLB42")

# Logistic model
final_model_combined <- glm(
  group ~ AGEYR + SEX + APOE4_positive + log_ORRESRAW + log_TPP181 + log_AMYLB42,
  data = all_data,
  family = binomial(link = "logit")
)

all_data$pred_combined <- predict(final_model_combined, type = "response")

# Combined ROC
roc_combined <- roc(all_data$group, all_data$pred_combined)
auc_combined <- auc(roc_combined)
ci_combined <- ci.auc(roc_combined)

# Compare with single markers
biomarkers <- c("log_ab4240", "log_ptau217_42", "log_AMYLB42",
                "log_ORRESRAW", "log_TPP181", "log_TPP181_42")

compare_results_combined <- biomarkers %>%
  map_df(function(bm) {
    roc_bm <- roc(all_data$group, all_data[[bm]])
    p_value <- roc.test(roc_combined, roc_bm, method = "delong")$p.value

    tibble(
      biomarker = bm,
      AUC_bm = as.numeric(auc(roc_bm)),
      AUC_combined = as.numeric(auc_combined),
      p_value = p_value
    )
  })

write_csv(compare_results_combined, "output/delong_compare_results.csv")

# Combined ROC summary
model_results <- data.frame(
  Biomarker = c(biomarkers, "Combined Model"),
  AUC = c(
    sapply(biomarkers, function(bm) auc(roc(all_data$group, all_data[[bm]]))),
    auc_combined
  ),
  AUC_lower = c(
    sapply(biomarkers, function(bm) ci.auc(roc(all_data$group, all_data[[bm]]))[1]),
    ci_combined[1]
  ),
  AUC_upper = c(
    sapply(biomarkers, function(bm) ci.auc(roc(all_data$group, all_data[[bm]]))[3]),
    ci_combined[3]
  ),
  Model_Type = c(rep("Biomarker", length(biomarkers)), "Biomarker + covariates")
)

write_csv(model_results, "output/auc_ci_results.csv")
