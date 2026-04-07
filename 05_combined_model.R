# Combined model variables
base_vars <- c("AGEYR", "SEX", "APOE4_positive")
biomarker_vars <- c("log_ORRESRAW", "log_TPP181", "log_AMYLB42")

# Outcome as 0/1
all_data <- all_data %>%
  mutate(
    ab_status = ifelse(group %in% c("positive", 1, "Aβ+"), 1, 0)
  )

# Biomarker list for single-marker models
biomarkers <- c(
  "log_ab4240",
  "log_ptau217_42",
  "log_AMYLB42",
  "log_ORRESRAW",
  "log_TPP181",
  "log_TPP181_42"
)

# ---------------------------
# Combined logistic model
# ---------------------------
final_model_combined <- glm(
  ab_status ~ AGEYR + SEX + APOE4_positive + log_ORRESRAW + log_TPP181 + log_AMYLB42,
  data = all_data,
  family = binomial(link = "logit")
)

all_data$pred_combined <- predict(final_model_combined, type = "response")

roc_combined <- roc(all_data$ab_status, all_data$pred_combined, quiet = TRUE)
auc_combined <- as.numeric(auc(roc_combined))
ci_combined <- ci.auc(roc_combined)

# ---------------------------
# DeLong test:
# combined model vs covariate-adjusted single-biomarker models
# ---------------------------
compare_results_combined <- map_df(biomarkers, function(bm) {

  form_bm <- as.formula(
    paste("ab_status ~ AGEYR + SEX + APOE4_positive +", bm)
  )

  fit_bm <- glm(
    form_bm,
    data = all_data,
    family = binomial(link = "logit")
  )

  pred_bm <- predict(fit_bm, type = "response")
  roc_bm  <- roc(all_data$ab_status, pred_bm, quiet = TRUE)

  delong_test <- roc.test(roc_combined, roc_bm, method = "delong")

  tibble(
    biomarker = bm,
    model_bm = paste0(bm, " + AGEYR + SEX + APOE4_positive"),
    AUC_bm = as.numeric(auc(roc_bm)),
    AUC_combined = auc_combined,
    p_value = delong_test$p.value
  )
})

compare_results_combined <- compare_results_combined %>%
  mutate(
    AUC_bm = round(AUC_bm, 3),
    AUC_combined = round(AUC_combined, 3),
    p_value_fmt = case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  )

write_csv(compare_results_combined, "output/delong_compare_results.csv")

# ---------------------------
# AUC summary:
# covariate-adjusted single-biomarker models + combined model
# ---------------------------
single_model_results <- map_df(biomarkers, function(bm) {

  form_bm <- as.formula(
    paste("ab_status ~ AGEYR + SEX + APOE4_positive +", bm)
  )

  fit_bm <- glm(
    form_bm,
    data = all_data,
    family = binomial(link = "logit")
  )

  pred_bm <- predict(fit_bm, type = "response")
  roc_bm  <- roc(all_data$ab_status, pred_bm, quiet = TRUE)
  ci_bm   <- ci.auc(roc_bm)

  tibble(
    Biomarker = bm,
    AUC = as.numeric(auc(roc_bm)),
    AUC_lower = as.numeric(ci_bm[1]),
    AUC_upper = as.numeric(ci_bm[3]),
    Model_Type = "Biomarker + covariates"
  )
})

combined_result <- tibble(
  Biomarker = "Combined Model",
  AUC = auc_combined,
  AUC_lower = as.numeric(ci_combined[1]),
  AUC_upper = as.numeric(ci_combined[3]),
  Model_Type = "Combined model"
)

model_results <- bind_rows(single_model_results, combined_result) %>%
  mutate(
    AUC = round(AUC, 3),
    AUC_lower = round(AUC_lower, 3),
    AUC_upper = round(AUC_upper, 3)
  )

write_csv(model_results, "output/auc_ci_results.csv")
