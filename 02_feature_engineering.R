# Reload merged data if needed
all_data <- demo %>%
  inner_join(amyloid_diagnose, by = c("BID", "SUBSTUDY")) %>%
  inner_join(plasma_wide, by = "BID") %>%
  inner_join(ptau217, by = "BID") %>%
  mutate(
    ab4240         = AMYLB42 / AMYLB40,
    ptau217_42     = ORRESRAW / AMYLB42,
    TPP181_42      = TPP181 / AMYLB42,
    log_AMYLB42    = log(AMYLB42),
    log_AMYLB40    = log(AMYLB40),
    log_ab4240     = log(ab4240),
    log_ORRESRAW   = log(ORRESRAW),
    log_ptau217_42 = log(ptau217_42),
    log_TPP181     = log(TPP181),
    log_TPP181_42  = log(TPP181_42),
    group          = factor(overall_score)
  )
write_csv(all_data, "output/all_data_ready.csv")
