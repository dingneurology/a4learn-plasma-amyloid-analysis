# Cutoff map
biomarker_cutoffs <- results_all %>%
  select(Biomarker, Cutoff) %>%
  deframe()

# PET vs biomarker plot
plot_biomarker_vs_pmod <- function(bm_name) {
  bm_cut <- biomarker_cutoffs[[bm_name]]

  plot_data <- all_data %>%
    transmute(
      pmod_suvr = pmod_suvr,
      biomarker = .data[[bm_name]],
      group_fac = case_when(
        group %in% c("negative", 0) ~ "Aβ-",
        group %in% c("positive", 1) ~ "Aβ+",
        TRUE ~ as.character(group)
      )
    ) %>%
    mutate(group_fac = factor(group_fac, levels = c("Aβ-", "Aβ+")))

  ggplot(plot_data, aes(x = pmod_suvr, y = biomarker)) +
    geom_point(aes(color = group_fac), alpha = 0.7, size = 2.2) +
    geom_smooth(method = "lm", se = FALSE, color = "black", size = 1.1) +
    geom_vline(xintercept = 1.15, linetype = "dashed", color = "#8A8A8A") +
    geom_hline(yintercept = bm_cut, linetype = "dashed", color = "#8A8A8A") +
    scale_color_manual(values = c("Aβ-" = "#4D4D4D", "Aβ+" = "#D55E00")) +
    labs(x = "Aβ-PET SUVR", y = bm_name) +
    theme_minimal(base_size = 14)
}

# Distribution plot
all_data2 <- all_data %>%
  mutate(
    ab_status = ifelse(group %in% c("positive", 1, "Aβ+"), 1, 0),
    group_fac = ifelse(ab_status == 1, "Aβ+", "Aβ-"),
    group_fac = factor(group_fac, levels = c("Aβ-", "Aβ+"))
  )

plot_biomarker_dist <- function(biomarker_name) {
  biomarker_cutoff <- biomarker_cutoffs[[biomarker_name]]

  plot_data <- all_data2 %>%
    transmute(
      biomarker = .data[[biomarker_name]],
      ab_status,
      group_fac
    )

  glm_fit <- glm(ab_status ~ biomarker, data = plot_data, family = binomial(link = "logit"))

  newdat <- tibble(
    biomarker = seq(min(plot_data$biomarker), max(plot_data$biomarker), length.out = 300)
  ) %>%
    mutate(prob_ab_pos = predict(glm_fit, newdata = ., type = "response"))

  ggplot(plot_data, aes(x = biomarker)) +
    geom_histogram(aes(fill = group_fac, y = ..count.. / max(..count..)),
                   bins = 40, alpha = 0.45, color = "white", position = "identity") +
    geom_density(aes(color = group_fac, y = ..scaled..), size = 1.1, alpha = 0.8) +
    geom_line(data = newdat, aes(x = biomarker, y = prob_ab_pos), color = "black", size = 1.2) +
    geom_vline(xintercept = biomarker_cutoff, color = "grey40", linetype = "dashed") +
    scale_fill_manual(values = c("Aβ-" = "#4D4D4D", "Aβ+" = "#D55E00")) +
    scale_color_manual(values = c("Aβ-" = "#4D4D4D", "Aβ+" = "#D55E00")) +
    theme_minimal(base_size = 14)
}

# Forest plot
forest_lm_biomarker <- function(data, biomarker, biomarker_label,
                                outcome = "pmod_suvr", age_var = "AGEYR", sex_var = "SEX") {

  fit <- lm(as.formula(paste0(outcome, " ~ ", biomarker, " + ", age_var, " + ", sex_var)), data = data)

  model_summary <- tidy(fit, conf.int = TRUE, conf.level = 0.95)

  coef_df <- model_summary |>
    filter(term != "(Intercept)") |>
    mutate(
      var_label = case_when(
        term == biomarker ~ biomarker_label,
        term == age_var ~ "Age (years)",
        grepl(paste0("^", sex_var), term) ~ "Sex (reference: male)",
        TRUE ~ term
      ),
      p_label = if_else(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
      beta_ci_label = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
    )

  tabletext <- rbind(
    c("Variable", "Estimate (95% CI)", "P value"),
    cbind(coef_df$var_label, coef_df$beta_ci_label, coef_df$p_label)
  )

  forestplot(
    labeltext = tabletext,
    mean = c(NA, coef_df$estimate),
    lower = c(NA, coef_df$conf.low),
    upper = c(NA, coef_df$conf.high),
    zero = 0,
    boxsize = 0.2,
    lineheight = unit(6, "mm"),
    colgap = unit(3, "mm"),
    lwd.zero = 1.2,
    lwd.ci = 2,
    col = fpColors(box = "#458B00", lines = "black", zero = "#7AC5CD"),
    xlab = "Estimate (95% CI)"
  )
}

# Example
forest_lm_biomarker(
  data = all_data,
  biomarker = "ptau217_42",
  biomarker_label = "Plasma p-tau217/Aβ42"
)
