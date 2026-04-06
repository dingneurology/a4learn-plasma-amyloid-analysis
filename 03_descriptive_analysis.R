# Theme and palette
theme_nature <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      legend.position = "none"
    )
}

pal_nc <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52")

# Half violin plot
plot_half_violin <- function(data, y_var, y_lab = NULL) {
  ggplot(data, aes(x = group, y = .data[[y_var]], fill = group)) +
    geom_half_violin(
      side = "l", scale = "width", width = 0.9,
      position = position_nudge(x = -0.15),
      color = NA, alpha = 0.7
    ) +
    geom_boxplot(
      width = 0.18, position = position_nudge(x = 0.05),
      outlier.shape = NA, fill = "white", color = "black", linewidth = 0.6
    ) +
    geom_jitter(
      aes(color = group), width = 0.08, size = 1.5, alpha = 0.35, stroke = 0
    ) +
    scale_fill_manual(values = pal_nc) +
    scale_color_manual(values = pal_nc) +
    labs(x = NULL, y = y_lab) +
    theme_nature()
}

# Example plots
p_ab4240 <- plot_half_violin(all_data, "log_ab4240", "log(Aβ42/Aβ40)")
p_ptau_ratio <- plot_half_violin(all_data, "log_ptau217_42", "log(p-tau217 / Aβ42)")
p_ab42 <- plot_half_violin(all_data, "log_AMYLB42", "log(Plasma Aβ42)")
p_ptau <- plot_half_violin(all_data, "log_ORRESRAW", "log(Plasma p-tau217)")
p_tpp181 <- plot_half_violin(all_data, "log_TPP181", "log(TPP181)")
p_tpp181_42 <- plot_half_violin(all_data, "log_TPP181_42", "log(TPP181 / Aβ42)")

# Unadjusted comparison
all_data$APOE4_positive <- factor(all_data$APOE4_positive)

perform_unadjusted_analysis <- function(data) {
  biomarkers <- c("log_ab4240", "log_ptau217_42", "log_AMYLB42",
                  "log_ORRESRAW", "log_TPP181", "log_TPP181_42")
  results <- list()

  for (biomarker in biomarkers) {
    summary_stats <- data %>%
      group_by(group) %>%
      summarise(
        mean_val = mean(get(biomarker), na.rm = TRUE),
        sd_val = sd(get(biomarker), na.rm = TRUE),
        .groups = "drop"
      )

    t_test <- t.test(as.formula(paste(biomarker, "~ group")), data = data, var.equal = TRUE)
    cohen_d <- cohens_d(as.formula(paste(biomarker, "~ group")), data = data)

    lm_model <- lm(as.formula(paste(biomarker, "~ group + AGEYR + SEX + APOE4_positive")), data = data)
    p_value_adjusted <- coef(summary(lm_model))[2, 4]

    results[[biomarker]] <- list(
      summary_stats = summary_stats,
      t_test = tidy(t_test),
      cohen_d = cohen_d$Cohens_d,
      p_value_adjusted = p_value_adjusted
    )
  }

  results
}

unadjusted_results <- perform_unadjusted_analysis(all_data)

# Table 1 helper
pvalue <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times = sapply(x, length)))
  if (is.numeric(y)) {
    p <- t.test(y ~ g)$p.value
  } else {
    p <- chisq.test(table(y, g))$p.value
  }
  c("", sub("<", "&lt;", format.pval(p, digits = 3, eps = 0.0001)))
}

table1(
  ~ AGEYR + factor(SEX) + EDCCNTU + factor(RACE) + APOE4_positive +
    ab4240 + ptau217_42 + TPP181 + AMYLB42 + ORRESRAW + TPP181_42 | group,
  data = all_data,
  overall = FALSE,
  extra.col = list(`P-value` = pvalue),
  test = "fisher"
)
