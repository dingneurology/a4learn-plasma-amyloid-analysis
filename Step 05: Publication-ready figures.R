############################################################
## Step 05: Publication-ready figures
##
## Project:
## ADNI-ASG: Autophagy–Synapse–Glia inflammatory gate
##
## Main figures:
## Figure 1. Study design and ASG framework
## Figure 2. Primary Aβ × ASG interaction on CSF p-tau / t-tau
## Figure 3. Four-group biochemical tau burden
## Figure 4. Secondary future tau PET validation
## Figure 5. Forest plot summary of key model terms
##
## Inputs:
## 03_merged_outputs/
##   Step03_analysis_ready_biochemical_fullcov.csv
##   Step03B_analysis_ready_future_tauPET_7y_fullcov.csv
##   Step03B_analysis_ready_future_tauPET_10y_fullcov.csv
##   Step03B_analysis_ready_future_tauPET_any_fullcov.csv
##
## 04_model_outputs/
##   Step04A_primary_interaction_prediction_grid.csv
##   Step04B_four_group_estimated_marginal_means.csv
##   Step04B_four_group_pairwise_contrasts_tukey.csv
##   Step04B_four_group_biochemical_descriptives.csv
##   Step04C_future_tauPET_7y_models_prediction_grid.csv
##   Step04C_future_tauPET_10y_models_prediction_grid.csv
##   Step04C_future_tauPET_any_models_prediction_grid.csv
##   Step04_manuscript_key_model_terms.csv
##
## Outputs:
## 05_figures/
##   PDF + PNG publication figures
##   Source plotting tables
############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(janitor)
  library(stringr)
  library(scales)
})

if (requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages(library(patchwork))
}

############################################################
## 0. Paths
############################################################

ROOT_DIR <- "/Users/dyr/Desktop/ADNI_ASG"

MERGED_DIR <- file.path(ROOT_DIR, "03_merged_outputs")
MODEL_DIR  <- file.path(ROOT_DIR, "04_model_outputs")
FIG_DIR    <- file.path(ROOT_DIR, "05_figures")
TAB_DIR    <- file.path(FIG_DIR, "source_tables")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

BIOCHEM_FILE <- file.path(
  MERGED_DIR,
  "Step03_analysis_ready_biochemical_fullcov.csv"
)

FUTURE_TAU_7Y_FILE <- file.path(
  MERGED_DIR,
  "Step03B_analysis_ready_future_tauPET_7y_fullcov.csv"
)

FUTURE_TAU_10Y_FILE <- file.path(
  MERGED_DIR,
  "Step03B_analysis_ready_future_tauPET_10y_fullcov.csv"
)

FUTURE_TAU_ANY_FILE <- file.path(
  MERGED_DIR,
  "Step03B_analysis_ready_future_tauPET_any_fullcov.csv"
)

PRIMARY_PRED_FILE <- file.path(
  MODEL_DIR,
  "Step04A_primary_interaction_prediction_grid.csv"
)

FOUR_GROUP_EMM_FILE <- file.path(
  MODEL_DIR,
  "Step04B_four_group_estimated_marginal_means.csv"
)

FOUR_GROUP_PAIR_FILE <- file.path(
  MODEL_DIR,
  "Step04B_four_group_pairwise_contrasts_tukey.csv"
)

FOUR_GROUP_DESC_FILE <- file.path(
  MODEL_DIR,
  "Step04B_four_group_biochemical_descriptives.csv"
)

KEY_TERMS_FILE <- file.path(
  MODEL_DIR,
  "Step04_manuscript_key_model_terms.csv"
)

FUTURE_PRED_7Y_FILE <- file.path(
  MODEL_DIR,
  "Step04C_future_tauPET_7y_models_prediction_grid.csv"
)

FUTURE_PRED_10Y_FILE <- file.path(
  MODEL_DIR,
  "Step04C_future_tauPET_10y_models_prediction_grid.csv"
)

FUTURE_PRED_ANY_FILE <- file.path(
  MODEL_DIR,
  "Step04C_future_tauPET_any_models_prediction_grid.csv"
)

############################################################
## 1. Helper functions
############################################################

safe_read_csv <- function(path) {
  if (!file.exists(path)) {
    message("File not found: ", path)
    return(tibble())
  }
  
  read_csv(path, show_col_types = FALSE) %>%
    clean_names()
}

make_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(x))
}

z_safe <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (sum(!is.na(x)) < 3) return(rep(NA_real_, length(x)))
  as.numeric(scale(x))
}

p_label <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "P < 0.001",
    p < 0.01 ~ paste0("P = ", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("P = ", formatC(p, format = "f", digits = 2))
  )
}

p_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

save_figure <- function(plot,
                        filename,
                        width = 8,
                        height = 6,
                        dpi = 320) {
  
  pdf_file <- file.path(FIG_DIR, paste0(filename, ".pdf"))
  png_file <- file.path(FIG_DIR, paste0(filename, ".png"))
  
  ggsave(
    pdf_file,
    plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  
  ggsave(
    png_file,
    plot,
    width = width,
    height = height,
    dpi = dpi
  )
  
  message("Saved: ", pdf_file)
  message("Saved: ", png_file)
}

theme_adni_asg <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.major.y = element_line(linewidth = 0.2, color = "grey88"),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(8, 10, 8, 10)
    )
}

group_levels <- c(
  "A- / Low ASG",
  "A- / High ASG",
  "A+ / Low ASG",
  "A+ / High ASG"
)

group_labels <- c(
  "Aβ− / Low ASG",
  "Aβ− / High ASG",
  "Aβ+ / Low ASG",
  "Aβ+ / High ASG"
)

outcome_labels <- c(
  z_csf_ptau = "CSF p-tau",
  z_csf_ttau = "CSF t-tau",
  z_tau_meta = "Future tau PET\nmeta-temporal",
  z_tau_entorhinal = "Future tau PET\nentorhinal",
  z_tau_inferior_temporal = "Future tau PET\ninferior temporal",
  z_tau_braak1 = "Future tau PET\nBraak I",
  z_tau_braak34 = "Future tau PET\nBraak III/IV",
  z_tau_braak56 = "Future tau PET\nBraak V/VI"
)

asg_level_levels <- c(
  "Low ASG (-1 SD)",
  "Mean ASG",
  "High ASG (+1 SD)"
)

asg_level_labels <- c(
  "Low ASG\n(-1 SD)",
  "Mean ASG",
  "High ASG\n(+1 SD)"
)

############################################################
## 2. Load data
############################################################

bio <- safe_read_csv(BIOCHEM_FILE)
future_tau_7y <- safe_read_csv(FUTURE_TAU_7Y_FILE)
future_tau_10y <- safe_read_csv(FUTURE_TAU_10Y_FILE)
future_tau_any <- safe_read_csv(FUTURE_TAU_ANY_FILE)

primary_pred <- safe_read_csv(PRIMARY_PRED_FILE)
four_group_emm <- safe_read_csv(FOUR_GROUP_EMM_FILE)
four_group_pair <- safe_read_csv(FOUR_GROUP_PAIR_FILE)
four_group_desc <- safe_read_csv(FOUR_GROUP_DESC_FILE)
key_terms <- safe_read_csv(KEY_TERMS_FILE)

future_pred_7y <- safe_read_csv(FUTURE_PRED_7Y_FILE) %>%
  mutate(dataset_window = "Within 7 years")

future_pred_10y <- safe_read_csv(FUTURE_PRED_10Y_FILE) %>%
  mutate(dataset_window = "Within 10 years")

future_pred_any <- safe_read_csv(FUTURE_PRED_ANY_FILE) %>%
  mutate(dataset_window = "Any future tau PET")

future_pred <- bind_rows(
  future_pred_7y,
  future_pred_10y,
  future_pred_any
)

if (nrow(bio) == 0) {
  stop("Primary biochemical dataset is empty or missing.")
}

############################################################
## 3. Standardize plotting variables
############################################################

bio <- bio %>%
  mutate(
    amyloid_asg_group = factor(
      amyloid_asg_group,
      levels = group_levels,
      labels = group_labels
    ),
    asg_group_median = factor(
      asg_group_median,
      levels = c("Low ASG", "High ASG")
    ),
    z_csf_ptau = make_numeric_safe(z_csf_ptau),
    z_csf_ttau = make_numeric_safe(z_csf_ttau),
    amyloid_centiloid = make_numeric_safe(amyloid_centiloid),
    asg_index = make_numeric_safe(asg_index),
    z_amyloid_centiloid = if ("z_amyloid_centiloid" %in% names(.)) {
      make_numeric_safe(z_amyloid_centiloid)
    } else {
      z_safe(amyloid_centiloid)
    },
    z_asg_index = if ("z_asg_index" %in% names(.)) {
      make_numeric_safe(z_asg_index)
    } else {
      z_safe(asg_index)
    }
  )

primary_pred <- primary_pred %>%
  mutate(
    outcome_label = recode(outcome, !!!outcome_labels),
    asg_level = factor(
      asg_level,
      levels = asg_level_levels,
      labels = asg_level_labels
    ),
    z_amyloid_centiloid = make_numeric_safe(z_amyloid_centiloid),
    predicted = make_numeric_safe(predicted),
    ci_low = make_numeric_safe(ci_low),
    ci_high = make_numeric_safe(ci_high)
  )

four_group_emm <- four_group_emm %>%
  mutate(
    amyloid_asg_group = factor(
      amyloid_asg_group,
      levels = group_levels,
      labels = group_labels
    ),
    outcome_label = recode(outcome, !!!outcome_labels),
    emmean = make_numeric_safe(emmean),
    se = make_numeric_safe(se),
    lower_cl = make_numeric_safe(lower_cl),
    upper_cl = make_numeric_safe(upper_cl)
  )

four_group_desc <- four_group_desc %>%
  mutate(
    amyloid_asg_group = factor(
      amyloid_asg_group,
      levels = group_levels,
      labels = group_labels
    )
  )

key_terms <- key_terms %>%
  mutate(
    estimate = make_numeric_safe(estimate),
    std_error = make_numeric_safe(std_error),
    p_value = make_numeric_safe(p_value),
    p_fdr_key_terms = make_numeric_safe(p_fdr_key_terms),
    n = make_numeric_safe(n),
    r_squared = make_numeric_safe(r_squared),
    adj_r_squared = make_numeric_safe(adj_r_squared),
    ci_low = estimate - 1.96 * std_error,
    ci_high = estimate + 1.96 * std_error,
    outcome_label = recode(outcome, !!!outcome_labels),
    term_clean = factor(
      term_clean,
      levels = c("Amyloid Centiloid", "ASG index", "Amyloid × ASG")
    )
  )

future_pred <- future_pred %>%
  mutate(
    outcome_label = recode(outcome, !!!outcome_labels),
    asg_level = factor(
      asg_level,
      levels = asg_level_levels,
      labels = asg_level_labels
    ),
    dataset_window = factor(
      dataset_window,
      levels = c("Within 7 years", "Within 10 years", "Any future tau PET")
    ),
    z_amyloid_centiloid = make_numeric_safe(z_amyloid_centiloid),
    predicted = make_numeric_safe(predicted),
    ci_low = make_numeric_safe(ci_low),
    ci_high = make_numeric_safe(ci_high)
  )

############################################################
## 4. Source tables for reproducibility
############################################################

write_csv(
  bio %>%
    select(
      rid,
      amyloid_centiloid,
      asg_index,
      amyloid_asg_group,
      z_csf_ptau,
      z_csf_ttau,
      age,
      sex,
      education,
      apoe4,
      dx_bl
    ),
  file.path(TAB_DIR, "Figure_source_primary_biochemical_data.csv")
)

write_csv(
  primary_pred,
  file.path(TAB_DIR, "Figure_source_primary_interaction_prediction_grid.csv")
)

write_csv(
  four_group_emm,
  file.path(TAB_DIR, "Figure_source_four_group_adjusted_means.csv")
)

write_csv(
  key_terms,
  file.path(TAB_DIR, "Figure_source_key_model_terms.csv")
)

############################################################
## Figure 1. Study design and ASG framework
############################################################

n_asg <- nrow(bio)
n_primary <- nrow(bio)
n_future_7y <- nrow(future_tau_7y)
n_future_10y <- nrow(future_tau_10y)

flow_nodes <- tibble(
  node = c(
    "ADNI CSF SomaScan 7K",
    "ASG index",
    "Amyloid PET Centiloid",
    "Primary biochemical outcome",
    "Secondary future tau PET",
    "Covariates"
  ),
  label = c(
    paste0("ADNI CSF SomaScan 7K\nASG proteins quantified\nn = ", n_asg),
    "Autophagy–Synapse–Glia\nimbalance index",
    paste0("Amyloid PET Centiloid\nmatched to ASG baseline\nn = ", n_primary),
    "CSF tau pathology\np-tau and t-tau\nprimary model",
    paste0("Future tau PET\n7y n = ", n_future_7y, "\n10y n = ", n_future_10y),
    "Age, sex, education,\nAPOE ε4, diagnosis"
  ),
  x = c(1, 2.5, 4, 5.5, 5.5, 4),
  y = c(3, 3, 3, 3.7, 2.3, 1.2),
  type = c("Input", "Index", "Exposure", "Outcome", "Outcome", "Covariates")
)

flow_edges <- tibble(
  x = c(1.55, 3.05, 4.55, 4.55, 4.25),
  y = c(3, 3, 3.1, 2.9, 1.5),
  xend = c(1.95, 3.45, 5.0, 5.0, 5.0),
  yend = c(3, 3, 3.55, 2.45, 2.9)
)

fig1 <- ggplot() +
  geom_segment(
    data = flow_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.55,
    color = "grey35"
  ) +
  geom_label(
    data = flow_nodes,
    aes(x = x, y = y, label = label, fill = type),
    label.size = 0.25,
    label.r = unit(0.15, "lines"),
    size = 3.6,
    lineheight = 0.96,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Input" = "#DDEAF6",
      "Index" = "#E5F3E1",
      "Exposure" = "#F7E8D0",
      "Outcome" = "#F5D6D6",
      "Covariates" = "#EEEEEE"
    )
  ) +
  coord_cartesian(xlim = c(0.4, 6.25), ylim = c(0.55, 4.35), clip = "off") +
  theme_void(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 11, hjust = 0),
    plot.margin = margin(12, 12, 12, 12)
  ) +
  labs(
    title = "Figure 1. Study design and ASG inflammatory–proteostatic framework",
    subtitle = "Baseline CSF ASG imbalance was tested as a modifier of amyloid-associated tau pathology."
  )

save_figure(
  fig1,
  "Figure1_study_design_ASG_framework",
  width = 10,
  height = 5.6
)

############################################################
## Figure 2. Primary Aβ × ASG interaction on CSF tau
############################################################

fig2 <- primary_pred %>%
  filter(outcome %in% c("z_csf_ptau", "z_csf_ttau")) %>%
  ggplot(
    aes(
      x = z_amyloid_centiloid,
      y = predicted,
      color = asg_level,
      fill = asg_level
    )
  ) +
  geom_ribbon(
    aes(ymin = ci_low, ymax = ci_high),
    alpha = 0.16,
    color = NA
  ) +
  geom_line(linewidth = 1.05) +
  facet_wrap(~ outcome_label, nrow = 1) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey55") +
  scale_color_manual(
    values = c(
      "Low ASG\n(-1 SD)" = "#2C7BB6",
      "Mean ASG" = "#555555",
      "High ASG\n(+1 SD)" = "#D7191C"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Low ASG\n(-1 SD)" = "#2C7BB6",
      "Mean ASG" = "#555555",
      "High ASG\n(+1 SD)" = "#D7191C"
    )
  ) +
  theme_adni_asg(base_size = 12) +
  labs(
    title = "Figure 2. ASG imbalance amplifies amyloid-associated CSF tau pathology",
    x = "Amyloid PET Centiloid, z-scored",
    y = "Predicted CSF tau biomarker, z-scored",
    color = "ASG level",
    fill = "ASG level"
  )

save_figure(
  fig2,
  "Figure2_primary_Amyloid_x_ASG_interaction_CSF_tau",
  width = 9.5,
  height = 5.2
)

############################################################
## Figure 3. Four-group CSF tau burden
############################################################

bio_long_tau <- bio %>%
  select(rid, amyloid_asg_group, z_csf_ptau, z_csf_ttau) %>%
  pivot_longer(
    cols = c(z_csf_ptau, z_csf_ttau),
    names_to = "outcome",
    values_to = "z_tau"
  ) %>%
  mutate(
    outcome_label = recode(outcome, !!!outcome_labels)
  ) %>%
  filter(!is.na(amyloid_asg_group), !is.na(z_tau))

fig3_raw <- bio_long_tau %>%
  ggplot(
    aes(
      x = amyloid_asg_group,
      y = z_tau,
      fill = amyloid_asg_group
    )
  ) +
  geom_violin(
    trim = FALSE,
    alpha = 0.35,
    color = NA,
    width = 0.9
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    alpha = 0.78,
    color = "black"
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.26,
    size = 0.8
  ) +
  facet_wrap(~ outcome_label, nrow = 1) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey55") +
  scale_fill_manual(
    values = c(
      "Aβ− / Low ASG" = "#9ECAE1",
      "Aβ− / High ASG" = "#4292C6",
      "Aβ+ / Low ASG" = "#FDAE6B",
      "Aβ+ / High ASG" = "#D94801"
    )
  ) +
  theme_adni_asg(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(
    title = "Figure 3A. CSF tau burden by amyloid / ASG group",
    x = NULL,
    y = "CSF tau biomarker, z-scored"
  )

fig3_emm <- four_group_emm %>%
  filter(outcome %in% c("z_csf_ptau", "z_csf_ttau")) %>%
  ggplot(
    aes(
      x = amyloid_asg_group,
      y = emmean,
      color = amyloid_asg_group
    )
  ) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey55") +
  geom_errorbar(
    aes(ymin = lower_cl, ymax = upper_cl),
    width = 0.12,
    linewidth = 0.7
  ) +
  geom_point(size = 3.0) +
  facet_wrap(~ outcome_label, nrow = 1) +
  scale_color_manual(
    values = c(
      "Aβ− / Low ASG" = "#9ECAE1",
      "Aβ− / High ASG" = "#4292C6",
      "Aβ+ / Low ASG" = "#FDAE6B",
      "Aβ+ / High ASG" = "#D94801"
    )
  ) +
  theme_adni_asg(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  ) +
  labs(
    title = "Figure 3B. Adjusted marginal means",
    x = NULL,
    y = "Adjusted CSF tau biomarker, z-scored"
  )

save_figure(
  fig3_raw,
  "Figure3A_four_group_raw_CSF_tau_distribution",
  width = 10,
  height = 5.5
)

save_figure(
  fig3_emm,
  "Figure3B_four_group_adjusted_CSF_tau_means",
  width = 10,
  height = 5.5
)

if (requireNamespace("patchwork", quietly = TRUE)) {
  fig3_combined <- fig3_raw / fig3_emm +
    patchwork::plot_annotation(
      title = "Figure 3. Aβ-positive individuals with high ASG imbalance show the highest CSF tau burden"
    )
  
  save_figure(
    fig3_combined,
    "Figure3_four_group_CSF_tau_combined",
    width = 10,
    height = 10
  )
}

############################################################
## Figure 4. Future tau PET validation
############################################################

future_interaction_terms <- key_terms %>%
  filter(
    analysis_block == "Secondary future tau PET",
    term_clean == "Amyloid × ASG",
    outcome %in% c(
      "z_tau_meta",
      "z_tau_entorhinal",
      "z_tau_inferior_temporal"
    )
  ) %>%
  mutate(
    dataset_label = factor(
      dataset_label,
      levels = c(
        "within 7 years fullcov",
        "within 10 years fullcov",
        "any future fullcov"
      ),
      labels = c(
        "Within 7 years",
        "Within 10 years",
        "Any future"
      )
    ),
    outcome_label = recode(outcome, !!!outcome_labels),
    row_label = paste0(outcome_label, " | ", dataset_label),
    significant_fdr = p_fdr_key_terms < 0.05
  ) %>%
  arrange(outcome_label, dataset_label) %>%
  mutate(
    row_label = factor(row_label, levels = rev(unique(row_label)))
  )

write_csv(
  future_interaction_terms,
  file.path(TAB_DIR, "Figure4_future_tauPET_interaction_forest_terms.csv")
)

fig4_forest <- future_interaction_terms %>%
  ggplot(
    aes(
      y = row_label,
      x = estimate,
      color = significant_fdr
    )
  ) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey50") +
  geom_segment(
    aes(x = ci_low, xend = ci_high, y = row_label, yend = row_label),
    linewidth = 0.75
  ) +
  geom_point(size = 2.8) +
  scale_color_manual(
    values = c(
      "TRUE" = "#D7191C",
      "FALSE" = "#555555"
    ),
    labels = c(
      "TRUE" = "FDR < 0.05",
      "FALSE" = "FDR ≥ 0.05"
    )
  ) +
  theme_adni_asg(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2, color = "grey88")
  ) +
  labs(
    title = "Figure 4A. Secondary prospective tau PET validation",
    subtitle = "Forest plot of Aβ × ASG interaction terms for future tau PET outcomes",
    x = "Standardized β for Aβ × ASG interaction",
    y = NULL,
    color = NULL
  )

save_figure(
  fig4_forest,
  "Figure4A_future_tauPET_interaction_forest",
  width = 9.5,
  height = 5.5
)

future_pred_it <- future_pred %>%
  filter(
    outcome == "z_tau_inferior_temporal",
    dataset_window %in% c("Within 7 years", "Within 10 years")
  )

if (nrow(future_pred_it) > 0) {
  
  fig4_pred_it <- future_pred_it %>%
    ggplot(
      aes(
        x = z_amyloid_centiloid,
        y = predicted,
        color = asg_level,
        fill = asg_level
      )
    ) +
    geom_ribbon(
      aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.15,
      color = NA
    ) +
    geom_line(linewidth = 1.0) +
    facet_wrap(~ dataset_window, nrow = 1) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey55") +
    scale_color_manual(
      values = c(
        "Low ASG\n(-1 SD)" = "#2C7BB6",
        "Mean ASG" = "#555555",
        "High ASG\n(+1 SD)" = "#D7191C"
      )
    ) +
    scale_fill_manual(
      values = c(
        "Low ASG\n(-1 SD)" = "#2C7BB6",
        "Mean ASG" = "#555555",
        "High ASG\n(+1 SD)" = "#D7191C"
      )
    ) +
    theme_adni_asg(base_size = 12) +
    labs(
      title = "Figure 4B. Aβ × ASG interaction on future inferior temporal tau PET",
      x = "Baseline amyloid PET Centiloid, z-scored",
      y = "Predicted future inferior temporal tau PET, z-scored",
      color = "ASG level",
      fill = "ASG level"
    )
  
  save_figure(
    fig4_pred_it,
    "Figure4B_future_inferior_temporal_tauPET_interaction",
    width = 9.5,
    height = 5.3
  )
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    fig4_combined <- fig4_forest / fig4_pred_it +
      patchwork::plot_annotation(
        title = "Figure 4. Secondary prospective tau PET validation"
      )
    
    save_figure(
      fig4_combined,
      "Figure4_future_tauPET_validation_combined",
      width = 10,
      height = 10
    )
  }
}

############################################################
## Figure 5. Forest plot summary of key model terms
############################################################

primary_key_forest <- key_terms %>%
  filter(
    analysis_block == "Primary biochemical",
    outcome %in% c("z_csf_ptau", "z_csf_ttau")
  ) %>%
  mutate(
    outcome_label = recode(outcome, !!!outcome_labels),
    row_label = paste0(outcome_label, " | ", term_clean),
    significant_fdr = p_fdr_key_terms < 0.05
  ) %>%
  arrange(outcome_label, term_clean) %>%
  mutate(
    row_label = factor(row_label, levels = rev(unique(row_label)))
  )

fig5_primary <- primary_key_forest %>%
  ggplot(
    aes(
      y = row_label,
      x = estimate,
      color = significant_fdr
    )
  ) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey50") +
  geom_segment(
    aes(x = ci_low, xend = ci_high, y = row_label, yend = row_label),
    linewidth = 0.75
  ) +
  geom_point(size = 3.0) +
  geom_text(
    aes(label = p_stars(p_fdr_key_terms)),
    nudge_x = 0.08,
    size = 4.2,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "TRUE" = "#D7191C",
      "FALSE" = "#555555"
    ),
    labels = c(
      "TRUE" = "FDR < 0.05",
      "FALSE" = "FDR ≥ 0.05"
    )
  ) +
  theme_adni_asg(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2, color = "grey88")
  ) +
  labs(
    title = "Figure 5A. Primary biochemical models",
    subtitle = "Standardized coefficients for key terms",
    x = "Standardized β",
    y = NULL,
    color = NULL
  )

save_figure(
  fig5_primary,
  "Figure5A_primary_biochemical_key_terms_forest",
  width = 8.8,
  height = 4.8
)

future_key_forest <- key_terms %>%
  filter(
    analysis_block == "Secondary future tau PET",
    term_clean == "Amyloid × ASG",
    outcome %in% c("z_tau_meta", "z_tau_entorhinal", "z_tau_inferior_temporal")
  ) %>%
  mutate(
    dataset_label = factor(
      dataset_label,
      levels = c(
        "within 7 years fullcov",
        "within 10 years fullcov",
        "any future fullcov"
      ),
      labels = c(
        "7y",
        "10y",
        "Any"
      )
    ),
    outcome_label = recode(outcome, !!!outcome_labels),
    row_label = paste0(outcome_label, " | ", dataset_label),
    significant_fdr = p_fdr_key_terms < 0.05
  ) %>%
  arrange(outcome_label, dataset_label) %>%
  mutate(
    row_label = factor(row_label, levels = rev(unique(row_label)))
  )

fig5_future <- future_key_forest %>%
  ggplot(
    aes(
      y = row_label,
      x = estimate,
      color = significant_fdr
    )
  ) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey50") +
  geom_segment(
    aes(x = ci_low, xend = ci_high, y = row_label, yend = row_label),
    linewidth = 0.75
  ) +
  geom_point(size = 2.8) +
  geom_text(
    aes(label = p_stars(p_fdr_key_terms)),
    nudge_x = 0.08,
    size = 4.0,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "TRUE" = "#D7191C",
      "FALSE" = "#555555"
    ),
    labels = c(
      "TRUE" = "FDR < 0.05",
      "FALSE" = "FDR ≥ 0.05"
    )
  ) +
  theme_adni_asg(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2, color = "grey88")
  ) +
  labs(
    title = "Figure 5B. Secondary future tau PET models",
    subtitle = "Aβ × ASG interaction terms",
    x = "Standardized β",
    y = NULL,
    color = NULL
  )

save_figure(
  fig5_future,
  "Figure5B_future_tauPET_interaction_terms_forest",
  width = 8.8,
  height = 5.2
)

if (requireNamespace("patchwork", quietly = TRUE)) {
  fig5_combined <- fig5_primary / fig5_future +
    patchwork::plot_annotation(
      title = "Figure 5. Summary of key model terms"
    )
  
  save_figure(
    fig5_combined,
    "Figure5_key_model_terms_forest_combined",
    width = 9.5,
    height = 10
  )
}

############################################################
## Supplementary Figure S1. Arm-specific ASG interaction forest
############################################################

ARM_FILE <- file.path(
  MODEL_DIR,
  "Step04A_sensitivity_ASG_arm_interaction_models_coefficients.csv"
)

arm_coef <- safe_read_csv(ARM_FILE)

if (nrow(arm_coef) > 0) {
  
  arm_interaction <- arm_coef %>%
    filter(str_detect(term, ":")) %>%
    mutate(
      estimate = make_numeric_safe(estimate),
      std_error = make_numeric_safe(std_error),
      p_value = make_numeric_safe(p_value),
      ci_low = estimate - 1.96 * std_error,
      ci_high = estimate + 1.96 * std_error,
      outcome_label = recode(outcome, !!!outcome_labels),
      arm = case_when(
        str_detect(term, "autophagy") ~ "Autophagy arm",
        str_detect(term, "synaptic") ~ "Synaptic arm",
        str_detect(term, "glial") ~ "Glial arm",
        TRUE ~ term
      ),
      p_fdr_arm = p.adjust(p_value, method = "BH"),
      significant_fdr = p_fdr_arm < 0.05,
      row_label = paste0(outcome_label, " | ", arm)
    ) %>%
    arrange(outcome_label, arm) %>%
    mutate(
      row_label = factor(row_label, levels = rev(unique(row_label)))
    )
  
  write_csv(
    arm_interaction,
    file.path(TAB_DIR, "FigureS1_ASG_arm_interaction_terms.csv")
  )
  
  fig_s1_arm <- arm_interaction %>%
    ggplot(
      aes(
        y = row_label,
        x = estimate,
        color = significant_fdr
      )
    ) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey50") +
    geom_segment(
      aes(x = ci_low, xend = ci_high, y = row_label, yend = row_label),
      linewidth = 0.75
    ) +
    geom_point(size = 2.8) +
    scale_color_manual(
      values = c(
        "TRUE" = "#D7191C",
        "FALSE" = "#555555"
      ),
      labels = c(
        "TRUE" = "FDR < 0.05",
        "FALSE" = "FDR ≥ 0.05"
      )
    ) +
    theme_adni_asg(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.2, color = "grey88")
    ) +
    labs(
      title = "Figure S1. Arm-specific Aβ × ASG component interactions",
      x = "Standardized β for interaction term",
      y = NULL,
      color = NULL
    )
  
  save_figure(
    fig_s1_arm,
    "FigureS1_ASG_arm_specific_interaction_forest",
    width = 9,
    height = 5.2
  )
}

############################################################
## Supplementary Figure S2. A+ subgroup key model terms
############################################################

APLUS_FILE <- file.path(
  MODEL_DIR,
  "Step04D_Aplus_subgroup_biochemical_models.csv"
)

aplus_coef <- safe_read_csv(APLUS_FILE)

if (nrow(aplus_coef) > 0) {
  
  aplus_key <- aplus_coef %>%
    filter(term %in% c("z_asg_index", "z_amyloid_centiloid")) %>%
    mutate(
      estimate = make_numeric_safe(estimate),
      std_error = make_numeric_safe(std_error),
      p_value = make_numeric_safe(p_value),
      ci_low = estimate - 1.96 * std_error,
      ci_high = estimate + 1.96 * std_error,
      outcome_label = recode(outcome, !!!outcome_labels),
      term_clean = case_when(
        term == "z_asg_index" ~ "ASG index",
        term == "z_amyloid_centiloid" ~ "Amyloid Centiloid",
        TRUE ~ term
      ),
      p_fdr_aplus = p.adjust(p_value, method = "BH"),
      significant_fdr = p_fdr_aplus < 0.05,
      row_label = paste0(outcome_label, " | ", term_clean)
    ) %>%
    arrange(outcome_label, term_clean) %>%
    mutate(
      row_label = factor(row_label, levels = rev(unique(row_label)))
    )
  
  write_csv(
    aplus_key,
    file.path(TAB_DIR, "FigureS2_Aplus_subgroup_key_terms.csv")
  )
  
  fig_s2_aplus <- aplus_key %>%
    ggplot(
      aes(
        y = row_label,
        x = estimate,
        color = significant_fdr
      )
    ) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey50") +
    geom_segment(
      aes(x = ci_low, xend = ci_high, y = row_label, yend = row_label),
      linewidth = 0.75
    ) +
    geom_point(size = 2.8) +
    scale_color_manual(
      values = c(
        "TRUE" = "#D7191C",
        "FALSE" = "#555555"
      ),
      labels = c(
        "TRUE" = "FDR < 0.05",
        "FALSE" = "FDR ≥ 0.05"
      )
    ) +
    theme_adni_asg(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.2, color = "grey88")
    ) +
    labs(
      title = "Figure S2. Aβ-positive subgroup biochemical models",
      x = "Standardized β",
      y = NULL,
      color = NULL
    )
  
  save_figure(
    fig_s2_aplus,
    "FigureS2_Aplus_subgroup_biochemical_forest",
    width = 8.5,
    height = 4.6
  )
}

############################################################
## 5. Figure index
############################################################

figure_index <- tibble(
  figure = c(
    "Figure 1",
    "Figure 2",
    "Figure 3A",
    "Figure 3B",
    "Figure 3 combined",
    "Figure 4A",
    "Figure 4B",
    "Figure 4 combined",
    "Figure 5A",
    "Figure 5B",
    "Figure 5 combined",
    "Figure S1",
    "Figure S2"
  ),
  file_prefix = c(
    "Figure1_study_design_ASG_framework",
    "Figure2_primary_Amyloid_x_ASG_interaction_CSF_tau",
    "Figure3A_four_group_raw_CSF_tau_distribution",
    "Figure3B_four_group_adjusted_CSF_tau_means",
    "Figure3_four_group_CSF_tau_combined",
    "Figure4A_future_tauPET_interaction_forest",
    "Figure4B_future_inferior_temporal_tauPET_interaction",
    "Figure4_future_tauPET_validation_combined",
    "Figure5A_primary_biochemical_key_terms_forest",
    "Figure5B_future_tauPET_interaction_terms_forest",
    "Figure5_key_model_terms_forest_combined",
    "FigureS1_ASG_arm_specific_interaction_forest",
    "FigureS2_Aplus_subgroup_biochemical_forest"
  ),
  purpose = c(
    "Study design and biological framework",
    "Primary Aβ × ASG interaction on CSF tau",
    "Raw four-group CSF tau distribution",
    "Adjusted marginal means for four-group CSF tau",
    "Combined four-group figure",
    "Future tau PET interaction forest",
    "Future inferior temporal tau PET interaction plot",
    "Combined future tau PET validation figure",
    "Primary biochemical forest plot",
    "Future tau PET interaction forest plot",
    "Combined key model terms forest plot",
    "Arm-specific ASG sensitivity",
    "Aβ-positive subgroup sensitivity"
  )
)

write_csv(
  figure_index,
  file.path(FIG_DIR, "Step05_figure_index.csv")
)

############################################################
## 6. Final messages
############################################################

message("\n============================================================")
message("Step 05 finished.")
message("Figures saved to:")
message(FIG_DIR)
message("Source plotting tables saved to:")
message(TAB_DIR)
message("============================================================")

message("\nMain figures to inspect first:")
message("1) Figure2_primary_Amyloid_x_ASG_interaction_CSF_tau.pdf")
message("2) Figure3_four_group_CSF_tau_combined.pdf")
message("3) Figure4_future_tauPET_validation_combined.pdf")
message("4) Figure5_key_model_terms_forest_combined.pdf")

message("\nFigure index:")
message(file.path(FIG_DIR, "Step05_figure_index.csv"))
