options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
output_root <- Sys.getenv("BMI_BP_OUTPUT_DIR", file.path(root, "outputs"))
out <- file.path(output_root, "hrs_analysis")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(out, "analysis_intervals.csv"), fileEncoding = "UTF-8", check.names = FALSE)

d$interval <- factor(d$interval, levels = c("2014-2018", "2018-2022"))
d$race_f <- factor(d$raracem)
d$region4 <- ifelse(d$region14 %in% c(1, 2), 1,
                    ifelse(d$region14 %in% c(3, 4), 2,
                           ifelse(d$region14 %in% c(5, 6, 7), 3,
                                  ifelse(d$region14 %in% c(8, 9), 4, NA))))
d$region4_f <- factor(d$region4)
d$rural_f <- factor(d$rural, levels = c(0, 1), labels = c("Non-rural", "Rural"))
d$female_f <- factor(d$female, levels = c(0, 1), labels = c("Men", "Women"))
d$age65_f <- factor(d$age65, levels = c(0, 1), labels = c("Age <65", "Age >=65"))
d$obese_f <- factor(d$obese_2014, levels = c(0, 1), labels = c("BMI <30", "BMI >=30"))
d$change_group <- factor(d$bmi_change_cat_code, levels = c(0, -1, 1), labels = c("Stable", "Loss >5%", "Gain >5%"))
d$burden_q <- factor(d$burden_q, levels = c("Q1", "Q2", "Q3", "Q4"))
d$vim_q <- factor(d$vim_q, levels = c("Q1", "Q2", "Q3", "Q4"))

covars <- c(
  "age_start10", "female", "race_f", "raedyrs", "married_partnered",
  "current_smoker", "current_drinker", "r12diabe", "cvd_history",
  "r12shlt", "region4_f", "rural", "interval"
)

endpoint_labels <- c(
  hypertension_onset = "Hypertension onset",
  treatment_initiation = "Treatment initiation",
  untreated_bp_improvement = "Untreated BP improvement",
  control_loss = "Loss of BP control",
  control_achievement = "Achievement of BP control",
  death_next_wave = "Death"
)

primary_endpoints <- names(endpoint_labels)[1:5]
all_endpoints <- names(endpoint_labels)

model_data <- function(endpoint, weight_type = "ipcw", extra_filter = rep(TRUE, nrow(d))) {
  weight_var <- if (endpoint == "death_next_wave" || weight_type == "base") "base_weight" else "model_weight"
  needed <- c(endpoint, weight_var, "secu", "stratum", covars,
              "excess_bmi25_z", "bmi_vim_z", "bmi_mean_z", "bmi_auc_z",
              "bmi_arv_z", "bmi_cv_z", "bmi_sd_z", "change_group",
              "burden_q", "vim_q", "female_f", "age65_f", "obese_f", "rural_f")
  keep <- !is.na(d[[endpoint]]) & !is.na(d[[weight_var]]) & d[[weight_var]] > 0 & extra_filter
  keep <- keep & complete.cases(d[, needed, drop = FALSE])
  z <- d[keep, , drop = FALSE]
  z$analysis_w <- if (weight_type == "none") 1 else z[[weight_var]]
  z
}

fit_svy_rr <- function(endpoint, exposure_terms, model_name,
                       weight_type = "ipcw", extra_filter = rep(TRUE, nrow(d)),
                       keep_terms = exposure_terms) {
  z <- model_data(endpoint, weight_type, extra_filter)
  if (nrow(z) < 80 || sum(z[[endpoint]]) < 20 || length(unique(z[[endpoint]])) < 2) return(NULL)
  des <- svydesign(ids = ~secu, strata = ~stratum, weights = ~analysis_w, data = z, nest = TRUE)
  f <- as.formula(paste(endpoint, "~", paste(c(exposure_terms, covars), collapse = " + ")))
  fit <- tryCatch(svyglm(f, design = des, family = quasipoisson(link = "log")), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$coefficients
  out_rows <- list()
  for (term in rownames(co)) {
    if (!any(vapply(keep_terms, function(k) startsWith(term, k), logical(1)))) next
    estimate <- unname(co[term, 1])
    se <- unname(co[term, 2])
    p <- unname(co[term, ncol(co)])
    out_rows[[length(out_rows) + 1]] <- data.frame(
      endpoint = endpoint,
      endpoint_label = endpoint_labels[[endpoint]] %||% endpoint,
      model = model_name,
      term = term,
      log_rr = estimate,
      se = se,
      rr = exp(estimate),
      ci_low = exp(estimate - 1.96 * se),
      ci_high = exp(estimate + 1.96 * se),
      p_value = p,
      n_intervals = nrow(z),
      events = sum(z[[endpoint]]),
      participants = length(unique(z$person_id)),
      psu = length(unique(z$secu)),
      strata = length(unique(z$stratum)),
      weight_type = weight_type,
      stringsAsFactors = FALSE
    )
  }
  if (length(out_rows) == 0) return(NULL)
  do.call(rbind, out_rows)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# Main models: cumulative excess BMI burden and VIM adjusted for mean BMI.
main_list <- list()
for (ep in all_endpoints) {
  main_list[[length(main_list) + 1]] <- fit_svy_rr(ep, "excess_bmi25_z", "Cumulative excess BMI burden", keep_terms = "excess_bmi25_z")
  main_list[[length(main_list) + 1]] <- fit_svy_rr(ep, c("bmi_vim_z", "bmi_mean_z"), "BMI variability (VIM)", keep_terms = "bmi_vim_z")
  main_list[[length(main_list) + 1]] <- fit_svy_rr(ep, "bmi_mean_z", "Mean BMI", keep_terms = "bmi_mean_z")
}
main_results <- do.call(rbind, Filter(Negate(is.null), main_list))
main_results$p_fdr <- NA_real_
primary_ix <- main_results$endpoint %in% primary_endpoints & main_results$model %in% c("Cumulative excess BMI burden", "BMI variability (VIM)")
main_results$p_fdr[primary_ix] <- p.adjust(main_results$p_value[primary_ix], method = "BH")
write.csv(main_results, file.path(out, "main_transition_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Alternative exposure definitions and analytic assumptions.
sens_list <- list()
for (ep in primary_endpoints) {
  sens_list[[length(sens_list) + 1]] <- fit_svy_rr(ep, "bmi_auc_z", "Alternative total BMI-years", keep_terms = "bmi_auc_z")
  sens_list[[length(sens_list) + 1]] <- fit_svy_rr(ep, c("bmi_arv_z", "bmi_mean_z"), "Alternative variability ARV", keep_terms = "bmi_arv_z")
  sens_list[[length(sens_list) + 1]] <- fit_svy_rr(ep, c("bmi_cv_z", "bmi_mean_z"), "Alternative variability CV", keep_terms = "bmi_cv_z")
  sens_list[[length(sens_list) + 1]] <- fit_svy_rr(ep, c("bmi_sd_z", "bmi_mean_z"), "Alternative variability SD", keep_terms = "bmi_sd_z")
  sens_list[[length(sens_list) + 1]] <- fit_svy_rr(ep, "excess_bmi25_z", "Base survey weight only", weight_type = "base", keep_terms = "excess_bmi25_z")
  sens_list[[length(sens_list) + 1]] <- fit_svy_rr(ep, "excess_bmi25_z", "Unweighted", weight_type = "none", keep_terms = "excess_bmi25_z")
  sens_list[[length(sens_list) + 1]] <- fit_svy_rr(ep, "excess_bmi25_z", "Exclude baseline CVD", extra_filter = d$cvd_history == 0, keep_terms = "excess_bmi25_z")
  sens_list[[length(sens_list) + 1]] <- fit_svy_rr(ep, "excess_bmi25_z", "Exclude >5% BMI loss", extra_filter = d$bmi_change_cat_code != -1, keep_terms = "excess_bmi25_z")
  sens_list[[length(sens_list) + 1]] <- fit_svy_rr(ep, c("bmi_vim_z", "bmi_mean_z"), "VIM excluding >5% BMI loss", extra_filter = d$bmi_change_cat_code != -1, keep_terms = "bmi_vim_z")
}

# ACC/AHA threshold sensitivity maps to the three comparable endpoints.
accaha_map <- c(
  hypertension_onset_accaha = "Hypertension onset",
  control_loss_accaha = "Loss of BP control",
  control_achievement_accaha = "Achievement of BP control"
)
for (ep in names(accaha_map)) {
  endpoint_labels[[ep]] <- accaha_map[[ep]]
  sens_list[[length(sens_list) + 1]] <- fit_svy_rr(ep, "excess_bmi25_z", "ACC/AHA threshold: burden", keep_terms = "excess_bmi25_z")
  sens_list[[length(sens_list) + 1]] <- fit_svy_rr(ep, c("bmi_vim_z", "bmi_mean_z"), "ACC/AHA threshold: VIM", keep_terms = "bmi_vim_z")
}
sensitivity_results <- do.call(rbind, Filter(Negate(is.null), sens_list))
write.csv(sensitivity_results, file.path(out, "sensitivity_transition_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Weight-change categories relative to stable BMI.
change_list <- list()
for (ep in primary_endpoints) {
  change_list[[length(change_list) + 1]] <- fit_svy_rr(
    ep, "change_group", "BMI change category", keep_terms = "change_group"
  )
}
change_results <- do.call(rbind, Filter(Negate(is.null), change_list))
change_results$p_fdr <- p.adjust(change_results$p_value, method = "BH")
write.csv(change_results, file.path(out, "weight_change_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Quartile models for non-linearity screening.
quartile_list <- list()
for (ep in primary_endpoints) {
  quartile_list[[length(quartile_list) + 1]] <- fit_svy_rr(ep, "burden_q", "Burden quartiles", keep_terms = "burden_q")
  quartile_list[[length(quartile_list) + 1]] <- fit_svy_rr(ep, c("vim_q", "bmi_mean_z"), "VIM quartiles", keep_terms = "vim_q")
}
quartile_results <- do.call(rbind, Filter(Negate(is.null), quartile_list))
quartile_results$p_fdr <- p.adjust(quartile_results$p_value, method = "BH")
write.csv(quartile_results, file.path(out, "quartile_transition_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Interaction tests and stratified estimates for prespecified subgroups.
interaction_test <- function(endpoint, exposure, subgroup, model_label) {
  z <- model_data(endpoint)
  if (nrow(z) < 150 || sum(z[[endpoint]]) < 30 || length(unique(z[[subgroup]])) < 2) return(NULL)
  des <- svydesign(ids = ~secu, strata = ~stratum, weights = ~analysis_w, data = z, nest = TRUE)
  extra <- if (exposure == "bmi_vim_z") "bmi_mean_z" else character(0)
  remove_covars <- character(0)
  if (subgroup == "female_f") remove_covars <- c(remove_covars, "female")
  if (subgroup == "rural_f") remove_covars <- c(remove_covars, "rural")
  main_covars <- setdiff(covars, remove_covars)
  f <- as.formula(paste(endpoint, "~", paste(c(paste0(exposure, "*", subgroup), extra, main_covars), collapse = " + ")))
  fit <- tryCatch(svyglm(f, design = des, family = quasipoisson(link = "log")), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  interaction_terms <- grep(paste0("^", exposure, ":", subgroup), names(coef(fit)), value = TRUE)
  if (length(interaction_terms) == 0) interaction_terms <- grep(paste0("^", subgroup, ".*:", exposure), names(coef(fit)), value = TRUE)
  if (length(interaction_terms) == 0) return(NULL)
  b <- coef(fit)[interaction_terms[1]]
  se <- sqrt(vcov(fit)[interaction_terms[1], interaction_terms[1]])
  p <- 2 * pnorm(abs(b / se), lower.tail = FALSE)
  data.frame(endpoint = endpoint, endpoint_label = endpoint_labels[[endpoint]], exposure = exposure,
             subgroup = subgroup, model = model_label, interaction_log_rr = b,
             interaction_rr = exp(b), ci_low = exp(b - 1.96 * se), ci_high = exp(b + 1.96 * se),
             p_interaction = p, n_intervals = nrow(z), events = sum(z[[endpoint]]))
}

interaction_list <- list()
subgroups <- c("female_f", "age65_f", "obese_f", "rural_f")
for (ep in c("hypertension_onset", "control_loss", "control_achievement")) {
  for (sg in subgroups) {
    interaction_list[[length(interaction_list) + 1]] <- interaction_test(ep, "excess_bmi25_z", sg, "Burden interaction")
    interaction_list[[length(interaction_list) + 1]] <- interaction_test(ep, "bmi_vim_z", sg, "VIM interaction")
  }
}
interaction_results <- do.call(rbind, Filter(Negate(is.null), interaction_list))
interaction_results$p_fdr <- p.adjust(interaction_results$p_interaction, method = "BH")
write.csv(interaction_results, file.path(out, "subgroup_interactions.csv"), row.names = FALSE, fileEncoding = "UTF-8")

fit_stratum_rr <- function(endpoint, exposure_terms, keep_term, subgroup_var, subgroup_value, subgroup_label) {
  z <- model_data(endpoint)
  z <- z[z[[subgroup_var]] == subgroup_value, , drop = FALSE]
  if (nrow(z) < 80 || sum(z[[endpoint]]) < 20 || length(unique(z[[endpoint]])) < 2) return(NULL)
  local_covars <- covars
  if (subgroup_var == "female") local_covars <- setdiff(local_covars, "female")
  if (subgroup_var == "rural") local_covars <- setdiff(local_covars, "rural")
  des <- svydesign(ids = ~secu, strata = ~stratum, weights = ~analysis_w, data = z, nest = TRUE)
  f <- as.formula(paste(endpoint, "~", paste(c(exposure_terms, local_covars), collapse = " + ")))
  fit <- tryCatch(svyglm(f, design = des, family = quasipoisson(link = "log")), error = function(e) NULL)
  if (is.null(fit) || !(keep_term %in% names(coef(fit)))) return(NULL)
  b <- coef(fit)[keep_term]
  se <- sqrt(vcov(fit)[keep_term, keep_term])
  p <- 2 * pnorm(abs(b / se), lower.tail = FALSE)
  data.frame(endpoint = endpoint, endpoint_label = endpoint_labels[[endpoint]],
             exposure = ifelse(keep_term == "excess_bmi25_z", "Cumulative excess BMI burden", "BMI variability (VIM)"),
             subgroup = subgroup_var, stratum = subgroup_label,
             rr = exp(b), ci_low = exp(b - 1.96 * se), ci_high = exp(b + 1.96 * se), p_value = p,
             n_intervals = nrow(z), events = sum(z[[endpoint]]), participants = length(unique(z$person_id)))
}

strata_spec <- list(
  list(var = "female", values = c(0, 1), labels = c("Men", "Women")),
  list(var = "age65", values = c(0, 1), labels = c("Age <65", "Age >=65")),
  list(var = "obese_2014", values = c(0, 1), labels = c("BMI <30", "BMI >=30")),
  list(var = "rural", values = c(0, 1), labels = c("Non-rural", "Rural"))
)
stratified_list <- list()
for (ep in c("hypertension_onset", "control_loss", "control_achievement")) {
  for (spec in strata_spec) {
    for (j in seq_along(spec$values)) {
      stratified_list[[length(stratified_list) + 1]] <- fit_stratum_rr(ep, "excess_bmi25_z", "excess_bmi25_z", spec$var, spec$values[j], spec$labels[j])
      stratified_list[[length(stratified_list) + 1]] <- fit_stratum_rr(ep, c("bmi_vim_z", "bmi_mean_z"), "bmi_vim_z", spec$var, spec$values[j], spec$labels[j])
    }
  }
}
stratified_results <- do.call(rbind, Filter(Negate(is.null), stratified_list))
write.csv(stratified_results, file.path(out, "subgroup_stratified_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Model-level decision summary.
decision <- data.frame(
  item = c("Effect scale", "Primary burden", "Primary variability", "Multiplicity", "Survey design", "Time interpretation"),
  decision = c(
    "Four-year transition risk ratio from survey-weighted modified Poisson models",
    "Cumulative excess BMI-years above BMI 25, per cohort SD",
    "VIM per cohort SD, adjusted for mean BMI",
    "Benjamini-Hochberg FDR across 10 primary burden/VIM transition tests",
    "HRS physical-measurement weight x exposure-completeness weight x IPCW; PSU and strata retained",
    "Panel-state transitions; exact event times are not observed"
  )
)
write.csv(decision, file.path(out, "statistical_decisions.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Figure 1: transition probability heatmap (descriptive evidence).
tr <- read.csv(file.path(out, "transition_counts.csv"), fileEncoding = "UTF-8-BOM")
tr$total_origin <- ave(tr$n, tr$interval, tr$origin, FUN = sum)
tr$probability <- ifelse(tr$total_origin > 0, tr$n / tr$total_origin, NA)
tr$origin_f <- factor(tr$origin, levels = 5:1,
                      labels = c("Treated uncontrolled", "Treated controlled", "Untreated hypertension", "Elevated untreated", "Normal untreated"))
tr$destination_f <- factor(tr$destination, levels = 1:5,
                           labels = c("Normal untreated", "Elevated untreated", "Untreated hypertension", "Treated controlled", "Treated uncontrolled"))

theme_set(theme_classic(base_size = 7, base_family = "sans") +
            theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
                  plot.title = element_text(size = 8, face = "bold"), strip.text = element_text(size = 7, face = "bold"),
                  legend.title = element_text(size = 6.5), legend.text = element_text(size = 6), panel.grid = element_blank()))

p_heat <- ggplot(tr, aes(destination_f, origin_f, fill = probability)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n > 0, sprintf("%d\n%.1f%%", n, 100 * probability), "")), size = 2.0, lineheight = 0.9) +
  facet_wrap(~interval, nrow = 1) +
  scale_fill_gradient(low = "#F3F6F8", high = "#3B6C8E", labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Destination state", y = "Origin state", fill = "Row probability",
       title = "Observed four-year transitions among participants with blood pressure states at both waves") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

# Figure 2: main effect forest plot.
forest <- subset(main_results, model %in% c("Cumulative excess BMI burden", "BMI variability (VIM)") & endpoint %in% primary_endpoints)
forest$model <- factor(forest$model, levels = c("Cumulative excess BMI burden", "BMI variability (VIM)"),
                       labels = c("Excess BMI burden", "BMI variability (VIM)"))
forest$endpoint_label <- factor(forest$endpoint_label, levels = rev(unname(endpoint_labels[primary_endpoints])))

p_forest <- ggplot(forest, aes(rr, endpoint_label, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = "#777777", linewidth = 0.35, linetype = 2) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.18,
                position = position_dodge(width = 0.45), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.45), size = 2.0) +
  scale_x_log10() +
  scale_colour_manual(values = c("#34688A", "#B0684B")) +
  labs(x = "Four-year transition risk ratio per 1 SD (log scale)", y = NULL, colour = NULL, shape = NULL,
       title = "HRS transition estimates for BMI burden and variability") +
  theme(legend.position = "top")

save_pub <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(file.path(out, paste0(stem, ".svg")), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(file.path(out, paste0(stem, ".pdf")), width = w, height = h, family = "sans")
  print(plot)
  dev.off()
  ragg::agg_tiff(file.path(out, paste0(stem, ".tiff")), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(file.path(out, paste0(stem, "_preview.png")), width = w, height = h, units = "in", res = 180)
  print(plot)
  dev.off()
}

save_pub(p_heat, "figure1_transition_heatmap", 183, 105)
save_pub(p_forest, "figure2_main_forest", 150, 100)

figure_qa <- data.frame(
  figure = c("Figure 1", "Figure 2"),
  core_conclusion = c("Blood pressure and treatment states are dynamic but state persistence is common.",
                      "Associations of BMI burden and variability differ across transition types."),
  archetype = c("Quantitative grid", "Quantitative grid"),
  source_data = c("transition_counts.csv", "main_transition_models.csv"),
  statistics = c("Observed counts and row percentages", "Survey-weighted modified Poisson RR with 95% CI"),
  exports = c("SVG/PDF/TIFF/PNG", "SVG/PDF/TIFF/PNG")
)
write.csv(figure_qa, file.path(out, "figure_qa_contract.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("Main models:\n")
print(main_results[, c("endpoint", "model", "rr", "ci_low", "ci_high", "p_value", "p_fdr", "n_intervals", "events")])
cat("\nOutputs:", out, "\n")
