options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
output_root <- Sys.getenv("BMI_BP_OUTPUT_DIR", file.path(root, "outputs"))
out <- file.path(output_root, "elsa_analysis")
d <- read.csv(file.path(out, "elsa_intervals.csv"), fileEncoding = "UTF-8", check.names = FALSE)

d$interval <- factor(d$interval, levels = c("W6-W8", "W6-W9"))
d$education_f <- factor(d$education)
d$region_f <- factor(d$region)
d$physical_activity_f <- factor(d$physical_activity)
d$female_f <- factor(d$female, levels = c(0, 1), labels = c("Men", "Women"))
d$age65_f <- factor(d$age65, levels = c(0, 1), labels = c("Age <65", "Age >=65"))
d$obese_f <- factor(d$obese_w6, levels = c(0, 1), labels = c("BMI <30", "BMI >=30"))
d$rural_f <- factor(d$rural, levels = c(0, 1), labels = c("Urban", "Rural"))
d$change_group <- factor(d$bmi_change_cat_code, levels = c(0, -1, 1), labels = c("Stable", "Loss >5%", "Gain >5%"))
d$burden_q <- factor(d$burden_q, levels = c("Q1", "Q2", "Q3", "Q4"))
d$vim_q <- factor(d$vim_q, levels = c("Q1", "Q2", "Q3", "Q4"))

endpoint_labels <- c(
  hypertension_onset = "Hypertension onset",
  treatment_initiation = "Treatment initiation",
  untreated_bp_improvement = "Untreated BP improvement",
  control_loss = "Loss of BP control",
  control_achievement = "Achievement of BP control"
)
endpoints <- names(endpoint_labels)

full_covars <- c(
  "age_start10", "female", "education_f", "married_partnered", "current_smoker", "current_drinker",
  "diabetes", "cvd_history", "self_health", "region_f", "rural", "interval"
)

covars_for <- function(endpoint, add_physical_activity = FALSE) {
  base_endpoint <- sub("_accaha$", "", endpoint)
  result <- if (base_endpoint == "control_loss") {
    c("age_start10", "female", "region_f", "rural", "interval")
  } else if (base_endpoint == "control_achievement") {
    c("age_start10", "female", "education_f", "diabetes", "region_f", "rural", "interval")
  } else {
    full_covars
  }
  if (add_physical_activity) result <- unique(c(result, "physical_activity_f"))
  result
}

model_data <- function(endpoint, exposure_terms, weight_type = "ipcw", extra_filter = rep(TRUE, nrow(d)),
                       add_physical_activity = FALSE, additional = character(0)) {
  weight_var <- switch(weight_type, ipcw = "model_weight", base = "base_weight",
                       outcome_nurse = "outcome_nurse_weight", none = NA_character_)
  needed <- unique(c(endpoint, "idahhw6", covars_for(endpoint, add_physical_activity), exposure_terms, additional))
  if (!is.na(weight_var)) needed <- unique(c(needed, weight_var))
  keep <- !is.na(d[[endpoint]]) & extra_filter & complete.cases(d[, needed, drop = FALSE])
  if (!is.na(weight_var)) keep <- keep & !is.na(d[[weight_var]]) & d[[weight_var]] > 0
  z <- d[keep, , drop = FALSE]
  z$analysis_w <- if (weight_type == "none") rep(1, nrow(z)) else z[[weight_var]]
  z
}

fit_rr <- function(endpoint, exposure_terms, model_name, keep_terms = exposure_terms,
                   weight_type = "ipcw", extra_filter = rep(TRUE, nrow(d)), add_physical_activity = FALSE) {
  z <- model_data(endpoint, exposure_terms, weight_type, extra_filter, add_physical_activity)
  if (nrow(z) < 60 || sum(z[[endpoint]]) < 20 || length(unique(z[[endpoint]])) < 2) return(NULL)
  des <- svydesign(ids = ~idahhw6, weights = ~analysis_w, data = z)
  formula <- as.formula(paste(endpoint, "~", paste(c(exposure_terms, covars_for(endpoint, add_physical_activity)), collapse = " + ")))
  fit <- tryCatch(svyglm(formula, design = des, family = quasipoisson(link = "log")), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  coefficients <- summary(fit)$coefficients
  rows <- list()
  for (term in rownames(coefficients)) {
    if (!any(vapply(keep_terms, function(key) startsWith(term, key), logical(1)))) next
    beta <- unname(coefficients[term, 1]); se <- unname(coefficients[term, 2]); p <- unname(coefficients[term, ncol(coefficients)])
    rows[[length(rows) + 1]] <- data.frame(
      cohort = "ELSA", endpoint = endpoint, endpoint_label = endpoint_labels[[sub("_accaha$", "", endpoint)]],
      model = model_name, term = term, log_rr = beta, se = se, rr = exp(beta),
      ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se), p_value = p,
      n_intervals = nrow(z), events = sum(z[[endpoint]]), participants = length(unique(z$person_id)),
      households = length(unique(z$idahhw6)), weight_type = weight_type
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

fit_tir <- function(endpoint, exposure_terms, model_name, keep_term) {
  z <- model_data(endpoint, exposure_terms, "ipcw", additional = "interval_years")
  if (nrow(z) < 60 || sum(z[[endpoint]]) < 20 || length(unique(z[[endpoint]])) < 2) return(NULL)
  des <- svydesign(ids = ~idahhw6, weights = ~analysis_w, data = z)
  formula <- as.formula(paste(endpoint, "~", paste(c(exposure_terms, covars_for(endpoint), "offset(log(interval_years))"), collapse = " + ")))
  fit <- tryCatch(svyglm(formula, design = des, family = quasibinomial(link = "cloglog")), error = function(e) NULL)
  if (is.null(fit) || !(keep_term %in% names(coef(fit)))) return(NULL)
  beta <- coef(fit)[keep_term]
  se <- sqrt(vcov(fit)[keep_term, keep_term])
  p <- 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  data.frame(
    cohort = "ELSA", endpoint = endpoint, endpoint_label = endpoint_labels[[endpoint]], model = model_name,
    term = keep_term, log_tir = beta, se = se, tir = exp(beta), ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se), p_value = p, n_intervals = nrow(z), events = sum(z[[endpoint]]),
    participants = length(unique(z$person_id)), interval_years = "4-6"
  )
}

# Primary HRS/CHNS-compatible modified Poisson models.
main_list <- list()
tir_list <- list()
for (endpoint in endpoints) {
  main_list[[length(main_list) + 1]] <- fit_rr(endpoint, "excess_bmi25_z", "Cumulative excess BMI burden", "excess_bmi25_z")
  main_list[[length(main_list) + 1]] <- fit_rr(endpoint, c("bmi_vim_z", "bmi_mean_z"), "BMI variability (VIM)", "bmi_vim_z")
  main_list[[length(main_list) + 1]] <- fit_rr(endpoint, "bmi_mean_z", "Mean BMI", "bmi_mean_z")
  tir_list[[length(tir_list) + 1]] <- fit_tir(endpoint, "excess_bmi25_z", "Cumulative excess BMI burden", "excess_bmi25_z")
  tir_list[[length(tir_list) + 1]] <- fit_tir(endpoint, c("bmi_vim_z", "bmi_mean_z"), "BMI variability (VIM)", "bmi_vim_z")
}
main <- do.call(rbind, Filter(Negate(is.null), main_list))
primary <- main$model %in% c("Cumulative excess BMI burden", "BMI variability (VIM)")
main$p_fdr <- NA_real_
main$p_fdr[primary] <- p.adjust(main$p_value[primary], method = "BH")
write.csv(main, file.path(out, "elsa_main_transition_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

tir <- do.call(rbind, Filter(Negate(is.null), tir_list))
tir$p_fdr <- p.adjust(tir$p_value, method = "BH")
write.csv(tir, file.path(out, "elsa_harmonized_tir.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Sensitivity analyses matching the prior cohorts.
sensitivity_list <- list()
for (endpoint in endpoints) {
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, "bmi_auc_z", "Alternative total BMI-years", "bmi_auc_z")
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, c("bmi_arv_z", "bmi_mean_z"), "Alternative variability ARV", "bmi_arv_z")
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, c("bmi_cv_z", "bmi_mean_z"), "Alternative variability CV", "bmi_cv_z")
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, c("bmi_sd_z", "bmi_mean_z"), "Alternative variability SD", "bmi_sd_z")
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, "excess_bmi25_z", "Baseline/exposure weight only", "excess_bmi25_z", weight_type = "base")
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, "excess_bmi25_z", "Outcome nurse weight", "excess_bmi25_z", weight_type = "outcome_nurse")
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, "excess_bmi25_z", "Unweighted", "excess_bmi25_z", weight_type = "none")
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, "excess_bmi25_z", "Exclude baseline CVD", "excess_bmi25_z", extra_filter = d$cvd_history == 0)
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, "excess_bmi25_z", "Exclude >5% BMI loss", "excess_bmi25_z", extra_filter = d$bmi_change_cat_code != -1)
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, c("bmi_vim_z", "bmi_mean_z"), "VIM excluding >5% BMI loss", "bmi_vim_z", extra_filter = d$bmi_change_cat_code != -1)
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, "excess_bmi25_z", "Physical activity adjusted", "excess_bmi25_z", add_physical_activity = TRUE)
}

accaha_map <- setNames(paste0(endpoints, "_accaha"), endpoints)
for (base_endpoint in names(accaha_map)) {
  endpoint <- accaha_map[[base_endpoint]]
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, "excess_bmi25_z", "ACC/AHA threshold: burden", "excess_bmi25_z")
  sensitivity_list[[length(sensitivity_list) + 1]] <- fit_rr(endpoint, c("bmi_vim_z", "bmi_mean_z"), "ACC/AHA threshold: VIM", "bmi_vim_z")
}
sensitivity <- do.call(rbind, Filter(Negate(is.null), sensitivity_list))
write.csv(sensitivity, file.path(out, "elsa_sensitivity_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

quartile_list <- list()
change_list <- list()
for (endpoint in endpoints) {
  quartile_list[[length(quartile_list) + 1]] <- fit_rr(endpoint, "burden_q", "Burden quartiles", "burden_q")
  quartile_list[[length(quartile_list) + 1]] <- fit_rr(endpoint, c("vim_q", "bmi_mean_z"), "VIM quartiles", "vim_q")
  change_list[[length(change_list) + 1]] <- fit_rr(endpoint, "change_group", "BMI change category", "change_group")
}
quartiles <- do.call(rbind, Filter(Negate(is.null), quartile_list))
write.csv(quartiles, file.path(out, "elsa_quartile_transition_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")
changes <- do.call(rbind, Filter(Negate(is.null), change_list))
changes$p_fdr <- p.adjust(changes$p_value, method = "BH")
write.csv(changes, file.path(out, "elsa_weight_change_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Prespecified effect-modification tests.
interaction_test <- function(endpoint, exposure, subgroup, model_label) {
  extra <- if (exposure == "bmi_vim_z") "bmi_mean_z" else character(0)
  z <- model_data(endpoint, c(exposure, extra), "ipcw", additional = subgroup)
  if (nrow(z) < 150 || sum(z[[endpoint]]) < 30 || length(unique(z[[subgroup]])) < 2) return(NULL)
  covars <- covars_for(endpoint)
  remove_map <- c(female_f = "female", age65_f = "age_start10", obese_f = "", rural_f = "rural")
  covars <- setdiff(covars, remove_map[[subgroup]])
  des <- svydesign(ids = ~idahhw6, weights = ~analysis_w, data = z)
  formula <- as.formula(paste(endpoint, "~", paste(c(paste0(exposure, "*", subgroup), extra, covars), collapse = " + ")))
  fit <- tryCatch(svyglm(formula, design = des, family = quasipoisson(link = "log")), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  terms <- grep(paste0("^", exposure, ":", subgroup), names(coef(fit)), value = TRUE)
  if (length(terms) != 1) return(NULL)
  term <- terms[1]; beta <- coef(fit)[term]; se <- sqrt(vcov(fit)[term, term])
  data.frame(
    endpoint = endpoint, endpoint_label = endpoint_labels[[endpoint]], exposure = exposure, subgroup = subgroup,
    model = model_label, interaction_log_rr = beta, interaction_ratio = exp(beta),
    ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se),
    p_interaction = 2 * pnorm(abs(beta / se), lower.tail = FALSE), n_intervals = nrow(z), events = sum(z[[endpoint]])
  )
}

interaction_list <- list()
for (endpoint in endpoints) {
  for (subgroup in c("female_f", "age65_f", "obese_f", "rural_f")) {
    interaction_list[[length(interaction_list) + 1]] <- interaction_test(endpoint, "excess_bmi25_z", subgroup, "Burden interaction")
    interaction_list[[length(interaction_list) + 1]] <- interaction_test(endpoint, "bmi_vim_z", subgroup, "VIM interaction")
  }
}
interactions <- do.call(rbind, Filter(Negate(is.null), interaction_list))
write.csv(interactions, file.path(out, "elsa_subgroup_interactions.csv"), row.names = FALSE, fileEncoding = "UTF-8")

decisions <- data.frame(
  item = c("Effect scale", "Primary burden", "Primary variability", "Multiplicity", "Weighting", "Correlation", "Time interpretation"),
  decision = c(
    "Transition risk ratio from household-clustered survey-weighted modified Poisson models",
    "Cumulative excess BMI above 25 kg/m2 over W2-W6, per cohort SD",
    "VIM per cohort SD, adjusted for mean BMI",
    "Benjamini-Hochberg FDR across 10 primary burden/VIM tests",
    "W6 nurse weight x stabilized exposure-history weight x stabilized outcome-observation weight",
    "Household-clustered robust standard errors",
    "Panel-state transitions; W6-W8 and W6-W9 durations are adjusted by interval indicator; TIR sensitivity uses log-time offset"
  )
)
write.csv(decisions, file.path(out, "elsa_statistical_decisions.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Publication figures.
theme_set(theme_classic(base_size = 7, base_family = "sans") +
            theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
                  plot.title = element_text(size = 8, face = "bold"), strip.text = element_text(size = 7, face = "bold"),
                  legend.title = element_text(size = 6.5), legend.text = element_text(size = 6), panel.grid = element_blank()))

transitions <- read.csv(file.path(out, "elsa_transition_counts.csv"), fileEncoding = "UTF-8-BOM")
transitions$total_origin <- ave(transitions$n, transitions$interval, transitions$origin, FUN = sum)
transitions$probability <- ifelse(transitions$total_origin > 0, transitions$n / transitions$total_origin, NA)
transitions$origin_f <- factor(transitions$origin, levels = 5:1,
                               labels = c("Treated uncontrolled", "Treated controlled", "Untreated hypertension", "Elevated untreated", "Normal untreated"))
transitions$destination_f <- factor(transitions$destination, levels = 1:5,
                                    labels = c("Normal untreated", "Elevated untreated", "Untreated hypertension", "Treated controlled", "Treated uncontrolled"))
heatmap <- ggplot(transitions, aes(destination_f, origin_f, fill = probability)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n > 0, sprintf("%d\n%.1f%%", n, 100 * probability), "")), size = 2.0, lineheight = 0.9) +
  facet_wrap(~interval, nrow = 1) +
  scale_fill_gradient(low = "#F3F6F8", high = "#3B6C8E", labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Destination state", y = "Origin state", fill = "Row probability", title = "ELSA observed blood pressure and treatment-state transitions") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

forest <- subset(main, model %in% c("Cumulative excess BMI burden", "BMI variability (VIM)"))
forest$model <- factor(forest$model, levels = c("Cumulative excess BMI burden", "BMI variability (VIM)"),
                       labels = c("Excess BMI burden", "BMI variability (VIM)"))
forest$endpoint_label <- factor(forest$endpoint_label, levels = rev(unname(endpoint_labels)))
forest_plot <- ggplot(forest, aes(rr, endpoint_label, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = "#777777", linewidth = 0.35, linetype = 2) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.18,
                position = position_dodge(width = 0.45), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.45), size = 2.0) +
  scale_x_log10() +
  scale_colour_manual(values = c("#34688A", "#B0684B")) +
  labs(x = "Transition risk ratio per 1 SD (log scale)", y = NULL, colour = NULL, shape = NULL,
       title = "ELSA validation estimates for cumulative BMI burden and variability") +
  theme(legend.position = "top")

save_pub <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  width <- width_mm / 25.4; height <- height_mm / 25.4
  svglite::svglite(file.path(out, paste0(stem, ".svg")), width = width, height = height); print(plot); dev.off()
  grDevices::cairo_pdf(file.path(out, paste0(stem, ".pdf")), width = width, height = height, family = "sans"); print(plot); dev.off()
  ragg::agg_tiff(file.path(out, paste0(stem, ".tiff")), width = width, height = height, units = "in", res = dpi, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(file.path(out, paste0(stem, "_preview.png")), width = width, height = height, units = "in", res = 180); print(plot); dev.off()
}
save_pub(heatmap, "elsa_figure1_transition_heatmap", 183, 105)
save_pub(forest_plot, "elsa_figure2_main_forest", 150, 100)

cat("ELSA main models:\n")
print(main[, c("endpoint", "model", "rr", "ci_low", "ci_high", "p_value", "p_fdr", "n_intervals", "events")])
cat("\nELSA harmonized TIR models:\n")
print(tir[, c("endpoint", "model", "tir", "ci_low", "ci_high", "p_value", "p_fdr", "n_intervals", "events")])
