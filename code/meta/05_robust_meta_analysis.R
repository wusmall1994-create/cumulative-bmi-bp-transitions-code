options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
output_root <- Sys.getenv("BMI_BP_OUTPUT_DIR", file.path(root, "outputs"))
dynamic_out <- file.path(output_root, "dynamic_analysis")
primary_out <- file.path(output_root, "meta_three_cohort")

reml_tau2 <- function(yi, vi) {
  objective <- function(tau2) {
    weights <- 1 / (vi + tau2)
    mean <- sum(weights * yi) / sum(weights)
    sum(log(vi + tau2)) + log(sum(weights)) + sum(weights * (yi - mean)^2)
  }
  upper <- max(c(var(yi, na.rm = TRUE) * 20, max(vi) * 100, 1))
  result <- optimize(objective, interval = c(0, upper), tol = 1e-12)
  if (result$minimum < 1e-8) 0 else result$minimum
}

meta_estimates <- function(data, yi = "yi", sei = "sei", transform = "exp") {
  y <- data[[yi]]
  se <- data[[sei]]
  keep <- is.finite(y) & is.finite(se) & se > 0
  y <- y[keep]
  se <- se[keep]
  k <- length(y)
  if (k < 2) return(NULL)
  vi <- se^2
  fixed_weights <- 1 / vi
  fixed_mean <- sum(fixed_weights * y) / sum(fixed_weights)
  fixed_se <- sqrt(1 / sum(fixed_weights))
  q_fixed <- sum(fixed_weights * (y - fixed_mean)^2)
  i2 <- ifelse(q_fixed > 0, max(0, (q_fixed - (k - 1)) / q_fixed) * 100, 0)
  tau2 <- reml_tau2(y, vi)
  random_weights <- 1 / (vi + tau2)
  random_mean <- sum(random_weights * y) / sum(random_weights)
  random_se <- sqrt(1 / sum(random_weights))
  q_random <- sum(random_weights * (y - random_mean)^2)
  hk_scale <- q_random / (k - 1)
  hk_se <- sqrt(hk_scale / sum(random_weights))
  modified_hk_se <- sqrt(max(1, hk_scale) / sum(random_weights))
  t_critical <- qt(0.975, df = k - 1)
  normal_critical <- qnorm(0.975)
  prediction_critical <- ifelse(k > 2, qt(0.975, df = k - 2), NA_real_)

  convert <- if (transform == "exp") exp else identity
  data.frame(
    k = k,
    fixed_estimate = convert(fixed_mean),
    fixed_ci_low = convert(fixed_mean - normal_critical * fixed_se),
    fixed_ci_high = convert(fixed_mean + normal_critical * fixed_se),
    fixed_p = 2 * pnorm(abs(fixed_mean / fixed_se), lower.tail = FALSE),
    reml_estimate = convert(random_mean),
    reml_ci_low = convert(random_mean - normal_critical * random_se),
    reml_ci_high = convert(random_mean + normal_critical * random_se),
    reml_p = 2 * pnorm(abs(random_mean / random_se), lower.tail = FALSE),
    hk_ci_low = convert(random_mean - t_critical * hk_se),
    hk_ci_high = convert(random_mean + t_critical * hk_se),
    hk_p = 2 * pt(abs(random_mean / hk_se), df = k - 1, lower.tail = FALSE),
    modified_hk_ci_low = convert(random_mean - t_critical * modified_hk_se),
    modified_hk_ci_high = convert(random_mean + t_critical * modified_hk_se),
    modified_hk_p = 2 * pt(abs(random_mean / modified_hk_se), df = k - 1, lower.tail = FALSE),
    prediction_interval_low = ifelse(k > 2, convert(random_mean - prediction_critical * sqrt(tau2 + modified_hk_se^2)), NA_real_),
    prediction_interval_high = ifelse(k > 2, convert(random_mean + prediction_critical * sqrt(tau2 + modified_hk_se^2)), NA_real_),
    tau2_reml = tau2,
    q = q_fixed,
    q_p = pchisq(q_fixed, df = k - 1, lower.tail = FALSE),
    i2_percent = i2,
    hk_scale = hk_scale,
    prediction_interval_caution = ifelse(k <= 3, "VERY_UNCERTAIN_WITH_K_LE_3", "STANDARD_CAUTION")
  )
}

run_grouped_meta <- function(data, group_columns, yi, sei, transform = "exp", source) {
  keys <- interaction(data[, group_columns, drop = FALSE], drop = TRUE, lex.order = TRUE)
  groups <- split(data, keys)
  rows <- list()
  for (subset in groups) {
    estimate <- meta_estimates(subset, yi, sei, transform)
    if (is.null(estimate)) next
    identifiers <- subset[1, group_columns, drop = FALSE]
    identifiers$source <- source
    rows[[length(rows) + 1]] <- cbind(identifiers, estimate)
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

leave_one_out <- function(data, group_columns, cohort_column, yi, sei, transform = "exp", source) {
  keys <- interaction(data[, group_columns, drop = FALSE], drop = TRUE, lex.order = TRUE)
  groups <- split(data, keys)
  rows <- list()
  for (subset in groups) {
    if (nrow(subset) < 3) next
    for (cohort in unique(subset[[cohort_column]])) {
      reduced <- subset[subset[[cohort_column]] != cohort, ]
      estimate <- meta_estimates(reduced, yi, sei, transform)
      identifiers <- subset[1, group_columns, drop = FALSE]
      identifiers$source <- source
      identifiers$omitted_cohort <- cohort
      rows[[length(rows) + 1]] <- cbind(identifiers, estimate)
    }
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

weight_influence <- function(data, group_columns, cohort_column, yi, sei, source) {
  keys <- interaction(data[, group_columns, drop = FALSE], drop = TRUE, lex.order = TRUE)
  groups <- split(data, keys)
  rows <- list()
  for (subset in groups) {
    keep <- is.finite(subset[[yi]]) & is.finite(subset[[sei]]) & subset[[sei]] > 0
    subset <- subset[keep, ]
    if (nrow(subset) < 2) next
    variance <- subset[[sei]]^2
    tau2 <- reml_tau2(subset[[yi]], variance)
    weights <- 1 / (variance + tau2)
    pooled <- sum(weights * subset[[yi]]) / sum(weights)
    standardized_residual <- (subset[[yi]] - pooled) / sqrt(variance + tau2)
    identifiers <- subset[, c(cohort_column, group_columns), drop = FALSE]
    identifiers$source <- source
    identifiers$random_weight_percent <- 100 * weights / sum(weights)
    identifiers$standardized_residual <- standardized_residual
    identifiers$absolute_standardized_residual <- abs(standardized_residual)
    rows[[length(rows) + 1]] <- identifiers
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

bind_fill <- function(frames) {
  frames <- Filter(Negate(is.null), frames)
  columns <- unique(unlist(lapply(frames, names)))
  frames <- lapply(frames, function(frame) {
    missing <- setdiff(columns, names(frame))
    for (column in missing) frame[[column]] <- NA
    frame[, columns, drop = FALSE]
  })
  do.call(rbind, frames)
}

primary <- read.csv(file.path(primary_out, "three_cohort_harmonized_tir.csv"), check.names = FALSE)
primary$yi <- primary$log_tir
primary$sei <- primary$se
primary_meta <- run_grouped_meta(primary, c("endpoint", "endpoint_label", "model"), "yi", "sei", "exp", "primary_transition_models")
primary_loo <- leave_one_out(primary, c("endpoint", "endpoint_label", "model"), "cohort", "yi", "sei", "exp", "primary_transition_models")
primary_influence <- weight_influence(primary, c("endpoint", "endpoint_label", "model"), "cohort", "yi", "sei", "primary_transition_models")

contrasts <- read.csv(file.path(dynamic_out, "joint_transition_contrasts_per_sd.csv"), check.names = FALSE)
contrasts <- contrasts[contrasts$contrast == "directional_difference", ]
contrast_meta <- run_grouped_meta(contrasts, c("contrast", "label", "exposure"), "log_ratio", "se", "exp", "joint_directional_contrast")
contrast_loo <- leave_one_out(contrasts, c("contrast", "label", "exposure"), "cohort", "log_ratio", "se", "exp", "joint_directional_contrast")
contrast_influence <- weight_influence(contrasts, c("contrast", "label", "exposure"), "cohort", "log_ratio", "se", "joint_directional_contrast")

per10 <- read.csv(file.path(dynamic_out, "joint_transition_effects_per_10_bmi_years.csv"), check.names = FALSE)
per10_meta <- run_grouped_meta(per10, c("transition", "exposure"), "log_tir", "se", "exp", "per_10_bmi_years")
per10_loo <- leave_one_out(per10, c("transition", "exposure"), "cohort", "log_tir", "se", "exp", "per_10_bmi_years")
per10_influence <- weight_influence(per10, c("transition", "exposure"), "cohort", "log_tir", "se", "per_10_bmi_years")

sensitivity <- read.csv(file.path(dynamic_out, "core3_and_buffered_sensitivity_models.csv"), check.names = FALSE)
sensitivity <- sensitivity[sensitivity$status == "ESTIMATED", ]
sensitivity_meta <- run_grouped_meta(sensitivity, c("endpoint", "exposure"), "log_tir", "se", "exp", "core3_and_buffered")
sensitivity_loo <- leave_one_out(sensitivity, c("endpoint", "exposure"), "cohort", "log_tir", "se", "exp", "core3_and_buffered")
sensitivity_influence <- weight_influence(sensitivity, c("endpoint", "exposure"), "cohort", "log_tir", "se", "core3_and_buffered")

absolute <- read.csv(file.path(dynamic_out, "standardized_4year_risk_differences.csv"), check.names = FALSE)
absolute$rd_se <- (absolute$ci_high - absolute$ci_low) / (2 * qnorm(0.975))
absolute_meta <- run_grouped_meta(absolute, c("transition", "comparison"), "risk_difference", "rd_se", "identity", "absolute_risk_difference")
absolute_loo <- leave_one_out(absolute, c("transition", "comparison"), "cohort", "risk_difference", "rd_se", "identity", "absolute_risk_difference")
absolute_influence <- weight_influence(absolute, c("transition", "comparison"), "cohort", "risk_difference", "rd_se", "absolute_risk_difference")

all_meta <- bind_fill(list(primary_meta, contrast_meta, per10_meta, sensitivity_meta, absolute_meta))
all_loo <- bind_fill(list(primary_loo, contrast_loo, per10_loo, sensitivity_loo, absolute_loo))
all_influence <- bind_fill(list(primary_influence, contrast_influence, per10_influence,
                                sensitivity_influence, absolute_influence))

write.csv(primary_meta, file.path(dynamic_out, "advanced_meta_primary_results.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(all_meta, file.path(dynamic_out, "advanced_meta_all_results.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(all_loo, file.path(dynamic_out, "advanced_meta_leave_one_out.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(all_influence, file.path(dynamic_out, "advanced_meta_weights_influence.csv"), row.names = FALSE, fileEncoding = "UTF-8")

method_notes <- data.frame(
  item = c("tau2", "hartung_knapp", "modified_hartung_knapp", "prediction_interval", "interpretation"),
  decision = c(
    "Restricted maximum likelihood (REML) for the intercept-only meta-analysis.",
    "Conventional Hartung-Knapp uses the residual scale q/(k-1) and t(k-1) critical values.",
    "Modified Hartung-Knapp constrains the residual scale to at least 1 to prevent paradoxically narrower intervals.",
    "Uses t(k-2); with k=3 it is intentionally very uncertain and is exploratory only.",
    "Report 'no clear heterogeneity detected'; do not equate I2=0 with homogeneity."
  )
)
write.csv(method_notes, file.path(dynamic_out, "advanced_meta_method_notes.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("\nPRIMARY REML + MODIFIED HARTUNG-KNAPP\n")
print(primary_meta[, c("endpoint", "model", "reml_estimate", "reml_ci_low", "reml_ci_high",
                       "modified_hk_ci_low", "modified_hk_ci_high", "modified_hk_p", "tau2_reml", "i2_percent")])
cat("\nJOINT DIRECTIONAL CONTRAST\n")
print(contrast_meta)
cat("\nPER 10 BMI-YEARS\n")
print(per10_meta[, c("transition", "reml_estimate", "reml_ci_low", "reml_ci_high",
                     "modified_hk_ci_low", "modified_hk_ci_high", "modified_hk_p")])
