options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(ggplot2)
  library(MASS)
})

set.seed(20260831)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
output_root <- Sys.getenv("BMI_BP_OUTPUT_DIR", file.path(root, "outputs"))
out <- file.path(output_root, "dynamic_analysis")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

endpoint_labels <- c(
  hypertension_onset = "Hypertension onset",
  untreated_bp_improvement = "Transition to a lower BP state"
)

prepare_data <- function(cohort) {
  file <- file.path(out, paste0(tolower(cohort), "_extended_intervals.csv"))
  data <- read.csv(file, check.names = FALSE, fileEncoding = "UTF-8")
  data$interval <- factor(data$interval)
  data$burden10 <- data$excess_bmi25 / 10
  if (cohort == "HRS") {
    data$race_f <- factor(data$raracem)
    data$region4 <- ifelse(data$region14 %in% c(1, 2), 1,
                           ifelse(data$region14 %in% c(3, 4), 2,
                                  ifelse(data$region14 %in% c(5, 6, 7), 3,
                                         ifelse(data$region14 %in% c(8, 9), 4, NA))))
    data$region_f <- factor(data$region4)
  } else if (cohort == "CHNS") {
    data$province_f <- factor(data$province)
  } else {
    data$education_f <- factor(data$education)
    data$region_f <- factor(data$region)
  }
  data
}

covars_for <- function(cohort) {
  if (cohort == "HRS") {
    return(c("age_start10", "female", "race_f", "raedyrs", "married_partnered", "current_smoker",
             "current_drinker", "r12diabe", "cvd_history", "r12shlt", "region_f", "rural", "interval"))
  }
  if (cohort == "CHNS") {
    return(c("age_start10", "female", "educ_years", "current_smoker", "current_drinker", "diabetes",
             "cvd_history", "self_health", "province_f", "urban", "interval"))
  }
  c("age_start10", "female", "education_f", "married_partnered", "current_smoker", "current_drinker",
    "diabetes", "cvd_history", "self_health", "region_f", "rural", "interval")
}

make_design <- function(data, cohort) {
  if (cohort == "HRS") {
    svydesign(ids = ~cluster, strata = ~strata_value, weights = ~model_weight, data = data, nest = TRUE)
  } else {
    svydesign(ids = ~cluster, weights = ~model_weight, data = data)
  }
}

stack_transitions <- function(data, cohort, exposure) {
  progression <- data[!is.na(data$hypertension_onset), ]
  progression$transition <- "Progression"
  progression$event <- progression$hypertension_onset
  improvement <- data[!is.na(data$untreated_bp_improvement), ]
  improvement$transition <- "Improvement"
  improvement$event <- improvement$untreated_bp_improvement
  z <- rbind(progression, improvement)
  z$transition <- factor(z$transition, levels = c("Progression", "Improvement"))
  covars <- covars_for(cohort)
  needed <- unique(c("event", "model_weight", "cluster", "interval_years", exposure, covars))
  if (cohort == "HRS") needed <- c(needed, "strata_value")
  keep <- z$model_weight > 0 & is.finite(z$model_weight) & complete.cases(z[, needed, drop = FALSE])
  z[keep, ]
}

joint_formula <- function(exposure, covars) {
  # Transition-specific coefficients reproduce the logic of fitting separate
  # risk-set models while retaining their covariance for formal contrasts.
  transition_terms <- paste0("transition:", c(exposure, covars))
  as.formula(paste("event ~ 0 + transition +", paste(transition_terms, collapse = " + "),
                   "+ offset(log(interval_years))"))
}

find_transition_term <- function(coef_names, transition, exposure) {
  candidates <- c(paste0("transition", transition, ":", exposure),
                  paste0(exposure, ":transition", transition))
  hit <- candidates[candidates %in% coef_names]
  if (length(hit) != 1) stop("Could not identify ", transition, " exposure term: ", paste(candidates, collapse = ", "))
  hit
}

contrast_row <- function(beta, variance, label, cohort, exposure, contrast_type) {
  se <- sqrt(max(variance, 0))
  data.frame(
    cohort = cohort, exposure = exposure, contrast = contrast_type, label = label,
    log_ratio = beta, se = se, ratio = exp(beta), ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se), p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE)
  )
}

fit_joint <- function(data, cohort, exposure) {
  z <- stack_transitions(data, cohort, exposure)
  formula <- joint_formula(exposure, covars_for(cohort))
  fit <- svyglm(formula, design = make_design(z, cohort), family = quasibinomial(link = "cloglog"))
  coefficients <- coef(fit)
  covariance <- vcov(fit)
  progression_term <- find_transition_term(names(coefficients), "Progression", exposure)
  improvement_term <- find_transition_term(names(coefficients), "Improvement", exposure)
  terms <- c(progression_term, improvement_term)
  b <- coefficients[terms]
  v <- covariance[terms, terms, drop = FALSE]
  effect_rows <- lapply(seq_along(terms), function(index) {
    se <- sqrt(v[index, index])
    transition <- c("Progression", "Improvement")[index]
    data.frame(
      cohort = cohort, exposure = exposure, transition = transition, term = terms[index],
      log_tir = b[index], se = se, tir = exp(b[index]), ci_low = exp(b[index] - 1.96 * se),
      ci_high = exp(b[index] + 1.96 * se), p_value = 2 * pnorm(abs(b[index] / se), lower.tail = FALSE),
      n_intervals = sum(z$transition == transition), events = sum(z$event[z$transition == transition]),
      participants = length(unique(z$person_id[z$transition == transition]))
    )
  })
  effects <- do.call(rbind, effect_rows)
  difference <- b[1] - b[2]
  difference_variance <- v[1, 1] + v[2, 2] - 2 * v[1, 2]
  sum_effect <- b[1] + b[2]
  sum_variance <- v[1, 1] + v[2, 2] + 2 * v[1, 2]
  contrasts <- rbind(
    contrast_row(difference, difference_variance,
                 "Progression TIR divided by improvement TIR", cohort, exposure, "directional_difference"),
    contrast_row(sum_effect, sum_variance,
                 "Magnitude asymmetry; null means equal and opposite log effects", cohort, exposure, "magnitude_asymmetry")
  )
  list(fit = fit, data = z, effects = effects, contrasts = contrasts,
       terms = terms, covariance = v, coefficients = b)
}

draw_coefficients <- function(fit, simulations = 1200) {
  beta <- coef(fit)
  keep <- is.finite(beta)
  beta <- beta[keep]
  covariance <- vcov(fit)[keep, keep, drop = FALSE]
  covariance <- (covariance + t(covariance)) / 2
  decomposition <- eigen(covariance, symmetric = TRUE)
  values <- pmax(decomposition$values, 0)
  root_matrix <- decomposition$vectors %*% diag(sqrt(values), nrow = length(values))
  noise <- matrix(rnorm(simulations * length(beta)), nrow = simulations) %*% t(root_matrix)
  draws <- sweep(noise, 2, beta, "+")
  list(beta = beta, covariance = covariance, draws = draws)
}

prediction_matrix <- function(fit, newdata, beta_names, horizon_years = 4) {
  terms_no_response <- delete.response(terms(fit))
  matrix <- model.matrix(terms_no_response, newdata)
  matrix <- matrix[, beta_names, drop = FALSE]
  list(matrix = matrix, offset = rep(log(horizon_years), nrow(newdata)))
}

marginal_absolute_probabilities <- function(joint, cohort, exposure_levels = c(0, 10, 20),
                                            horizon_years = 4, simulations = 1200) {
  fit <- joint$fit
  z <- joint$data
  burden_mean <- mean(z$excess_bmi25, na.rm = TRUE)
  burden_sd <- sd(z$excess_bmi25, na.rm = TRUE)
  coefficient_draws <- draw_coefficients(fit, simulations)
  probability_rows <- list()
  draw_store <- list()
  for (transition_value in levels(z$transition)) {
    base <- z[z$transition == transition_value, ]
    weights <- base$model_weight / sum(base$model_weight)
    for (level in exposure_levels) {
      newdata <- base
      newdata$excess_bmi25_z <- (level - burden_mean) / burden_sd
      newdata$burden10 <- level / 10
      newdata$interval_years <- horizon_years
      design <- prediction_matrix(fit, newdata, names(coefficient_draws$beta), horizon_years)
      eta <- design$matrix %*% coefficient_draws$beta + design$offset
      point_probability <- weighted.mean(1 - exp(-exp(eta)), weights)
      eta_draws <- design$matrix %*% t(coefficient_draws$draws) + design$offset
      probability_draws <- colSums((1 - exp(-exp(eta_draws))) * weights)
      key <- paste(transition_value, level, sep = "_")
      draw_store[[key]] <- probability_draws
      probability_rows[[length(probability_rows) + 1]] <- data.frame(
        cohort = cohort, transition = transition_value, burden_bmi_years = level,
        horizon_years = horizon_years, adjusted_probability = point_probability,
        ci_low = quantile(probability_draws, 0.025), ci_high = quantile(probability_draws, 0.975)
      )
    }
  }
  differences <- list()
  for (transition_value in levels(z$transition)) {
    low_key <- paste(transition_value, min(exposure_levels), sep = "_")
    high_key <- paste(transition_value, max(exposure_levels), sep = "_")
    difference_draws <- draw_store[[high_key]] - draw_store[[low_key]]
    point_rows <- do.call(rbind, probability_rows)
    low <- point_rows[point_rows$transition == transition_value & point_rows$burden_bmi_years == min(exposure_levels), "adjusted_probability"]
    high <- point_rows[point_rows$transition == transition_value & point_rows$burden_bmi_years == max(exposure_levels), "adjusted_probability"]
    difference <- high - low
    differences[[length(differences) + 1]] <- data.frame(
      cohort = cohort, transition = transition_value,
      comparison = paste0(max(exposure_levels), " vs ", min(exposure_levels), " BMI-years"),
      risk_difference = difference, ci_low = quantile(difference_draws, 0.025),
      ci_high = quantile(difference_draws, 0.975), events_per_1000 = 1000 * difference,
      events_per_1000_ci_low = 1000 * quantile(difference_draws, 0.025),
      events_per_1000_ci_high = 1000 * quantile(difference_draws, 0.975)
    )
  }
  list(probabilities = do.call(rbind, probability_rows), differences = do.call(rbind, differences))
}

fit_single_endpoint <- function(data, cohort, endpoint, exposure = "excess_bmi25_z", minimum_events = 20) {
  covars <- covars_for(cohort)
  needed <- unique(c(endpoint, exposure, "model_weight", "cluster", "interval_years", covars))
  if (cohort == "HRS") needed <- c(needed, "strata_value")
  keep <- !is.na(data[[endpoint]]) & data$model_weight > 0 & complete.cases(data[, needed, drop = FALSE])
  z <- data[keep, ]
  events <- sum(z[[endpoint]])
  if (nrow(z) < 60 || events < minimum_events || length(unique(z[[endpoint]])) < 2) {
    return(data.frame(cohort = cohort, endpoint = endpoint, exposure = exposure, status = "EVENT_GATE_FAIL",
                      log_tir = NA_real_, se = NA_real_, tir = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                      p_value = NA_real_, n_intervals = nrow(z), events = events,
                      participants = length(unique(z$person_id))))
  }
  formula <- as.formula(paste(endpoint, "~", paste(c(exposure, covars, "offset(log(interval_years))"), collapse = " + ")))
  fit <- svyglm(formula, design = make_design(z, cohort), family = quasibinomial(link = "cloglog"))
  beta <- coef(fit)[exposure]
  se <- sqrt(vcov(fit)[exposure, exposure])
  data.frame(
    cohort = cohort, endpoint = endpoint, exposure = exposure, status = "ESTIMATED",
    log_tir = beta, se = se, tir = exp(beta), ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se), p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
    n_intervals = nrow(z), events = events, participants = length(unique(z$person_id))
  )
}

rcs_basis <- function(x, knots) {
  stopifnot(length(knots) == 4)
  truncated <- function(value, knot) pmax(value - knot, 0)^3
  last <- knots[4]
  penultimate <- knots[3]
  scale <- (last - knots[1])^2
  result <- sapply(1:2, function(index) {
    knot <- knots[index]
    (truncated(x, knot) - ((last - knot) / (last - penultimate)) * truncated(x, penultimate) +
       ((penultimate - knot) / (last - penultimate)) * truncated(x, last)) / scale
  })
  colnames(result) <- c("rcs1", "rcs2")
  result
}

wald_test <- function(fit, terms) {
  available <- terms[terms %in% names(coef(fit))]
  if (length(available) == 0) return(c(statistic = NA, df = 0, p_value = NA))
  beta <- coef(fit)[available]
  variance <- vcov(fit)[available, available, drop = FALSE]
  statistic <- as.numeric(t(beta) %*% MASS::ginv(variance) %*% beta)
  c(statistic = statistic, df = length(available), p_value = pchisq(statistic, df = length(available), lower.tail = FALSE))
}

fit_two_part_dose_response <- function(data, cohort, endpoint) {
  covars <- covars_for(cohort)
  data$any_burden <- as.integer(data$excess_bmi25 > 0)
  positive <- data$burden10[data$burden10 > 0 & is.finite(data$burden10)]
  knots <- as.numeric(quantile(positive, c(0.05, 0.35, 0.65, 0.95), na.rm = TRUE, names = FALSE))
  basis <- rcs_basis(data$burden10, knots)
  data$rcs1 <- basis[, 1]
  data$rcs2 <- basis[, 2]
  terms_of_interest <- c("any_burden", "burden10", "rcs1", "rcs2")
  needed <- unique(c(endpoint, terms_of_interest, "model_weight", "cluster", "interval_years", covars))
  if (cohort == "HRS") needed <- c(needed, "strata_value")
  keep <- !is.na(data[[endpoint]]) & data$model_weight > 0 & complete.cases(data[, needed, drop = FALSE])
  z <- data[keep, ]
  if (sum(z[[endpoint]]) < 20) return(NULL)
  formula <- as.formula(paste(endpoint, "~", paste(c(terms_of_interest, covars, "offset(log(interval_years))"), collapse = " + ")))
  fit <- svyglm(formula, design = make_design(z, cohort), family = quasibinomial(link = "cloglog"))
  nonlinear <- wald_test(fit, c("rcs1", "rcs2"))
  overall <- wald_test(fit, terms_of_interest)
  tests <- data.frame(
    cohort = cohort, endpoint = endpoint, n_intervals = nrow(z), events = sum(z[[endpoint]]),
    knots_bmi_years = paste(round(knots * 10, 2), collapse = ";"),
    p_overall = overall["p_value"], p_nonlinearity = nonlinear["p_value"]
  )
  list(fit = fit, data = z, knots = knots, tests = tests)
}

marginal_dose_curve <- function(result, cohort, endpoint, horizon_years = 4, simulations = 1000) {
  fit <- result$fit
  z <- result$data
  knots <- result$knots
  grid10 <- c(0, seq(knots[1], knots[4], length.out = 31))
  coefficient_draws <- draw_coefficients(fit, simulations)
  weights <- z$model_weight / sum(z$model_weight)
  rows <- list()
  for (level10 in grid10) {
    newdata <- z
    newdata$burden10 <- level10
    newdata$excess_bmi25 <- level10 * 10
    newdata$any_burden <- as.integer(level10 > 0)
    basis <- rcs_basis(rep(level10, nrow(newdata)), knots)
    newdata$rcs1 <- basis[, 1]
    newdata$rcs2 <- basis[, 2]
    newdata$interval_years <- horizon_years
    design <- prediction_matrix(fit, newdata, names(coefficient_draws$beta), horizon_years)
    eta <- design$matrix %*% coefficient_draws$beta + design$offset
    point <- weighted.mean(1 - exp(-exp(eta)), weights)
    eta_draws <- design$matrix %*% t(coefficient_draws$draws) + design$offset
    probability_draws <- colSums((1 - exp(-exp(eta_draws))) * weights)
    rows[[length(rows) + 1]] <- data.frame(
      cohort = cohort, endpoint = endpoint, burden_bmi_years = level10 * 10,
      horizon_years = horizon_years, adjusted_probability = point,
      ci_low = quantile(probability_draws, 0.025), ci_high = quantile(probability_draws, 0.975)
    )
  }
  do.call(rbind, rows)
}

fit_sustained <- function(cohort, endpoint, model = "age_sex") {
  file <- file.path(out, paste0(tolower(cohort), "_sustained_improvement.csv"))
  data <- read.csv(file, check.names = FALSE, fileEncoding = "UTF-8")
  if (cohort == "HRS") {
    data$age10_dynamic <- data$r12agey_e / 10
  } else {
    data$age10_dynamic <- data$age2006 / 10
  }
  covars <- c("age10_dynamic", "female")
  if (model == "full" && cohort == "CHNS") {
    data$province_f <- factor(data$province)
    covars <- c("age10_dynamic", "female", "educ_years", "current_smoker", "current_drinker", "diabetes",
                "cvd_history", "self_health", "province_f", "urban")
  }
  needed <- unique(c(endpoint, "excess_bmi25_z", "model_weight", "cluster", covars))
  if (cohort == "HRS") needed <- c(needed, "strata_value")
  keep <- !is.na(data[[endpoint]]) & data$model_weight > 0 & complete.cases(data[, needed, drop = FALSE])
  z <- data[keep, ]
  events <- sum(z[[endpoint]])
  if (nrow(z) < 30 || events < 15) {
    return(data.frame(cohort = cohort, endpoint = endpoint, model = model, status = "EVENT_GATE_FAIL",
                      log_rr = NA_real_, se = NA_real_, rr = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                      p_value = NA_real_, n = nrow(z), events = events))
  }
  design <- if (cohort == "HRS") {
    svydesign(ids = ~cluster, strata = ~strata_value, weights = ~model_weight, data = z, nest = TRUE)
  } else {
    svydesign(ids = ~cluster, weights = ~model_weight, data = z)
  }
  formula <- as.formula(paste(endpoint, "~ excess_bmi25_z +", paste(covars, collapse = " + ")))
  fit <- svyglm(formula, design = design, family = quasipoisson(link = "log"))
  beta <- coef(fit)["excess_bmi25_z"]
  se <- sqrt(vcov(fit)["excess_bmi25_z", "excess_bmi25_z"])
  data.frame(
    cohort = cohort, endpoint = endpoint, model = model, status = ifelse(events < 20, "EXPLORATORY_SPARSE", "ESTIMATED"),
    log_rr = beta, se = se, rr = exp(beta), ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se),
    p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE), n = nrow(z), events = events
  )
}

datasets <- lapply(c("HRS", "CHNS", "ELSA"), prepare_data)
names(datasets) <- c("HRS", "CHNS", "ELSA")

joint_z <- lapply(names(datasets), function(cohort) fit_joint(datasets[[cohort]], cohort, "excess_bmi25_z"))
names(joint_z) <- names(datasets)
joint_10 <- lapply(names(datasets), function(cohort) fit_joint(datasets[[cohort]], cohort, "burden10"))
names(joint_10) <- names(datasets)

joint_effects <- do.call(rbind, lapply(joint_z, `[[`, "effects"))
joint_contrasts <- do.call(rbind, lapply(joint_z, `[[`, "contrasts"))
per10_effects <- do.call(rbind, lapply(joint_10, `[[`, "effects"))
per10_contrasts <- do.call(rbind, lapply(joint_10, `[[`, "contrasts"))
write.csv(joint_effects, file.path(out, "joint_transition_effects_per_sd.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(joint_contrasts, file.path(out, "joint_transition_contrasts_per_sd.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(per10_effects, file.path(out, "joint_transition_effects_per_10_bmi_years.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(per10_contrasts, file.path(out, "joint_transition_contrasts_per_10_bmi_years.csv"), row.names = FALSE, fileEncoding = "UTF-8")

absolute <- lapply(names(joint_z), function(cohort) marginal_absolute_probabilities(joint_z[[cohort]], cohort))
absolute_probabilities <- do.call(rbind, lapply(absolute, `[[`, "probabilities"))
absolute_differences <- do.call(rbind, lapply(absolute, `[[`, "differences"))
write.csv(absolute_probabilities, file.path(out, "standardized_4year_transition_probabilities.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(absolute_differences, file.path(out, "standardized_4year_risk_differences.csv"), row.names = FALSE, fileEncoding = "UTF-8")

single_endpoints <- c("core3_hypertension_onset", "core3_hypertension_improvement",
                      "buffer5_hypertension_onset", "buffer5_untreated_improvement",
                      "buffer10_hypertension_onset", "buffer10_untreated_improvement")
sensitivity_rows <- list()
for (cohort in names(datasets)) {
  for (endpoint in single_endpoints) {
    sensitivity_rows[[length(sensitivity_rows) + 1]] <- fit_single_endpoint(datasets[[cohort]], cohort, endpoint)
  }
}
sensitivity <- do.call(rbind, sensitivity_rows)
write.csv(sensitivity, file.path(out, "core3_and_buffered_sensitivity_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

dose_results <- list()
for (cohort in names(datasets)) {
  for (endpoint in c("hypertension_onset", "untreated_bp_improvement")) {
    result <- fit_two_part_dose_response(datasets[[cohort]], cohort, endpoint)
    if (!is.null(result)) dose_results[[paste(cohort, endpoint, sep = "_")]] <- result
  }
}
dose_tests <- do.call(rbind, lapply(dose_results, `[[`, "tests"))
dose_tests$p_nonlinearity_fdr <- p.adjust(dose_tests$p_nonlinearity, method = "BH")
dose_tests$p_overall_fdr <- p.adjust(dose_tests$p_overall, method = "BH")
write.csv(dose_tests, file.path(out, "two_part_rcs_tests.csv"), row.names = FALSE, fileEncoding = "UTF-8")
dose_curves <- do.call(rbind, lapply(names(dose_results), function(name) {
  parts <- strsplit(name, "_", fixed = TRUE)[[1]]
  cohort <- parts[1]
  endpoint <- paste(parts[-1], collapse = "_")
  marginal_dose_curve(dose_results[[name]], cohort, endpoint)
}))
write.csv(dose_curves, file.path(out, "two_part_rcs_standardized_curves.csv"), row.names = FALSE, fileEncoding = "UTF-8")

sustained <- rbind(
  fit_sustained("HRS", "initial_improvement", "age_sex"),
  fit_sustained("HRS", "sustained_improvement", "age_sex"),
  fit_sustained("HRS", "maintenance_among_initial_improvers", "age_sex"),
  fit_sustained("CHNS", "initial_improvement", "age_sex"),
  fit_sustained("CHNS", "sustained_improvement", "age_sex"),
  fit_sustained("CHNS", "maintenance_among_initial_improvers", "age_sex"),
  fit_sustained("CHNS", "initial_improvement", "full"),
  fit_sustained("CHNS", "sustained_improvement", "full"),
  fit_sustained("CHNS", "maintenance_among_initial_improvers", "full")
)
write.csv(sustained, file.path(out, "sustained_improvement_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Figure contract: quantitative grid.  The figure defends one conclusion:
# higher cumulative excess BMI burden raises progression probability and lowers
# improvement probability.  Points show marginal 4-year probabilities; the
# lower panel translates 20 versus 0 BMI-years into events per 1000.
transition_labels <- c(Progression = "Hypertension onset", Improvement = "Transition to a lower BP state")
palette <- c(HRS = "#3B6FB6", CHNS = "#D17C2F", ELSA = "#4D9A78")
absolute_probabilities$transition_label <- factor(transition_labels[absolute_probabilities$transition],
                                                  levels = unname(transition_labels))
absolute_differences$transition_label <- factor(transition_labels[absolute_differences$transition],
                                                levels = unname(transition_labels))

p_probability <- ggplot(absolute_probabilities,
                        aes(burden_bmi_years, 100 * adjusted_probability, color = cohort, group = cohort)) +
  geom_ribbon(aes(ymin = 100 * ci_low, ymax = 100 * ci_high, fill = cohort), alpha = 0.10, color = NA) +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.8, shape = 21, fill = "white", stroke = 0.6) +
  facet_wrap(~transition_label, scales = "free_y", nrow = 1) +
  scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
  scale_x_continuous(breaks = c(0, 10, 20)) +
  labs(x = "Cumulative excess BMI burden (BMI-years)", y = "Adjusted 4-year probability (%)",
       color = NULL, fill = NULL) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7),
        legend.position = "top", legend.justification = "left", panel.spacing = unit(10, "pt"))

p_difference <- ggplot(absolute_differences,
                       aes(events_per_1000, cohort, color = cohort)) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.35, linetype = "dashed") +
  geom_errorbarh(aes(xmin = events_per_1000_ci_low, xmax = events_per_1000_ci_high),
                 height = 0.14, linewidth = 0.55) +
  geom_point(size = 2) +
  facet_wrap(~transition_label, scales = "free_x", nrow = 1) +
  scale_color_manual(values = palette, guide = "none") +
  labs(x = "Adjusted difference per 1000 people (20 vs 0 BMI-years)", y = NULL) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7),
        panel.spacing = unit(10, "pt"))

if (requireNamespace("patchwork", quietly = TRUE)) {
  figure <- p_probability / p_difference + patchwork::plot_layout(heights = c(1.55, 1)) +
    patchwork::plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 8)))
} else {
  figure <- p_probability
}

save_publication_figure <- function(plot, stem, width_mm = 183, height_mm = 128, dpi = 600) {
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = width, height = height)
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width, height = height, family = "Arial")
  print(plot); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width, height = height, units = "in", res = dpi)
  print(plot); dev.off()
  ragg::agg_png(paste0(stem, "_preview.png"), width = width, height = height, units = "in", res = 180)
  print(plot); dev.off()
}
save_publication_figure(figure, file.path(out, "figure_dynamic_absolute_transitions"))

dose_curves$endpoint_label <- endpoint_labels[dose_curves$endpoint]
dose_curves$facet_label <- paste(dose_curves$cohort, dose_curves$endpoint_label, sep = " | ")
dose_tests$facet_label <- paste(dose_tests$cohort, endpoint_labels[dose_tests$endpoint], sep = " | ")
dose_annotations <- data.frame(
  facet_label = dose_tests$facet_label,
  label = paste0("P-nonlinearity = ", formatC(dose_tests$p_nonlinearity, digits = 2, format = "f")),
  stringsAsFactors = FALSE
)
p_dose <- ggplot(dose_curves, aes(burden_bmi_years, 100 * adjusted_probability)) +
  geom_ribbon(aes(ymin = 100 * ci_low, ymax = 100 * ci_high), fill = "#8EBAD9", alpha = 0.25) +
  geom_line(color = "#275D8C", linewidth = 0.65) +
  geom_text(data = dose_annotations, aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1.4, size = 2.1, color = "#333333") +
  facet_wrap(~facet_label, scales = "free", ncol = 2) +
  labs(x = "Cumulative excess BMI burden (BMI-years)", y = "Adjusted 4-year probability (%)") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 6.5),
        panel.spacing = unit(7, "pt"))
save_publication_figure(p_dose, file.path(out, "figure_two_part_dose_response"), height_mm = 150)

figure_contract <- data.frame(
  item = c("core_conclusion", "evidence_panel_a", "evidence_panel_b", "archetype", "backend", "export_contract", "review_risks"),
  value = c(
    "Higher cumulative excess BMI burden raises hypertension-onset probability and lowers probability of transition to a lower BP state.",
    "Marginal standardized 4-year probabilities at 0, 10, and 20 BMI-years in each cohort.",
    "Adjusted event differences per 1000 people comparing 20 with 0 BMI-years.",
    "Quantitative grid with a probability hero panel and a clinical-difference support panel.",
    "R/ggplot2 exclusively.",
    "183 x 128 mm; editable SVG/PDF; 600-dpi TIFF; PNG preview; source CSVs retained.",
    "Model-based 4-year standardization; structural zero burden; sparse improvement risk set; observational interpretation."
  )
)
write.csv(figure_contract, file.path(out, "figure_dynamic_absolute_transitions_contract.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("\nJOINT EFFECTS PER SD\n")
print(joint_effects[, c("cohort", "transition", "tir", "ci_low", "ci_high", "p_value", "events")])
cat("\nJOINT DIRECTIONAL CONTRASTS\n")
print(joint_contrasts[joint_contrasts$contrast == "directional_difference", ])
cat("\nABSOLUTE DIFFERENCES PER 1000\n")
print(absolute_differences)
cat("\nDOSE-RESPONSE TESTS\n")
print(dose_tests)
cat("\nSUSTAINED IMPROVEMENT\n")
print(sustained)
