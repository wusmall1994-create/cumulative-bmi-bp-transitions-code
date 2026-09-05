options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
  library(MASS)
})

set.seed(20260902)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
output_root <- Sys.getenv("BMI_BP_OUTPUT_DIR", file.path(root, "outputs"))
source_dir <- file.path(output_root, "dynamic_analysis")
out <- file.path(output_root, "methodological_sensitivity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

endpoint_labels <- c(
  hypertension_onset = "Hypertension onset",
  untreated_bp_improvement = "Transition to a lower untreated BP state"
)

zscore <- function(x) {
  value <- as.numeric(x)
  (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)
}

effective_n <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (!length(w)) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

prepare_data <- function(cohort) {
  path <- file.path(source_dir, paste0(tolower(cohort), "_extended_intervals.csv"))
  d <- read.csv(path, check.names = FALSE, fileEncoding = "UTF-8")
  d$interval <- factor(d$interval)
  d$cohort <- cohort
  d$person_cluster <- factor(d$person_id)
  d$start_sbp10 <- (d$start_sbp - 130) / 10
  d$start_dbp10 <- (d$start_dbp - 80) / 10
  d$delta_sbp <- d$end_sbp - d$start_sbp
  d$delta_dbp <- d$end_dbp - d$start_dbp
  if (cohort == "HRS") {
    d$race_f <- factor(d$raracem)
    d$region4 <- ifelse(d$region14 %in% c(1, 2), 1,
                        ifelse(d$region14 %in% c(3, 4), 2,
                               ifelse(d$region14 %in% c(5, 6, 7), 3,
                                      ifelse(d$region14 %in% c(8, 9), 4, NA))))
    d$region_f <- factor(d$region4)
    d$psu_id <- interaction(d$strata_value, d$secu, drop = TRUE)
  } else if (cohort == "CHNS") {
    d$province_f <- factor(d$province)
    d$psu_id <- factor(d$commid)
    baseline <- read.csv(file.path(output_root, "chns_analysis", "chns_baseline_2006.csv"),
                         check.names = FALSE, fileEncoding = "UTF-8")
    burden_at <- function(threshold) {
      e00 <- pmax(baseline$bmi_2000 - threshold, 0)
      e04 <- pmax(baseline$bmi_2004 - threshold, 0)
      e06 <- pmax(baseline$bmi_2006 - threshold, 0)
      (e00 + e04) / 2 * 4 + (e04 + e06) / 2 * 2
    }
    baseline$excess_bmi24 <- burden_at(24)
    baseline$excess_bmi23 <- burden_at(23)
    threshold_map <- baseline[, c("person_id", "excess_bmi24", "excess_bmi23")]
    d <- merge(d, threshold_map, by = "person_id", all.x = TRUE, sort = FALSE)
    d$excess_bmi24_z <- zscore(d$excess_bmi24)
    d$excess_bmi23_z <- zscore(d$excess_bmi23)
  } else {
    d$education_f <- factor(d$education)
    d$region_f <- factor(d$region)
    d$psu_id <- factor(d$cluster)
  }
  d
}

minimal_covars <- function(cohort) {
  if (cohort == "HRS") {
    return(c("age_start10", "female", "race_f", "raedyrs", "married_partnered",
             "current_smoker", "current_drinker", "region_f", "rural", "interval"))
  }
  if (cohort == "CHNS") {
    return(c("age_start10", "female", "educ_years", "current_smoker", "current_drinker",
             "province_f", "urban", "interval"))
  }
  c("age_start10", "female", "education_f", "married_partnered", "current_smoker",
    "current_drinker", "region_f", "rural", "interval")
}

health_covars <- function(cohort) {
  if (cohort == "HRS") return(c("r12diabe", "cvd_history", "r12shlt"))
  c("diabetes", "cvd_history", "self_health")
}

make_design <- function(d, cohort, weights = "model_weight", cluster_mode = "sampling") {
  if (weights == "none") {
    d$analysis_weight <- 1
  } else {
    d$analysis_weight <- d[[weights]]
  }
  if (cluster_mode == "person") {
    return(svydesign(ids = ~person_cluster, weights = ~analysis_weight, data = d))
  }
  if (cohort == "HRS") {
    return(svydesign(ids = ~psu_id, strata = ~strata_value, weights = ~analysis_weight,
                     data = d, nest = TRUE))
  }
  svydesign(ids = ~psu_id, weights = ~analysis_weight, data = d)
}

complete_model_data <- function(d, cohort, endpoints, exposure, covars,
                                weights = "model_weight", extra = character(0)) {
  needed <- unique(c(endpoints, exposure, covars, "interval_years", "person_id", "psu_id", extra))
  if (weights != "none") needed <- c(needed, weights)
  if (cohort == "HRS") needed <- c(needed, "strata_value")
  keep <- complete.cases(d[, needed, drop = FALSE])
  if (weights != "none") keep <- keep & is.finite(d[[weights]]) & d[[weights]] > 0
  d[keep, , drop = FALSE]
}

fit_endpoint <- function(d, cohort, endpoint, exposure = "excess_bmi25_z",
                         covar_set = "minimal", bp_adjust = FALSE,
                         weights = "model_weight", cluster_mode = "sampling") {
  covars <- minimal_covars(cohort)
  if (covar_set == "expanded") covars <- c(covars, health_covars(cohort))
  if (bp_adjust) covars <- c(covars, "start_sbp10", "start_dbp10")
  z <- d[!is.na(d[[endpoint]]), , drop = FALSE]
  z <- complete_model_data(z, cohort, endpoint, exposure, covars, weights)
  z <- droplevels(z)
  covars <- covars[vapply(covars, function(variable) {
    !is.factor(z[[variable]]) || nlevels(z[[variable]]) > 1
  }, logical(1))]
  if (nrow(z) < 60 || sum(z[[endpoint]]) < 20 || length(unique(z[[endpoint]])) < 2) return(NULL)
  rhs <- c(exposure, covars, "offset(log(interval_years))")
  form <- as.formula(paste(endpoint, "~", paste(rhs, collapse = " + ")))
  fit <- svyglm(form, design = make_design(z, cohort, weights, cluster_mode),
                family = quasibinomial(link = "cloglog"))
  beta <- coef(fit)[exposure]
  se <- sqrt(vcov(fit)[exposure, exposure])
  data.frame(
    cohort = cohort, endpoint = endpoint, endpoint_label = endpoint_labels[[endpoint]],
    model = paste(covar_set, ifelse(bp_adjust, "plus_baseline_bp", ""),
                  ifelse(cluster_mode == "person", "person_clustered", "sampling_clustered"), sep = ";"),
    exposure = exposure, log_tir = beta, se = se, tir = exp(beta),
    ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se),
    p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
    n_intervals = nrow(z), events = sum(z[[endpoint]]), participants = length(unique(z$person_id)),
    sampling_clusters = length(unique(z$psu_id)), person_clusters = length(unique(z$person_id))
  )
}

joint_formula <- function(exposure, covars) {
  transition_terms <- paste0("transition:", c(exposure, covars))
  as.formula(paste("event ~ 0 + transition +", paste(transition_terms, collapse = " + "),
                   "+ offset(log(interval_years))"))
}

find_transition_term <- function(names_vector, transition, exposure) {
  candidates <- c(paste0("transition", transition, ":", exposure),
                  paste0(exposure, ":transition", transition))
  hit <- candidates[candidates %in% names_vector]
  if (length(hit) != 1) stop("Unable to identify transition-specific exposure term")
  hit
}

fit_joint <- function(d, cohort, exposure = "excess_bmi25_z", covar_set = "minimal",
                      bp_adjust = FALSE, weights = "model_weight", cluster_mode = "sampling") {
  progression <- d[!is.na(d$hypertension_onset), , drop = FALSE]
  progression$transition <- "Progression"
  progression$event <- progression$hypertension_onset
  improvement <- d[!is.na(d$untreated_bp_improvement), , drop = FALSE]
  improvement$transition <- "Improvement"
  improvement$event <- improvement$untreated_bp_improvement
  z <- rbind(progression, improvement)
  z$transition <- factor(z$transition, levels = c("Progression", "Improvement"))
  covars <- minimal_covars(cohort)
  if (covar_set == "expanded") covars <- c(covars, health_covars(cohort))
  if (bp_adjust) covars <- c(covars, "start_sbp10", "start_dbp10")
  z <- complete_model_data(z, cohort, "event", exposure, covars, weights)
  z <- droplevels(z)
  covars <- covars[vapply(covars, function(variable) {
    !is.factor(z[[variable]]) || nlevels(z[[variable]]) > 1
  }, logical(1))]
  fit <- svyglm(joint_formula(exposure, covars), design = make_design(z, cohort, weights, cluster_mode),
                family = quasibinomial(link = "cloglog"))
  coefs <- coef(fit)
  covariance <- vcov(fit)
  pterm <- find_transition_term(names(coefs), "Progression", exposure)
  iterm <- find_transition_term(names(coefs), "Improvement", exposure)
  terms <- c(pterm, iterm)
  b <- coefs[terms]
  v <- covariance[terms, terms, drop = FALSE]
  joint_stat <- as.numeric(t(b) %*% MASS::ginv(v) %*% b)
  diff <- b[1] - b[2]
  diff_var <- v[1, 1] + v[2, 2] - 2 * v[1, 2]
  effect_rows <- lapply(1:2, function(index) {
    transition <- c("Progression", "Improvement")[index]
    se <- sqrt(v[index, index])
    data.frame(
      cohort = cohort, covariate_set = covar_set, bp_adjust = bp_adjust,
      cluster_mode = cluster_mode, exposure = exposure, transition = transition,
      log_tir = b[index], se = se, tir = exp(b[index]),
      ci_low = exp(b[index] - 1.96 * se), ci_high = exp(b[index] + 1.96 * se),
      p_value = 2 * pnorm(abs(b[index] / se), lower.tail = FALSE),
      n_intervals = sum(z$transition == transition), events = sum(z$event[z$transition == transition]),
      participants = length(unique(z$person_id[z$transition == transition])),
      sampling_clusters = length(unique(z$psu_id[z$transition == transition]))
    )
  })
  tests <- data.frame(
    cohort = cohort, covariate_set = covar_set, bp_adjust = bp_adjust,
    cluster_mode = cluster_mode, exposure = exposure,
    joint_wald_chisq = joint_stat, joint_wald_df = 2,
    joint_wald_p = pchisq(joint_stat, df = 2, lower.tail = FALSE),
    directional_log_ratio = diff, directional_se = sqrt(diff_var),
    directional_ratio = exp(diff),
    directional_ci_low = exp(diff - 1.96 * sqrt(diff_var)),
    directional_ci_high = exp(diff + 1.96 * sqrt(diff_var)),
    directional_p = 2 * pnorm(abs(diff / sqrt(diff_var)), lower.tail = FALSE),
    covariance_progression_improvement = v[1, 2],
    total_participants = length(unique(z$person_id)), total_sampling_clusters = length(unique(z$psu_id))
  )
  list(fit = fit, data = z, effects = do.call(rbind, effect_rows), tests = tests,
       coefficient_names = names(coefs), coefficients = coefs, covariance = covariance)
}

draw_coefficients <- function(fit, simulations = 1500) {
  beta <- coef(fit)
  keep <- is.finite(beta)
  beta <- beta[keep]
  covariance <- vcov(fit)[keep, keep, drop = FALSE]
  covariance <- (covariance + t(covariance)) / 2
  eig <- eigen(covariance, symmetric = TRUE)
  root_matrix <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = length(eig$values))
  noise <- matrix(rnorm(simulations * length(beta)), nrow = simulations) %*% t(root_matrix)
  list(beta = beta, draws = sweep(noise, 2, beta, "+"))
}

marginal_percentile_probabilities <- function(joint, cohort, horizon = 4, simulations = 1500) {
  z <- joint$data
  person_exposure <- aggregate(excess_bmi25 ~ person_id, data = z, FUN = function(x) x[1])
  levels <- as.numeric(quantile(person_exposure$excess_bmi25, c(0.25, 0.75), na.rm = TRUE, names = FALSE))
  names(levels) <- c("P25", "P75")
  exposure_mean <- mean(z$excess_bmi25, na.rm = TRUE)
  exposure_sd <- sd(z$excess_bmi25, na.rm = TRUE)
  coefficient_draws <- draw_coefficients(joint$fit, simulations)
  rows <- list()
  difference_rows <- list()
  draw_store <- list()
  terms_no_response <- delete.response(terms(joint$fit))
  for (transition in levels(z$transition)) {
    base <- z[z$transition == transition, , drop = FALSE]
    weights <- base$model_weight / sum(base$model_weight)
    for (label in names(levels)) {
      value <- levels[[label]]
      newdata <- base
      newdata$excess_bmi25 <- value
      newdata$excess_bmi25_z <- (value - exposure_mean) / exposure_sd
      newdata$interval_years <- horizon
      matrix <- model.matrix(terms_no_response, newdata)
      matrix <- matrix[, names(coefficient_draws$beta), drop = FALSE]
      eta <- as.vector(matrix %*% coefficient_draws$beta) + log(horizon)
      probability <- weighted.mean(1 - exp(-exp(eta)), weights)
      eta_draws <- matrix %*% t(coefficient_draws$draws) + log(horizon)
      probability_draws <- colSums((1 - exp(-exp(eta_draws))) * weights)
      draw_store[[paste(transition, label, sep = "_")]] <- probability_draws
      rows[[length(rows) + 1]] <- data.frame(
        cohort = cohort, transition = transition, percentile = label,
        burden_bmi_years = value, horizon_years = horizon,
        adjusted_probability = probability,
        ci_low = quantile(probability_draws, 0.025), ci_high = quantile(probability_draws, 0.975)
      )
    }
    diff_draws <- draw_store[[paste(transition, "P75", sep = "_")]] -
      draw_store[[paste(transition, "P25", sep = "_")]]
    p_rows <- do.call(rbind, rows)
    low <- tail(p_rows$adjusted_probability[p_rows$cohort == cohort & p_rows$transition == transition & p_rows$percentile == "P25"], 1)
    high <- tail(p_rows$adjusted_probability[p_rows$cohort == cohort & p_rows$transition == transition & p_rows$percentile == "P75"], 1)
    difference_rows[[length(difference_rows) + 1]] <- data.frame(
      cohort = cohort, transition = transition,
      p25_bmi_years = levels[["P25"]], p75_bmi_years = levels[["P75"]],
      risk_difference = high - low,
      ci_low = quantile(diff_draws, 0.025), ci_high = quantile(diff_draws, 0.975),
      events_per_1000 = 1000 * (high - low),
      events_per_1000_ci_low = 1000 * quantile(diff_draws, 0.025),
      events_per_1000_ci_high = 1000 * quantile(diff_draws, 0.975)
    )
  }
  list(probabilities = do.call(rbind, rows), differences = do.call(rbind, difference_rows))
}

fit_continuous_change <- function(d, cohort, outcome, exposure = "excess_bmi25_z") {
  covars <- minimal_covars(cohort)
  baseline <- if (outcome == "delta_sbp") "start_sbp10" else "start_dbp10"
  z <- complete_model_data(d, cohort, outcome, exposure, c(covars, baseline), "model_weight")
  form <- as.formula(paste(outcome, "~", paste(c(exposure, baseline, covars), collapse = " + ")))
  fit <- svyglm(form, design = make_design(z, cohort), family = gaussian())
  beta <- coef(fit)[exposure]
  se <- sqrt(vcov(fit)[exposure, exposure])
  data.frame(
    cohort = cohort, outcome = outcome, exposure = exposure, beta_mmHg = beta, se = se,
    ci_low = beta - 1.96 * se, ci_high = beta + 1.96 * se,
    p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
    n_intervals = nrow(z), participants = length(unique(z$person_id)), sampling_clusters = length(unique(z$psu_id))
  )
}

weight_sensitivity <- function(d, cohort, endpoint) {
  z <- d[!is.na(d[[endpoint]]) & is.finite(d$model_weight) & d$model_weight > 0, , drop = FALSE]
  q01 <- quantile(z$model_weight, c(0.01, 0.99), na.rm = TRUE)
  q025 <- quantile(z$model_weight, c(0.025, 0.975), na.rm = TRUE)
  z$weight_trim_1_99 <- pmin(pmax(z$model_weight, q01[1]), q01[2])
  z$weight_trim_2_5_97_5 <- pmin(pmax(z$model_weight, q025[1]), q025[2])
  rows <- list()
  for (spec in list(c("Primary combined weight", "model_weight"),
                    c("Combined weight trimmed 1st-99th", "weight_trim_1_99"),
                    c("Combined weight trimmed 2.5th-97.5th", "weight_trim_2_5_97_5"),
                    c("Unweighted", "none"))) {
    row <- fit_endpoint(z, cohort, endpoint, covar_set = "minimal", weights = spec[2])
    if (!is.null(row)) {
      row$weight_specification <- spec[1]
      rows[[length(rows) + 1]] <- row
    }
  }
  do.call(rbind, rows)
}

weight_distribution <- function(d, cohort) {
  rows <- list()
  for (endpoint in names(endpoint_labels)) {
    z <- d[!is.na(d[[endpoint]]) & is.finite(d$model_weight) & d$model_weight > 0, , drop = FALSE]
    qs <- quantile(z$model_weight, c(0, 0.01, 0.5, 0.99, 1), names = FALSE)
    rows[[length(rows) + 1]] <- data.frame(
      cohort = cohort, endpoint = endpoint, n = nrow(z),
      weight_min = qs[1], weight_p01 = qs[2], weight_median = qs[3], weight_p99 = qs[4], weight_max = qs[5],
      weight_mean = mean(z$model_weight), weight_cv = sd(z$model_weight) / mean(z$model_weight),
      effective_n = effective_n(z$model_weight), effective_n_fraction = effective_n(z$model_weight) / nrow(z)
    )
  }
  do.call(rbind, rows)
}

reml_tau2 <- function(yi, vi) {
  objective <- function(tau2) {
    weights <- 1 / (vi + tau2)
    mean_value <- sum(weights * yi) / sum(weights)
    sum(log(vi + tau2)) + log(sum(weights)) + sum(weights * (yi - mean_value)^2)
  }
  upper <- max(c(var(yi, na.rm = TRUE) * 20, max(vi) * 100, 1))
  result <- optimize(objective, interval = c(0, upper), tol = 1e-12)
  if (result$minimum < 1e-8) 0 else result$minimum
}

meta_estimates <- function(d) {
  y <- d$log_tir
  vi <- d$se^2
  k <- length(y)
  fw <- 1 / vi
  fixed_mean <- sum(fw * y) / sum(fw)
  fixed_se <- sqrt(1 / sum(fw))
  q <- sum(fw * (y - fixed_mean)^2)
  i2 <- ifelse(q > 0, max(0, (q - (k - 1)) / q) * 100, 0)
  tau2 <- reml_tau2(y, vi)
  rw <- 1 / (vi + tau2)
  random_mean <- sum(rw * y) / sum(rw)
  random_se <- sqrt(1 / sum(rw))
  scale <- sum(rw * (y - random_mean)^2) / (k - 1)
  mhk_se <- sqrt(max(1, scale) / sum(rw))
  crit_t <- qt(0.975, k - 1)
  data.frame(
    k = k,
    fixed_tir = exp(fixed_mean), fixed_ci_low = exp(fixed_mean - 1.96 * fixed_se),
    fixed_ci_high = exp(fixed_mean + 1.96 * fixed_se), fixed_p = 2 * pnorm(abs(fixed_mean / fixed_se), lower.tail = FALSE),
    reml_tir = exp(random_mean), reml_ci_low = exp(random_mean - 1.96 * random_se),
    reml_ci_high = exp(random_mean + 1.96 * random_se), reml_p = 2 * pnorm(abs(random_mean / random_se), lower.tail = FALSE),
    modified_hk_ci_low = exp(random_mean - crit_t * mhk_se),
    modified_hk_ci_high = exp(random_mean + crit_t * mhk_se),
    modified_hk_p = 2 * pt(abs(random_mean / mhk_se), df = k - 1, lower.tail = FALSE),
    tau2_reml = tau2, q = q, q_p = pchisq(q, k - 1, lower.tail = FALSE), i2_percent = i2
  )
}

meta_by_transition <- function(effects) {
  rows <- list()
  for (transition in unique(effects$transition)) {
    z <- effects[effects$transition == transition, ]
    row <- meta_estimates(z)
    row$transition <- transition
    rows[[length(rows) + 1]] <- row
  }
  do.call(rbind, rows)
}

meta_by_sensitivity <- function(effects) {
  rows <- list()
  for (sensitivity in unique(effects$sensitivity)) {
    for (transition in unique(effects$transition)) {
      z <- effects[effects$sensitivity == sensitivity & effects$transition == transition, ]
      if (nrow(z) < 2) next
      row <- meta_estimates(z)
      row$sensitivity <- sensitivity
      row$transition <- transition
      rows[[length(rows) + 1]] <- row
    }
  }
  do.call(rbind, rows)
}

bivariate_fixed_joint <- function(joints) {
  precision_sum <- matrix(0, 2, 2)
  score_sum <- matrix(0, 2, 1)
  for (joint in joints) {
    effects <- joint$effects
    b <- matrix(effects$log_tir[match(c("Progression", "Improvement"), effects$transition)], 2, 1)
    # Recover the 2x2 covariance for transition-specific exposure terms.
    cn <- names(joint$coefficients)
    pterm <- find_transition_term(cn, "Progression", "excess_bmi25_z")
    iterm <- find_transition_term(cn, "Improvement", "excess_bmi25_z")
    v <- joint$covariance[c(pterm, iterm), c(pterm, iterm), drop = FALSE]
    precision <- MASS::ginv(v)
    precision_sum <- precision_sum + precision
    score_sum <- score_sum + precision %*% b
  }
  covariance <- MASS::ginv(precision_sum)
  beta <- covariance %*% score_sum
  statistic <- as.numeric(t(beta) %*% MASS::ginv(covariance) %*% beta)
  data.frame(
    progression_tir = exp(beta[1]), improvement_tir = exp(beta[2]),
    joint_wald_chisq = statistic, joint_wald_df = 2,
    joint_wald_p = pchisq(statistic, 2, lower.tail = FALSE),
    method = "Bivariate fixed-effect GLS using cohort-specific within-study covariance"
  )
}

datasets <- lapply(c("HRS", "CHNS", "ELSA"), prepare_data)
names(datasets) <- c("HRS", "CHNS", "ELSA")

main_joints <- lapply(names(datasets), function(cohort) fit_joint(datasets[[cohort]], cohort))
names(main_joints) <- names(datasets)
main_effects <- do.call(rbind, lapply(main_joints, `[[`, "effects"))
joint_tests <- do.call(rbind, lapply(main_joints, `[[`, "tests"))
meta_results <- meta_by_transition(main_effects)
bivariate_joint <- bivariate_fixed_joint(main_joints)

expanded_joints <- lapply(names(datasets), function(cohort) fit_joint(datasets[[cohort]], cohort, covar_set = "expanded"))
bp_joints <- lapply(names(datasets), function(cohort) fit_joint(datasets[[cohort]], cohort, bp_adjust = TRUE))
person_joints <- lapply(names(datasets), function(cohort) fit_joint(datasets[[cohort]], cohort, cluster_mode = "person"))

model_comparisons <- rbind(
  transform(main_effects, sensitivity = "Minimal prespecified adjustment"),
  transform(do.call(rbind, lapply(expanded_joints, `[[`, "effects")), sensitivity = "Extended health-status adjustment"),
  transform(do.call(rbind, lapply(bp_joints, `[[`, "effects")), sensitivity = "Additional baseline SBP and DBP adjustment"),
  transform(do.call(rbind, lapply(person_joints, `[[`, "effects")), sensitivity = "Person-clustered standard errors")
)
model_comparison_meta <- meta_by_sensitivity(model_comparisons)

chns_threshold_rows <- list()
for (exposure in c("excess_bmi25_z", "excess_bmi24_z", "excess_bmi23_z")) {
  joint <- fit_joint(datasets$CHNS, "CHNS", exposure = exposure)
  row <- joint$effects
  row$threshold_kg_m2 <- sub("excess_bmi", "", sub("_z", "", exposure))
  chns_threshold_rows[[length(chns_threshold_rows) + 1]] <- row
}
chns_thresholds <- do.call(rbind, chns_threshold_rows)

absolute <- lapply(names(main_joints), function(cohort) marginal_percentile_probabilities(main_joints[[cohort]], cohort))
absolute_probabilities <- do.call(rbind, lapply(absolute, `[[`, "probabilities"))
absolute_differences <- do.call(rbind, lapply(absolute, `[[`, "differences"))

continuous_changes <- do.call(rbind, lapply(names(datasets), function(cohort) {
  rbind(fit_continuous_change(datasets[[cohort]], cohort, "delta_sbp"),
        fit_continuous_change(datasets[[cohort]], cohort, "delta_dbp"))
}))

weight_sensitivities <- do.call(rbind, lapply(names(datasets), function(cohort) {
  do.call(rbind, lapply(names(endpoint_labels), function(endpoint) weight_sensitivity(datasets[[cohort]], cohort, endpoint)))
}))
weight_distributions <- do.call(rbind, lapply(names(datasets), function(cohort) weight_distribution(datasets[[cohort]], cohort)))

additional_robustness <- do.call(rbind, lapply(names(datasets), function(cohort) {
  d <- datasets[[cohort]]
  cvd_var <- ifelse(cohort == "HRS", "cvd_history", "cvd_history")
  filters <- list(
    "Exclude baseline cardiovascular disease" = !is.na(d[[cvd_var]]) & d[[cvd_var]] == 0,
    "Exclude BMI loss greater than 5 percent" = !is.na(d$bmi_change_pct) & d$bmi_change_pct >= -5,
    "Exclude BMI loss greater than 10 percent" = !is.na(d$bmi_change_pct) & d$bmi_change_pct >= -10,
    "Exclude mean BMI below 18.5 kg/m2" = !is.na(d$bmi_mean) & d$bmi_mean >= 18.5
  )
  if (cohort %in% c("HRS", "CHNS")) {
    first_end <- min(d$end_year, na.rm = TRUE)
    filters[["Exclude first outcome interval"]] <- !is.na(d$end_year) & d$end_year > first_end
  }
  rows <- list()
  for (label in names(filters)) {
    joint <- fit_joint(d[filters[[label]], , drop = FALSE], cohort)
    row <- joint$effects
    row$sensitivity <- label
    rows[[length(rows) + 1]] <- row
  }
  do.call(rbind, rows)
}))
additional_robustness_meta <- meta_by_sensitivity(additional_robustness)

cluster_summary <- do.call(rbind, lapply(names(datasets), function(cohort) {
  d <- datasets[[cohort]]
  data.frame(
    cohort = cohort, participants = length(unique(d$person_id)), intervals = nrow(d),
    sampling_cluster_definition = ifelse(cohort == "HRS", "stratum-by-SECU PSU",
                                         ifelse(cohort == "CHNS", "community", "baseline household")),
    sampling_clusters = length(unique(d$psu_id)), strata = ifelse(cohort == "HRS", length(unique(d$strata_value)), NA),
    maximum_intervals_per_person = max(table(d$person_id))
  )
}))

write.csv(main_effects, file.path(out, "primary_minimal_joint_effects.csv"), row.names = FALSE)
write.csv(joint_tests, file.path(out, "primary_joint_wald_and_directional_contrasts.csv"), row.names = FALSE)
write.csv(meta_results, file.path(out, "primary_meta_results.csv"), row.names = FALSE)
write.csv(bivariate_joint, file.path(out, "bivariate_fixed_joint_test.csv"), row.names = FALSE)
write.csv(model_comparisons, file.path(out, "minimal_expanded_bp_cluster_sensitivity.csv"), row.names = FALSE)
write.csv(model_comparison_meta, file.path(out, "minimal_expanded_bp_cluster_sensitivity_meta.csv"), row.names = FALSE)
write.csv(chns_thresholds, file.path(out, "chns_bmi_threshold_sensitivity.csv"), row.names = FALSE)
write.csv(absolute_probabilities, file.path(out, "percentile_standardized_probabilities.csv"), row.names = FALSE)
write.csv(absolute_differences, file.path(out, "percentile_standardized_risk_differences.csv"), row.names = FALSE)
write.csv(continuous_changes, file.path(out, "continuous_sbp_dbp_change_models.csv"), row.names = FALSE)
write.csv(weight_sensitivities, file.path(out, "weighting_sensitivity_models.csv"), row.names = FALSE)
write.csv(weight_distributions, file.path(out, "combined_weight_distribution_and_ess.csv"), row.names = FALSE)
write.csv(cluster_summary, file.path(out, "clustering_units_summary.csv"), row.names = FALSE)
write.csv(additional_robustness, file.path(out, "disease_weightloss_early_interval_sensitivity.csv"), row.names = FALSE)
write.csv(additional_robustness_meta, file.path(out, "disease_weightloss_early_interval_sensitivity_meta.csv"), row.names = FALSE)

cat("Primary minimal-adjustment effects\n")
print(main_effects[, c("cohort", "transition", "tir", "ci_low", "ci_high", "p_value", "events")])
cat("\nJoint Wald and directional contrasts\n")
print(joint_tests[, c("cohort", "joint_wald_chisq", "joint_wald_p", "directional_ratio", "directional_ci_low", "directional_ci_high", "directional_p")])
cat("\nMeta-analysis\n")
print(meta_results[, c("transition", "reml_tir", "reml_ci_low", "reml_ci_high", "modified_hk_ci_low", "modified_hk_ci_high", "modified_hk_p", "i2_percent")])
cat("\nBivariate fixed-effect joint test\n")
print(bivariate_joint)
cat("\nCHNS BMI threshold sensitivity\n")
print(chns_thresholds[, c("threshold_kg_m2", "transition", "tir", "ci_low", "ci_high", "p_value")])
