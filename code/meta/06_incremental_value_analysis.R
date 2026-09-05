options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

suppressPackageStartupMessages({
  library(survey)
})

set.seed(20260831)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
output_root <- Sys.getenv("BMI_BP_OUTPUT_DIR", file.path(root, "outputs"))
out <- file.path(output_root, "incremental_value")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

endpoints <- c("hypertension_onset", "untreated_bp_improvement")
endpoint_labels <- c(
  hypertension_onset = "Hypertension onset",
  untreated_bp_improvement = "Untreated BP improvement"
)

standardize <- function(x) as.numeric((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

# HRS
hrs_root <- file.path(output_root, "hrs_analysis")
hrs <- read.csv(file.path(hrs_root, "analysis_intervals.csv"), fileEncoding = "UTF-8", check.names = FALSE)
hrs_base <- read.csv(file.path(hrs_root, "analysis_baseline.csv"), fileEncoding = "UTF-8", check.names = FALSE)
hrs <- merge(hrs, hrs_base[, c("person_id", "r8pmbmi", "r12pmbmi")], by = "person_id", all.x = TRUE)
hrs$cohort <- "HRS"
hrs$baseline_bmi <- hrs$r8pmbmi
hrs$last_bmi <- hrs$r12pmbmi
hrs$interval_years <- hrs$end_year - hrs$start_year
hrs$interval <- factor(hrs$interval, levels = c("2014-2018", "2018-2022"))
hrs$race_f <- factor(hrs$raracem)
hrs$region4 <- ifelse(hrs$region14 %in% c(1, 2), 1,
                      ifelse(hrs$region14 %in% c(3, 4), 2,
                             ifelse(hrs$region14 %in% c(5, 6, 7), 3,
                                    ifelse(hrs$region14 %in% c(8, 9), 4, NA))))
hrs$region_f <- factor(hrs$region4)
hrs$cluster <- hrs$secu
hrs$strata_value <- hrs$stratum

# CHNS
chns_root <- file.path(output_root, "chns_analysis")
chns <- read.csv(file.path(chns_root, "chns_intervals.csv"), fileEncoding = "UTF-8", check.names = FALSE)
chns_base <- read.csv(file.path(chns_root, "chns_baseline_2006.csv"), fileEncoding = "UTF-8", check.names = FALSE)
chns <- merge(chns, chns_base[, c("person_id", "bmi_2000", "bmi_2006")], by = "person_id", all.x = TRUE)
chns$cohort <- "CHNS"
chns$baseline_bmi <- chns$bmi_2000
chns$last_bmi <- chns$bmi_2006
# The CHNS interval file stores the trapezoidal BMI area but not its
# time-standardized counterpart.  The exposure window is fixed at six years
# (2000-2006), so the time-weighted mean is exactly AUC / 6.
chns$bmi_twmean <- chns$bmi_auc / 6
chns$interval_years <- chns$end_year - chns$start_year
chns$interval <- factor(chns$interval, levels = c("2006-2009", "2009-2011"))
chns$province_f <- factor(chns$province)
chns$cluster <- chns$commid
chns$strata_value <- NA

# ELSA
elsa_root <- file.path(output_root, "elsa_analysis")
elsa <- read.csv(file.path(elsa_root, "elsa_intervals.csv"), fileEncoding = "UTF-8", check.names = FALSE)
elsa$cohort <- "ELSA"
elsa$baseline_bmi <- elsa$bmi_w2
elsa$last_bmi <- elsa$bmi_w6
elsa$interval <- factor(elsa$interval, levels = c("W6-W8", "W6-W9"))
elsa$education_f <- factor(elsa$education)
elsa$region_f <- factor(elsa$region)
elsa$cluster <- elsa$idahhw6
elsa$strata_value <- NA

datasets <- list(HRS = hrs, CHNS = chns, ELSA = elsa)
for (cohort in names(datasets)) {
  data <- datasets[[cohort]]
  for (variable in c("baseline_bmi", "last_bmi", "bmi_mean", "bmi_twmean", "excess_bmi25", "bmi_vim")) {
    data[[paste0(variable, "_cmp_z")]] <- standardize(data[[variable]])
  }
  datasets[[cohort]] <- data
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

model_specs <- list(
  covariates_only = character(0),
  baseline_bmi = "baseline_bmi_cmp_z",
  last_bmi = "last_bmi_cmp_z",
  mean_bmi = "bmi_mean_cmp_z",
  time_weighted_mean_bmi = "bmi_twmean_cmp_z",
  cumulative_excess_burden = "excess_bmi25_cmp_z",
  burden_beyond_baseline = c("baseline_bmi_cmp_z", "excess_bmi25_cmp_z"),
  burden_beyond_mean = c("bmi_mean_cmp_z", "excess_bmi25_cmp_z"),
  vim_beyond_mean = c("bmi_mean_cmp_z", "bmi_vim_cmp_z")
)

model_labels <- c(
  covariates_only = "Covariates only",
  baseline_bmi = "Baseline BMI",
  last_bmi = "Last BMI",
  mean_bmi = "Mean BMI",
  time_weighted_mean_bmi = "Time-weighted mean BMI",
  cumulative_excess_burden = "Cumulative excess BMI burden",
  burden_beyond_baseline = "Burden beyond baseline BMI",
  burden_beyond_mean = "Burden beyond mean BMI",
  vim_beyond_mean = "VIM beyond mean BMI"
)

weighted_auc <- function(y, score, weight) {
  data <- data.frame(y = y, score = score, weight = weight)
  data <- data[complete.cases(data) & data$weight > 0, ]
  positive_total <- sum(data$weight[data$y == 1])
  negative_total <- sum(data$weight[data$y == 0])
  if (positive_total <= 0 || negative_total <= 0) return(NA_real_)
  groups <- aggregate(weight ~ score + y, data = data, sum)
  scores <- sort(unique(groups$score))
  cumulative_negative <- 0
  concordant <- 0
  for (score_value in scores) {
    at_score <- groups[groups$score == score_value, ]
    positive_weight <- sum(at_score$weight[at_score$y == 1])
    negative_weight <- sum(at_score$weight[at_score$y == 0])
    concordant <- concordant + positive_weight * (cumulative_negative + 0.5 * negative_weight)
    cumulative_negative <- cumulative_negative + negative_weight
  }
  concordant / (positive_total * negative_total)
}

fit_survey_effects <- function(data, cohort, endpoint, predictors, model_id) {
  covars <- covars_for(cohort)
  needed <- unique(c(endpoint, "model_weight", "cluster", "interval_years", covars, predictors))
  if (cohort == "HRS") needed <- unique(c(needed, "strata_value"))
  keep <- !is.na(data[[endpoint]]) & !is.na(data$model_weight) & data$model_weight > 0 & complete.cases(data[, needed, drop = FALSE])
  z <- data[keep, ]
  if (nrow(z) < 60 || sum(z[[endpoint]]) < 20 || length(unique(z[[endpoint]])) < 2) return(NULL)
  design <- if (cohort == "HRS") {
    svydesign(ids = ~cluster, strata = ~strata_value, weights = ~model_weight, data = z, nest = TRUE)
  } else {
    svydesign(ids = ~cluster, weights = ~model_weight, data = z)
  }
  rhs <- c(predictors, covars, "offset(log(interval_years))")
  formula <- as.formula(paste(endpoint, "~", paste(rhs, collapse = " + ")))
  fit <- tryCatch(svyglm(formula, design = design, family = quasibinomial(link = "cloglog")), error = function(e) NULL)
  if (is.null(fit) || length(predictors) == 0) return(NULL)
  rows <- list()
  for (term in predictors) {
    if (!(term %in% names(coef(fit)))) next
    beta <- coef(fit)[term]
    se <- sqrt(vcov(fit)[term, term])
    rows[[length(rows) + 1]] <- data.frame(
      cohort = cohort, endpoint = endpoint, endpoint_label = endpoint_labels[[endpoint]], model_id = model_id,
      model = model_labels[[model_id]], term = term, log_tir = beta, se = se, tir = exp(beta),
      ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se),
      p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE), n = nrow(z), events = sum(z[[endpoint]])
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

cross_validated_metrics <- function(data, cohort, endpoint, predictors, model_id, folds = 5) {
  covars <- covars_for(cohort)
  all_exposures <- unique(unlist(model_specs))
  needed <- unique(c(endpoint, "model_weight", "person_id", "interval_years", covars, all_exposures))
  keep <- !is.na(data[[endpoint]]) & !is.na(data$model_weight) & data$model_weight > 0 & complete.cases(data[, needed, drop = FALSE])
  z <- data[keep, ]
  if (nrow(z) < 100 || sum(z[[endpoint]]) < 20) return(NULL)
  people <- unique(z$person_id)
  # Reuse an identical person-level fold allocation for every candidate model
  # within a cohort-endpoint comparison.  Otherwise delta-AUC/Brier estimates
  # would partly reflect different random splits rather than model content.
  cohort_seed <- match(cohort, c("HRS", "CHNS", "ELSA")) * 100L
  endpoint_seed <- match(endpoint, endpoints) * 10L
  set.seed(20260831L + cohort_seed + endpoint_seed)
  fold_map <- setNames(sample(rep(seq_len(folds), length.out = length(people))), people)
  z$fold <- unname(fold_map[as.character(z$person_id)])
  prediction <- rep(NA_real_, nrow(z))
  rhs <- c(predictors, covars, "offset(log(interval_years))")
  formula <- as.formula(paste(endpoint, "~", paste(rhs, collapse = " + ")))
  for (fold in seq_len(folds)) {
    train <- z$fold != fold
    test <- z$fold == fold
    training_weight <- z$model_weight[train] / mean(z$model_weight[train])
    fit <- suppressWarnings(tryCatch(glm(formula, data = z[train, ], family = binomial(link = "cloglog"),
                                         weights = training_weight), error = function(e) NULL))
    if (is.null(fit)) next
    prediction[test] <- tryCatch(as.numeric(predict(fit, newdata = z[test, ], type = "response")), error = function(e) rep(NA_real_, sum(test)))
  }
  valid <- is.finite(prediction)
  y <- z[[endpoint]][valid]
  weight <- z$model_weight[valid]
  prediction <- pmin(pmax(prediction[valid], 1e-6), 1 - 1e-6)
  brier <- sum(weight * (y - prediction)^2) / sum(weight)
  auc <- weighted_auc(y, prediction, weight)
  data.frame(cohort = cohort, endpoint = endpoint, model_id = model_id, model = model_labels[[model_id]],
             n = sum(valid), events = sum(y), cv_brier = brier, cv_auc = auc)
}

effect_rows <- list()
metric_rows <- list()
for (cohort in names(datasets)) {
  data <- datasets[[cohort]]
  for (endpoint in endpoints) {
    for (model_id in names(model_specs)) {
      predictors <- model_specs[[model_id]]
      if (length(predictors) > 0) {
        effect_rows[[length(effect_rows) + 1]] <- fit_survey_effects(data, cohort, endpoint, predictors, model_id)
      }
      metric_rows[[length(metric_rows) + 1]] <- cross_validated_metrics(data, cohort, endpoint, predictors, model_id)
    }
  }
}
effects <- do.call(rbind, Filter(Negate(is.null), effect_rows))
metrics <- do.call(rbind, Filter(Negate(is.null), metric_rows))
write.csv(effects, file.path(out, "cohort_exposure_comparison_effects.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(metrics, file.path(out, "cohort_cross_validated_metrics.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Incremental prediction metrics relative to baseline-BMI and mean-BMI models.
metric_delta_rows <- list()
for (cohort in names(datasets)) {
  for (endpoint in endpoints) {
    z <- metrics[metrics$cohort == cohort & metrics$endpoint == endpoint, ]
    for (comparison in list(c("baseline_bmi", "burden_beyond_baseline"), c("mean_bmi", "burden_beyond_mean"))) {
      base <- z[z$model_id == comparison[1], ]
      expanded <- z[z$model_id == comparison[2], ]
      if (nrow(base) == 1 && nrow(expanded) == 1) {
        metric_delta_rows[[length(metric_delta_rows) + 1]] <- data.frame(
          cohort = cohort, endpoint = endpoint, base_model = model_labels[[comparison[1]]],
          expanded_model = model_labels[[comparison[2]]], delta_auc = expanded$cv_auc - base$cv_auc,
          delta_brier = expanded$cv_brier - base$cv_brier
        )
      }
    }
  }
}
metric_deltas <- do.call(rbind, metric_delta_rows)
write.csv(metric_deltas, file.path(out, "incremental_prediction_metrics.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Exposure correlations, two-predictor VIFs and condition numbers.
collinearity_rows <- list()
for (cohort in names(datasets)) {
  data <- datasets[[cohort]]
  pairs <- list(
    c("baseline_bmi_cmp_z", "excess_bmi25_cmp_z"),
    c("bmi_mean_cmp_z", "excess_bmi25_cmp_z"),
    c("bmi_mean_cmp_z", "bmi_vim_cmp_z"),
    c("bmi_twmean_cmp_z", "bmi_mean_cmp_z")
  )
  for (pair in pairs) {
    complete <- complete.cases(data[, pair])
    x <- as.matrix(data[complete, pair])
    correlation <- cor(x[, 1], x[, 2])
    singular <- svd(scale(x, center = TRUE, scale = FALSE))$d
    condition <- max(singular) / min(singular)
    collinearity_rows[[length(collinearity_rows) + 1]] <- data.frame(
      cohort = cohort, variable_1 = pair[1], variable_2 = pair[2], n = sum(complete), correlation = correlation,
      vif_two_predictor = 1 / (1 - correlation^2), condition_number = condition
    )
  }
}
collinearity <- do.call(rbind, collinearity_rows)
write.csv(collinearity, file.path(out, "exposure_collinearity_diagnostics.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# DerSimonian-Laird pooling for each exposure parameterization and joint burden term.
meta_dl <- function(z) {
  weights <- 1 / z$se^2
  fixed <- sum(weights * z$log_tir) / sum(weights)
  q <- sum(weights * (z$log_tir - fixed)^2)
  c_term <- sum(weights) - sum(weights^2) / sum(weights)
  tau2 <- max(0, (q - (nrow(z) - 1)) / c_term)
  random_weights <- 1 / (z$se^2 + tau2)
  pooled <- sum(random_weights * z$log_tir) / sum(random_weights)
  pooled_se <- sqrt(1 / sum(random_weights))
  data.frame(tir = exp(pooled), ci_low = exp(pooled - 1.96 * pooled_se), ci_high = exp(pooled + 1.96 * pooled_se),
             p_value = 2 * pnorm(abs(pooled / pooled_se), lower.tail = FALSE), tau2 = tau2,
             i2_percent = ifelse(q > 0, max(0, (q - (nrow(z) - 1)) / q) * 100, 0))
}

meta_targets <- list(
  baseline_bmi = "baseline_bmi_cmp_z",
  last_bmi = "last_bmi_cmp_z",
  mean_bmi = "bmi_mean_cmp_z",
  time_weighted_mean_bmi = "bmi_twmean_cmp_z",
  cumulative_excess_burden = "excess_bmi25_cmp_z",
  burden_beyond_baseline = "excess_bmi25_cmp_z",
  burden_beyond_mean = "excess_bmi25_cmp_z",
  vim_beyond_mean = "bmi_vim_cmp_z"
)
meta_rows <- list()
for (endpoint in endpoints) {
  for (model_id in names(meta_targets)) {
    z <- effects[effects$endpoint == endpoint & effects$model_id == model_id & effects$term == meta_targets[[model_id]], ]
    if (nrow(z) != 3) next
    pooled <- meta_dl(z)
    pooled$endpoint <- endpoint
    pooled$endpoint_label <- endpoint_labels[[endpoint]]
    pooled$model_id <- model_id
    pooled$model <- model_labels[[model_id]]
    pooled$k <- nrow(z)
    meta_rows[[length(meta_rows) + 1]] <- pooled
  }
}
meta <- do.call(rbind, meta_rows)
meta$p_fdr <- p.adjust(meta$p_value, method = "BH")
write.csv(meta, file.path(out, "incremental_exposure_meta_results.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Leave-one-cohort-out diagnostics are especially important for the joint
# burden models: a pooled coefficient can otherwise look precise while being
# dominated by the only cohort with modest collinearity.
loo_rows <- list()
for (endpoint in endpoints) {
  for (model_id in names(meta_targets)) {
    z_all <- effects[effects$endpoint == endpoint & effects$model_id == model_id & effects$term == meta_targets[[model_id]], ]
    if (nrow(z_all) != 3) next
    for (omitted in z_all$cohort) {
      z <- z_all[z_all$cohort != omitted, ]
      pooled <- meta_dl(z)
      pooled$endpoint <- endpoint
      pooled$model_id <- model_id
      pooled$model <- model_labels[[model_id]]
      pooled$omitted_cohort <- omitted
      pooled$k <- nrow(z)
      loo_rows[[length(loo_rows) + 1]] <- pooled
    }
  }
}
loo <- do.call(rbind, loo_rows)
write.csv(loo, file.path(out, "leave_one_cohort_out_meta.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("Collinearity diagnostics:\n")
print(collinearity)
cat("\nIncremental prediction metrics:\n")
print(metric_deltas)
cat("\nIncremental pooled effects:\n")
print(meta[, c("endpoint", "model", "tir", "ci_low", "ci_high", "p_value", "p_fdr", "i2_percent")])
