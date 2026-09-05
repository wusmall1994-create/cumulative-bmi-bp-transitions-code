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
hrs_out <- file.path(output_root, "hrs_analysis")
chns_out <- file.path(output_root, "chns_analysis")
out <- file.path(output_root, "meta_hrs_chns")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

endpoint_labels <- c(
  hypertension_onset = "Hypertension onset",
  treatment_initiation = "Treatment initiation",
  untreated_bp_improvement = "Untreated BP improvement",
  control_loss = "Loss of BP control",
  control_achievement = "Achievement of BP control"
)
endpoints <- names(endpoint_labels)

hrs <- read.csv(file.path(hrs_out, "analysis_intervals.csv"), fileEncoding = "UTF-8", check.names = FALSE)
hrs$interval <- factor(hrs$interval, levels = c("2014-2018", "2018-2022"))
hrs$race_f <- factor(hrs$raracem)
hrs$region4 <- ifelse(hrs$region14 %in% c(1, 2), 1,
                      ifelse(hrs$region14 %in% c(3, 4), 2,
                             ifelse(hrs$region14 %in% c(5, 6, 7), 3,
                                    ifelse(hrs$region14 %in% c(8, 9), 4, NA))))
hrs$region4_f <- factor(hrs$region4)
hrs$interval_years <- hrs$end_year - hrs$start_year

chns <- read.csv(file.path(chns_out, "chns_intervals.csv"), fileEncoding = "UTF-8", check.names = FALSE)
chns$interval <- factor(chns$interval, levels = c("2006-2009", "2009-2011"))
chns$province_f <- factor(chns$province)
chns$interval_years <- chns$end_year - chns$start_year

hrs_covars <- c(
  "age_start10", "female", "race_f", "raedyrs", "married_partnered", "current_smoker",
  "current_drinker", "r12diabe", "cvd_history", "r12shlt", "region4_f", "rural", "interval"
)
chns_full_covars <- c(
  "age_start10", "female", "educ_years", "current_smoker", "current_drinker", "diabetes",
  "cvd_history", "self_health", "province_f", "urban", "interval"
)
chns_covars_for <- function(endpoint) {
  if (endpoint == "control_loss") return(c("age_start10", "female", "urban", "interval"))
  if (endpoint == "control_achievement") return(c("age_start10", "female", "educ_years", "diabetes", "urban", "interval"))
  chns_full_covars
}

fit_tir <- function(data, cohort, endpoint, exposure_terms, keep_term, model_label) {
  covars <- if (cohort == "HRS") hrs_covars else chns_covars_for(endpoint)
  cluster_var <- if (cohort == "HRS") "secu" else "commid"
  needed <- unique(c(endpoint, "model_weight", cluster_var, "interval_years", covars, exposure_terms))
  if (cohort == "HRS") needed <- unique(c(needed, "stratum"))
  keep <- !is.na(data[[endpoint]]) & !is.na(data$model_weight) & data$model_weight > 0 & complete.cases(data[, needed, drop = FALSE])
  z <- data[keep, , drop = FALSE]
  if (nrow(z) < 60 || sum(z[[endpoint]]) < 20 || length(unique(z[[endpoint]])) < 2) return(NULL)
  if (cohort == "HRS") {
    des <- svydesign(ids = ~secu, strata = ~stratum, weights = ~model_weight, data = z, nest = TRUE)
  } else {
    des <- svydesign(ids = ~commid, weights = ~model_weight, data = z)
  }
  f <- as.formula(paste(endpoint, "~", paste(c(exposure_terms, covars, "offset(log(interval_years))"), collapse = " + ")))
  fit <- tryCatch(svyglm(f, design = des, family = quasibinomial(link = "cloglog")), error = function(e) NULL)
  if (is.null(fit) || !(keep_term %in% names(coef(fit)))) return(NULL)
  b <- coef(fit)[keep_term]
  se <- sqrt(vcov(fit)[keep_term, keep_term])
  p <- 2 * pnorm(abs(b / se), lower.tail = FALSE)
  data.frame(
    cohort = cohort, endpoint = endpoint, endpoint_label = endpoint_labels[[endpoint]], model = model_label,
    term = keep_term, log_tir = b, se = se, tir = exp(b), ci_low = exp(b - 1.96 * se), ci_high = exp(b + 1.96 * se),
    p_value = p, n_intervals = nrow(z), events = sum(z[[endpoint]]), participants = length(unique(z$person_id)),
    interval_years = ifelse(cohort == "HRS", "4", "2-3"), stringsAsFactors = FALSE
  )
}

cohort_results <- list()
for (ep in endpoints) {
  cohort_results[[length(cohort_results) + 1]] <- fit_tir(hrs, "HRS", ep, "excess_bmi25_z", "excess_bmi25_z", "Cumulative excess BMI burden")
  cohort_results[[length(cohort_results) + 1]] <- fit_tir(hrs, "HRS", ep, c("bmi_vim_z", "bmi_mean_z"), "bmi_vim_z", "BMI variability (VIM)")
  cohort_results[[length(cohort_results) + 1]] <- fit_tir(chns, "CHNS", ep, "excess_bmi25_z", "excess_bmi25_z", "Cumulative excess BMI burden")
  cohort_results[[length(cohort_results) + 1]] <- fit_tir(chns, "CHNS", ep, c("bmi_vim_z", "bmi_mean_z"), "bmi_vim_z", "BMI variability (VIM)")
}
cohort_results <- do.call(rbind, Filter(Negate(is.null), cohort_results))
write.csv(cohort_results, file.path(out, "harmonized_cohort_tir.csv"), row.names = FALSE, fileEncoding = "UTF-8")

meta_dl <- function(z) {
  k <- nrow(z)
  w <- 1 / z$se^2
  fixed <- sum(w * z$log_tir) / sum(w)
  q <- sum(w * (z$log_tir - fixed)^2)
  c_term <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (q - (k - 1)) / c_term)
  wr <- 1 / (z$se^2 + tau2)
  pooled <- sum(wr * z$log_tir) / sum(wr)
  pooled_se <- sqrt(1 / sum(wr))
  p <- 2 * pnorm(abs(pooled / pooled_se), lower.tail = FALSE)
  q_p <- pchisq(q, df = k - 1, lower.tail = FALSE)
  i2 <- ifelse(q > 0, max(0, (q - (k - 1)) / q) * 100, 0)
  data.frame(
    k = k, log_tir = pooled, se = pooled_se, tir = exp(pooled),
    ci_low = exp(pooled - 1.96 * pooled_se), ci_high = exp(pooled + 1.96 * pooled_se),
    p_value = p, tau2 = tau2, q = q, q_p = q_p, i2_percent = i2,
    fixed_tir = exp(fixed), fixed_ci_low = exp(fixed - 1.96 / sqrt(sum(w))), fixed_ci_high = exp(fixed + 1.96 / sqrt(sum(w)))
  )
}

meta_rows <- list()
for (ep in endpoints) {
  for (model in c("Cumulative excess BMI burden", "BMI variability (VIM)")) {
    z <- cohort_results[cohort_results$endpoint == ep & cohort_results$model == model, , drop = FALSE]
    if (nrow(z) != 2) next
    m <- meta_dl(z)
    m$endpoint <- ep
    m$endpoint_label <- endpoint_labels[[ep]]
    m$model <- model
    meta_rows[[length(meta_rows) + 1]] <- m
  }
}
meta <- do.call(rbind, meta_rows)
meta$p_fdr <- p.adjust(meta$p_value, method = "BH")
meta <- meta[, c("endpoint", "endpoint_label", "model", "k", "tir", "ci_low", "ci_high", "p_value", "p_fdr",
                 "tau2", "q", "q_p", "i2_percent", "fixed_tir", "fixed_ci_low", "fixed_ci_high", "log_tir", "se")]
write.csv(meta, file.path(out, "random_effects_meta_results.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Cohort-by-exposure heterogeneity is represented by Cochran Q (df=1) in the meta table.
decision <- data.frame(
  item = c("Harmonized effect", "Interval handling", "Pooling", "Interpretation boundary"),
  decision = c(
    "Transition intensity ratio from complementary log-log panel models",
    "Offset log(interval years): HRS 4 years; CHNS 3 and 2 years",
    "DerSimonian-Laird random effects with fixed-effect estimates also reported; k=2",
    "Exact transition times and within-interval intermediate states are unobserved"
  )
)
write.csv(decision, file.path(out, "meta_statistical_decisions.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Cohort plus pooled forest figure.
plot_cohort <- cohort_results[, c("endpoint", "endpoint_label", "model", "cohort", "tir", "ci_low", "ci_high")]
plot_meta <- data.frame(endpoint = meta$endpoint, endpoint_label = meta$endpoint_label, model = meta$model,
                        cohort = "Pooled", tir = meta$tir, ci_low = meta$ci_low, ci_high = meta$ci_high)
plot_data <- rbind(plot_cohort, plot_meta)
plot_data$model <- factor(plot_data$model, levels = c("Cumulative excess BMI burden", "BMI variability (VIM)"),
                          labels = c("Excess BMI burden", "BMI variability (VIM)"))
plot_data$cohort <- factor(plot_data$cohort, levels = c("HRS", "CHNS", "Pooled"))
plot_data$endpoint_label <- factor(plot_data$endpoint_label, levels = rev(unname(endpoint_labels)))

theme_set(theme_classic(base_size = 7, base_family = "sans") +
            theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
                  plot.title = element_text(size = 8, face = "bold"), strip.text = element_text(size = 7, face = "bold"),
                  legend.title = element_text(size = 6.5), legend.text = element_text(size = 6), panel.grid = element_blank()))

p <- ggplot(plot_data, aes(tir, endpoint_label, colour = cohort, shape = cohort)) +
  geom_vline(xintercept = 1, colour = "#777777", linewidth = 0.35, linetype = 2) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.15,
                position = position_dodge(width = 0.55), linewidth = 0.45) +
  geom_point(position = position_dodge(width = 0.55), size = 1.8) +
  facet_wrap(~model, nrow = 1, scales = "free_x") +
  scale_x_log10() +
  scale_colour_manual(values = c(HRS = "#34688A", CHNS = "#B0684B", Pooled = "#222222")) +
  labs(x = "Transition intensity ratio per 1 SD (log scale)", y = NULL, colour = NULL, shape = NULL,
       title = "Harmonized cohort-specific and pooled associations") +
  theme(legend.position = "top")

save_pub <- function(plot, stem, width_mm = 183, height_mm = 110, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(file.path(out, paste0(stem, ".svg")), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(file.path(out, paste0(stem, ".pdf")), width = w, height = h, family = "sans"); print(plot); dev.off()
  ragg::agg_tiff(file.path(out, paste0(stem, ".tiff")), width = w, height = h, units = "in", res = dpi, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(file.path(out, paste0(stem, "_preview.png")), width = w, height = h, units = "in", res = 180); print(plot); dev.off()
}
save_pub(p, "figure3_cross_cohort_meta")

cat("Harmonized cohort results:\n")
print(cohort_results[, c("cohort", "endpoint", "model", "tir", "ci_low", "ci_high", "p_value")])
cat("\nRandom-effects meta-analysis:\n")
print(meta[, c("endpoint", "model", "tir", "ci_low", "ci_high", "p_value", "p_fdr", "i2_percent")])
