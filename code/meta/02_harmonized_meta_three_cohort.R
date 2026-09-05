options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
  library(ragg)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
root <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)
output_root <- Sys.getenv("BMI_BP_OUTPUT_DIR", file.path(root, "outputs"))
two_out <- file.path(output_root, "meta_hrs_chns")
elsa_out <- file.path(output_root, "elsa_analysis")
out <- file.path(output_root, "meta_three_cohort")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

endpoint_labels <- c(
  hypertension_onset = "Hypertension onset",
  treatment_initiation = "Treatment initiation",
  untreated_bp_improvement = "Untreated BP improvement",
  control_loss = "Loss of BP control",
  control_achievement = "Achievement of BP control"
)
endpoints <- names(endpoint_labels)
models <- c("Cumulative excess BMI burden", "BMI variability (VIM)")

prior <- read.csv(file.path(two_out, "harmonized_cohort_tir.csv"), fileEncoding = "UTF-8", check.names = FALSE)
elsa <- read.csv(file.path(elsa_out, "elsa_harmonized_tir.csv"), fileEncoding = "UTF-8", check.names = FALSE)
keep <- c("cohort", "endpoint", "endpoint_label", "model", "term", "log_tir", "se", "tir", "ci_low", "ci_high",
          "p_value", "n_intervals", "events", "participants", "interval_years")
cohort_results <- rbind(prior[, keep], elsa[, keep])
cohort_results <- cohort_results[cohort_results$endpoint %in% endpoints & cohort_results$model %in% models, ]
cohort_results <- cohort_results[order(match(cohort_results$endpoint, endpoints), match(cohort_results$model, models),
                                       match(cohort_results$cohort, c("HRS", "CHNS", "ELSA"))), ]
cohort_results$p_fdr <- p.adjust(cohort_results$p_value, method = "BH")
write.csv(cohort_results, file.path(out, "three_cohort_harmonized_tir.csv"), row.names = FALSE, fileEncoding = "UTF-8")

meta_dl <- function(z) {
  k <- nrow(z)
  weights <- 1 / z$se^2
  fixed <- sum(weights * z$log_tir) / sum(weights)
  q <- sum(weights * (z$log_tir - fixed)^2)
  c_term <- sum(weights) - sum(weights^2) / sum(weights)
  tau2 <- ifelse(k > 1 && c_term > 0, max(0, (q - (k - 1)) / c_term), 0)
  random_weights <- 1 / (z$se^2 + tau2)
  pooled <- sum(random_weights * z$log_tir) / sum(random_weights)
  pooled_se <- sqrt(1 / sum(random_weights))
  p <- 2 * pnorm(abs(pooled / pooled_se), lower.tail = FALSE)
  q_p <- ifelse(k > 1, pchisq(q, df = k - 1, lower.tail = FALSE), NA)
  i2 <- ifelse(q > 0 && k > 1, max(0, (q - (k - 1)) / q) * 100, 0)
  data.frame(
    k = k, log_tir = pooled, se = pooled_se, tir = exp(pooled),
    ci_low = exp(pooled - 1.96 * pooled_se), ci_high = exp(pooled + 1.96 * pooled_se),
    p_value = p, tau2 = tau2, q = q, q_p = q_p, i2_percent = i2,
    fixed_tir = exp(fixed), fixed_ci_low = exp(fixed - 1.96 / sqrt(sum(weights))),
    fixed_ci_high = exp(fixed + 1.96 / sqrt(sum(weights)))
  )
}

meta_rows <- list()
direction_rows <- list()
for (endpoint in endpoints) {
  for (model in models) {
    z <- cohort_results[cohort_results$endpoint == endpoint & cohort_results$model == model, , drop = FALSE]
    if (nrow(z) < 2) next
    pooled <- meta_dl(z)
    pooled$endpoint <- endpoint
    pooled$endpoint_label <- endpoint_labels[[endpoint]]
    pooled$model <- model
    meta_rows[[length(meta_rows) + 1]] <- pooled
    direction_rows[[length(direction_rows) + 1]] <- data.frame(
      endpoint = endpoint, model = model, cohorts = nrow(z), positive = sum(z$log_tir > 0),
      negative = sum(z$log_tir < 0), direction_concordant = length(unique(sign(z$log_tir))) == 1,
      min_tir = min(z$tir), max_tir = max(z$tir)
    )
  }
}
meta <- do.call(rbind, meta_rows)
meta$p_fdr <- p.adjust(meta$p_value, method = "BH")
meta <- meta[, c("endpoint", "endpoint_label", "model", "k", "tir", "ci_low", "ci_high", "p_value", "p_fdr",
                 "tau2", "q", "q_p", "i2_percent", "fixed_tir", "fixed_ci_low", "fixed_ci_high", "log_tir", "se")]
write.csv(meta, file.path(out, "three_cohort_random_effects_meta.csv"), row.names = FALSE, fileEncoding = "UTF-8")
direction <- do.call(rbind, direction_rows)
write.csv(direction, file.path(out, "three_cohort_direction_concordance.csv"), row.names = FALSE, fileEncoding = "UTF-8")

decisions <- data.frame(
  item = c("Harmonized effect", "Interval handling", "Pooling", "Cohorts", "Interpretation boundary"),
  decision = c(
    "Transition intensity ratio from complementary log-log panel models",
    "Offset log(interval years): HRS 4 years; CHNS 3 and 2 years; ELSA 4 or 6 years",
    "DerSimonian-Laird random effects with fixed-effect estimates also reported",
    "HRS, CHNS and ELSA analysed independently before pooling",
    "Exact transition times and within-interval intermediate states are unobserved"
  )
)
write.csv(decisions, file.path(out, "three_cohort_meta_decisions.csv"), row.names = FALSE, fileEncoding = "UTF-8")

plot_cohort <- cohort_results[, c("endpoint", "endpoint_label", "model", "cohort", "tir", "ci_low", "ci_high")]
plot_meta <- data.frame(endpoint = meta$endpoint, endpoint_label = meta$endpoint_label, model = meta$model,
                        cohort = "Pooled", tir = meta$tir, ci_low = meta$ci_low, ci_high = meta$ci_high)
plot_data <- rbind(plot_cohort, plot_meta)
plot_data$model <- factor(plot_data$model, levels = models, labels = c("Excess BMI burden", "BMI variability (VIM)"))
plot_data$cohort <- factor(plot_data$cohort, levels = c("HRS", "CHNS", "ELSA", "Pooled"))
plot_data$endpoint_label <- factor(plot_data$endpoint_label, levels = rev(unname(endpoint_labels)))

theme_set(theme_classic(base_size = 7, base_family = "sans") +
            theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
                  plot.title = element_text(size = 8, face = "bold"), strip.text = element_text(size = 7, face = "bold"),
                  legend.title = element_text(size = 6.5), legend.text = element_text(size = 6), panel.grid = element_blank()))

plot <- ggplot(plot_data, aes(tir, endpoint_label, colour = cohort, shape = cohort)) +
  geom_vline(xintercept = 1, colour = "#777777", linewidth = 0.35, linetype = 2) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.15,
                position = position_dodge(width = 0.62), linewidth = 0.45) +
  geom_point(position = position_dodge(width = 0.62), size = 1.7) +
  facet_wrap(~model, nrow = 1, scales = "free_x") +
  scale_x_log10() +
  scale_colour_manual(values = c(HRS = "#34688A", CHNS = "#B0684B", ELSA = "#4E8B57", Pooled = "#222222")) +
  labs(x = "Transition intensity ratio per 1 SD (log scale)", y = NULL, colour = NULL, shape = NULL,
       title = "HRS, CHNS and ELSA cohort-specific and pooled associations") +
  theme(legend.position = "top")

save_pub <- function(plot, stem, width_mm = 183, height_mm = 112, dpi = 600) {
  width <- width_mm / 25.4; height <- height_mm / 25.4
  svglite::svglite(file.path(out, paste0(stem, ".svg")), width = width, height = height); print(plot); dev.off()
  grDevices::cairo_pdf(file.path(out, paste0(stem, ".pdf")), width = width, height = height, family = "sans"); print(plot); dev.off()
  ragg::agg_tiff(file.path(out, paste0(stem, ".tiff")), width = width, height = height, units = "in", res = dpi, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(file.path(out, paste0(stem, "_preview.png")), width = width, height = height, units = "in", res = 180); print(plot); dev.off()
}
save_pub(plot, "figure_three_cohort_meta")

cat("Three-cohort harmonized results:\n")
print(cohort_results[, c("cohort", "endpoint", "model", "tir", "ci_low", "ci_high", "p_value")])
cat("\nThree-cohort random-effects meta-analysis:\n")
print(meta[, c("endpoint", "model", "tir", "ci_low", "ci_high", "p_value", "p_fdr", "i2_percent")])
