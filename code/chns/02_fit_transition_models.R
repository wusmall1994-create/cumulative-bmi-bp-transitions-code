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
out <- file.path(output_root, "chns_analysis")
d <- read.csv(file.path(out, "chns_intervals.csv"), fileEncoding = "UTF-8", check.names = FALSE)

d$interval <- factor(d$interval, levels = c("2006-2009", "2009-2011"))
d$province_f <- factor(d$province)
d$change_group <- factor(d$bmi_change_cat_code, levels = c(0, -1, 1), labels = c("Stable", "Loss >5%", "Gain >5%"))

endpoint_labels <- c(
  hypertension_onset = "Hypertension onset",
  treatment_initiation = "Treatment initiation",
  untreated_bp_improvement = "Untreated BP improvement",
  control_loss = "Loss of BP control",
  control_achievement = "Achievement of BP control"
)
endpoints <- names(endpoint_labels)

full_covars <- c(
  "age_start10", "female", "educ_years", "current_smoker", "current_drinker",
  "diabetes", "cvd_history", "self_health", "province_f", "urban", "interval"
)

covars_for <- function(endpoint) {
  if (endpoint == "control_loss") return(c("age_start10", "female", "urban", "interval"))
  if (endpoint == "control_achievement") return(c("age_start10", "female", "educ_years", "diabetes", "urban", "interval"))
  full_covars
}

model_data <- function(endpoint, weight_type = "ipcw", extra_filter = rep(TRUE, nrow(d))) {
  weight_var <- if (weight_type == "base") "base_weight" else "model_weight"
  needed <- unique(c(endpoint, weight_var, "commid", covars_for(endpoint),
                     "excess_bmi25_z", "bmi_vim_z", "bmi_mean_z", "bmi_auc_z",
                     "bmi_arv_z", "bmi_cv_z", "bmi_sd_z", "change_group"))
  keep <- !is.na(d[[endpoint]]) & !is.na(d[[weight_var]]) & d[[weight_var]] > 0 & extra_filter
  keep <- keep & complete.cases(d[, needed, drop = FALSE])
  z <- d[keep, , drop = FALSE]
  z$analysis_w <- if (weight_type == "none") 1 else z[[weight_var]]
  z
}

fit_rr <- function(endpoint, exposure_terms, model_name, keep_terms = exposure_terms,
                   weight_type = "ipcw", extra_filter = rep(TRUE, nrow(d))) {
  z <- model_data(endpoint, weight_type, extra_filter)
  if (nrow(z) < 60 || sum(z[[endpoint]]) < 20 || length(unique(z[[endpoint]])) < 2) return(NULL)
  des <- svydesign(ids = ~commid, weights = ~analysis_w, data = z)
  f <- as.formula(paste(endpoint, "~", paste(c(exposure_terms, covars_for(endpoint)), collapse = " + ")))
  fit <- tryCatch(svyglm(f, design = des, family = quasipoisson(link = "log")), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  co <- summary(fit)$coefficients
  rows <- list()
  for (term in rownames(co)) {
    if (!any(vapply(keep_terms, function(k) startsWith(term, k), logical(1)))) next
    b <- unname(co[term, 1]); se <- unname(co[term, 2]); p <- unname(co[term, ncol(co)])
    rows[[length(rows) + 1]] <- data.frame(
      cohort = "CHNS", endpoint = endpoint, endpoint_label = endpoint_labels[[endpoint]], model = model_name,
      term = term, log_rr = b, se = se, rr = exp(b), ci_low = exp(b - 1.96 * se), ci_high = exp(b + 1.96 * se),
      p_value = p, n_intervals = nrow(z), events = sum(z[[endpoint]]), participants = length(unique(z$person_id)),
      communities = length(unique(z$commid)), weight_type = weight_type
    )
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

main_list <- list()
for (ep in endpoints) {
  main_list[[length(main_list) + 1]] <- fit_rr(ep, "excess_bmi25_z", "Cumulative excess BMI burden", "excess_bmi25_z")
  main_list[[length(main_list) + 1]] <- fit_rr(ep, c("bmi_vim_z", "bmi_mean_z"), "BMI variability (VIM)", "bmi_vim_z")
  main_list[[length(main_list) + 1]] <- fit_rr(ep, "bmi_mean_z", "Mean BMI", "bmi_mean_z")
}
main <- do.call(rbind, Filter(Negate(is.null), main_list))
ix <- main$model %in% c("Cumulative excess BMI burden", "BMI variability (VIM)")
main$p_fdr <- NA_real_
main$p_fdr[ix] <- p.adjust(main$p_value[ix], method = "BH")
write.csv(main, file.path(out, "chns_main_transition_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

sens_list <- list()
for (ep in endpoints) {
  sens_list[[length(sens_list) + 1]] <- fit_rr(ep, "bmi_auc_z", "Alternative total BMI-years", "bmi_auc_z")
  sens_list[[length(sens_list) + 1]] <- fit_rr(ep, c("bmi_arv_z", "bmi_mean_z"), "Alternative variability ARV", "bmi_arv_z")
  sens_list[[length(sens_list) + 1]] <- fit_rr(ep, c("bmi_cv_z", "bmi_mean_z"), "Alternative variability CV", "bmi_cv_z")
  sens_list[[length(sens_list) + 1]] <- fit_rr(ep, c("bmi_sd_z", "bmi_mean_z"), "Alternative variability SD", "bmi_sd_z")
  sens_list[[length(sens_list) + 1]] <- fit_rr(ep, "excess_bmi25_z", "Exposure weight only", "excess_bmi25_z", weight_type = "base")
  sens_list[[length(sens_list) + 1]] <- fit_rr(ep, "excess_bmi25_z", "Unweighted", "excess_bmi25_z", weight_type = "none")
  sens_list[[length(sens_list) + 1]] <- fit_rr(ep, "excess_bmi25_z", "Exclude baseline CVD", "excess_bmi25_z", extra_filter = d$cvd_history == 0)
  sens_list[[length(sens_list) + 1]] <- fit_rr(ep, "excess_bmi25_z", "Exclude >5% BMI loss", "excess_bmi25_z", extra_filter = d$bmi_change_cat_code != -1)
  sens_list[[length(sens_list) + 1]] <- fit_rr(ep, c("bmi_vim_z", "bmi_mean_z"), "VIM excluding >5% BMI loss", "bmi_vim_z", extra_filter = d$bmi_change_cat_code != -1)
}
sens <- do.call(rbind, Filter(Negate(is.null), sens_list))
write.csv(sens, file.path(out, "chns_sensitivity_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

change_list <- list()
for (ep in endpoints) change_list[[length(change_list) + 1]] <- fit_rr(ep, "change_group", "BMI change category", "change_group")
change <- do.call(rbind, Filter(Negate(is.null), change_list))
change$p_fdr <- p.adjust(change$p_value, method = "BH")
write.csv(change, file.path(out, "chns_weight_change_models.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Publication figures, rendered exclusively in R.
theme_set(theme_classic(base_size = 7, base_family = "sans") +
            theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
                  plot.title = element_text(size = 8, face = "bold"), strip.text = element_text(size = 7, face = "bold"),
                  legend.title = element_text(size = 6.5), legend.text = element_text(size = 6), panel.grid = element_blank()))

tr <- read.csv(file.path(out, "chns_transition_counts.csv"), fileEncoding = "UTF-8-BOM")
tr$total_origin <- ave(tr$n, tr$interval, tr$origin, FUN = sum)
tr$probability <- ifelse(tr$total_origin > 0, tr$n / tr$total_origin, NA)
tr$origin_f <- factor(tr$origin, levels = 5:1,
                      labels = c("Treated uncontrolled", "Treated controlled", "Untreated hypertension", "Elevated untreated", "Normal untreated"))
tr$destination_f <- factor(tr$destination, levels = 1:5,
                           labels = c("Normal untreated", "Elevated untreated", "Untreated hypertension", "Treated controlled", "Treated uncontrolled"))
p_heat <- ggplot(tr, aes(destination_f, origin_f, fill = probability)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(n > 0, sprintf("%d\n%.1f%%", n, 100 * probability), "")), size = 2.0, lineheight = 0.9) +
  facet_wrap(~interval, nrow = 1) +
  scale_fill_gradient(low = "#F3F6F8", high = "#3B6C8E", labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Destination state", y = "Origin state", fill = "Row probability",
       title = "CHNS observed blood pressure and treatment-state transitions") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

forest <- subset(main, model %in% c("Cumulative excess BMI burden", "BMI variability (VIM)"))
forest$model <- factor(forest$model, levels = c("Cumulative excess BMI burden", "BMI variability (VIM)"),
                       labels = c("Excess BMI burden", "BMI variability (VIM)"))
forest$endpoint_label <- factor(forest$endpoint_label, levels = rev(unname(endpoint_labels)))
p_forest <- ggplot(forest, aes(rr, endpoint_label, colour = model, shape = model)) +
  geom_vline(xintercept = 1, colour = "#777777", linewidth = 0.35, linetype = 2) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.18,
                position = position_dodge(width = 0.45), linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.45), size = 2.0) +
  scale_x_log10() +
  scale_colour_manual(values = c("#34688A", "#B0684B")) +
  labs(x = "Transition risk ratio per 1 SD (log scale)", y = NULL, colour = NULL, shape = NULL,
       title = "CHNS replication estimates for cumulative BMI burden and variability") +
  theme(legend.position = "top")

save_pub <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(file.path(out, paste0(stem, ".svg")), width = w, height = h); print(plot); dev.off()
  grDevices::cairo_pdf(file.path(out, paste0(stem, ".pdf")), width = w, height = h, family = "sans"); print(plot); dev.off()
  ragg::agg_tiff(file.path(out, paste0(stem, ".tiff")), width = w, height = h, units = "in", res = dpi, compression = "lzw"); print(plot); dev.off()
  ragg::agg_png(file.path(out, paste0(stem, "_preview.png")), width = w, height = h, units = "in", res = 180); print(plot); dev.off()
}
save_pub(p_heat, "chns_figure1_transition_heatmap", 183, 105)
save_pub(p_forest, "chns_figure2_main_forest", 150, 100)

cat("CHNS main models:\n")
print(main[, c("endpoint", "model", "rr", "ci_low", "ci_high", "p_value", "p_fdr", "n_intervals", "events")])
