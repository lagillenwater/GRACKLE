## Step 07: main and supplementary figures (spec Section 10).
source("analysis/hcmi_persistence/config.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
res_dir <- file.path(cfg$run, "results")
fig_dir <- file.path(res_dir, "figures"); dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
group_shapes <- c(persists = 16, drops_out = 17, reduced_concordant = 15, present_not_patient_specific = 1)
## Colors from the dataviz reference palette (categorical slots 1-2 validated with
## validate_palette.js, light mode: CVD dE 24.7, normal-vision dE 33.6); ink and grid from its chrome table.
compartment_cols <- c(`tumor-intrinsic` = "#2a78d6", microenvironment = "#eb6834")
diverging <- c(low = "#e34948", mid = "#f0efec", high = "#2a78d6")
theme_set(theme_minimal(base_size = 9) +
            theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e1e0d9", linewidth = 0.3),
                  text = element_text(colour = "#0b0b0b"), axis.text = element_text(colour = "#52514e"),
                  plot.background = element_rect(fill = "#fcfcfb", colour = NA)))
read_csv_df <- function(...) as.data.frame(fread(file.path(...)))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
mapPanel <- function(dir, title) {
  p <- read_csv_df(dir, "programs.csv")
  p <- p[p$stable %in% TRUE, , drop = FALSE]
  p$compartment <- ifelse(p$microenv %in% TRUE, "microenvironment", "tumor-intrinsic")
  p$label <- ifelse(is.na(p$top_enrichment) | p$top_enrichment == "", p$program, sub("^HALLMARK_", "", p$top_enrichment))
  ggplot(p, aes(ret_median_log2, conc_rho)) +
    geom_vline(xintercept = -1, linetype = "dashed", linewidth = 0.3, colour = "grey55") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
    geom_errorbar(aes(ymin = conc_ci_low, ymax = conc_ci_high), width = 0, linewidth = 0.25, colour = "grey60") +
    geom_point(aes(colour = compartment, shape = group), size = 2.4) +
    geom_text(aes(label = label), size = 2.2, vjust = -0.9, check_overlap = TRUE) +
    scale_colour_manual(values = compartment_cols) + scale_shape_manual(values = group_shapes, drop = FALSE) +
    labs(title = title, x = "Retention: median log2(organoid / tumor)", y = "Concordance (stratified Spearman)")
}
agreementPanel <- function() {
  a <- read_csv_df(res_dir, "compare", "map_agreement_pairs.csv")
  s <- read_csv_df(res_dir, "compare", "map_agreement.csv")
  a <- a[a$matched & a$stable_a & a$stable_b, , drop = FALSE]
  a$call <- ifelse(a$group_a == a$group_b, "same group", "different group")
  ggplot(a, aes(conc_a, conc_b, colour = call)) + geom_abline(linewidth = 0.3, colour = "grey75") + geom_point(size = 2) +
    scale_colour_manual(values = c(`same group` = "#2a78d6", `different group` = "#eb6834")) +
    labs(title = sprintf("2D vs 3D co-essentiality: kappa %.2f (n = %d)", s$kappa, s$n_matched_stable),
         x = "Concordance, coess_2d", y = "Concordance, coess_3d", colour = NULL)
}
heatmapPanel <- function(arm) {
  d <- read_csv_df(res_dir, arm, "per_type.csv")
  ggplot(d, aes(cancer_type, program, fill = rho)) + geom_tile() +
    scale_fill_gradient2(low = diverging[["low"]], mid = diverging[["mid"]], high = diverging[["high"]], midpoint = 0, limits = c(-1, 1)) +
    labs(title = paste(arm, "concordance by cancer type"), x = NULL, y = NULL)
}

save_fig((mapPanel(file.path(res_dir, "coess_2d"), "coess_2d") | mapPanel(file.path(res_dir, "coess_3d"), "coess_3d") |
            agreementPanel()) + plot_layout(guides = "collect"), "main_persistence_maps", 13, 4.5)
save_fig(heatmapPanel("coess_2d") | heatmapPanel("coess_3d"), "main_per_type_heatmaps", 10, 6)
save_fig(mapPanel(file.path(res_dir, "string_mech"), "string_mech") | mapPanel(file.path(res_dir, "string_all"), "string_all"),
         "supp_string_maps", 9, 4.5)
tr <- read_csv_df(res_dir, "compare", "transfer.csv")
save_fig(ggplot(tr, aes(arm, diff_mean)) + geom_hline(yintercept = 0, colour = "grey75") +
           geom_pointrange(aes(ymin = diff_mean - 1.96 * diff_se, ymax = diff_mean + 1.96 * diff_se)) +
           labs(x = NULL, y = "Tumor-to-organoid error: selected minus no gene prior"), "supp_transfer", 5, 3.5)
invisible(lapply(cfg$arms, function(arm) {
  m <- read_csv_df(res_dir, arm, "maps_all_settings.csv")
  save_fig(ggplot(m, aes(ret_median_log2, conc_rho, shape = group)) + geom_point(size = 1.1) + facet_wrap(~setting_id) +
             scale_shape_manual(values = group_shapes) + labs(title = paste(arm, ": every alpha setting")),
           paste0("supp_alpha_grid_", arm), 12, 10)
  sens <- file.path(res_dir, arm, c("sensitivity_k0.5", "sensitivity_k2"))
  sens <- sens[file.exists(file.path(sens, "programs.csv"))]
  if (length(sens)) save_fig(Reduce(`|`, lapply(sens, function(d) mapPanel(d, paste(arm, basename(d))))),
                             paste0("supp_rank_sensitivity_", arm), 4.5 * length(sens), 4.5)
}))
pw <- read_csv_df(res_dir, "compare", "power_curve.csv")
save_fig(ggplot(pw, aes(rho, power)) + geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey60") +
           geom_line() + geom_point() + labs(x = "True concordance (rho)", y = "Power, one-sided alpha 0.05",
                                             title = sprintf("Concordance test power at %d held-out pairs", pw$n_pairs[1])),
         "supp_power_curve", 5, 3.5)
invisible(lapply(cfg$arms, function(arm) {
  p <- read_csv_df(res_dir, arm, "programs.csv")
  p$program <- factor(p$program, levels = p$program[order(p$ref)])
  save_fig(ggplot(p, aes(program, n_matched, fill = stable)) + geom_col() + geom_hline(yintercept = 4, linetype = "dashed") +
             scale_fill_manual(values = c(`TRUE` = "#2a78d6", `FALSE` = "#c3c2b7")) +
             labs(title = paste(arm, ": folds matched per reference program"), x = NULL, y = "Folds matched (of 5)") +
             theme(axis.text.x = element_text(angle = 90, size = 5)),
           paste0("supp_stability_", arm), 10, 3.5)
}))
pa <- read_csv_df(cfg$processed, "k_parallel_analysis.csv")
save_fig(ggplot(pa, aes(component)) + geom_line(aes(y = observed, colour = "observed")) +
           geom_line(aes(y = threshold, colour = "null, 95th percentile")) + scale_y_log10() +
           labs(y = "Singular value", colour = NULL,
                title = sprintf("Parallel analysis: k = %d", readRDS(file.path(cfg$processed, "k.rds")))),
         "supp_parallel_analysis", 6, 4)
message("figures written to ", fig_dir)
