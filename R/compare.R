#' transferTable
#'
#' Tumor-to-organoid transfer reconstruction per arm (spec Section 8.4). The gene prior's
#' effect is the paired per-fold difference between the selected setting and alpha_2 = 0.
#' @param results Transfer-task results (meta: arm, variant, fold; fit: recon_test).
#' @return Data frame with arm, selected_mean/se, baseline_mean/se, diff_mean/se, n_folds.
#' @export
transferTable <- function(results) {
  d <- do.call(rbind, lapply(results, function(r) data.frame(
    arm = r$meta$arm, variant = r$meta$variant, fold = r$meta$fold,
    recon = if (is.null(r$fit)) NA_real_ else r$fit$recon_test, stringsAsFactors = FALSE)))
  paired <- merge(d[d$variant == "selected", ], d[d$variant == "baseline", ],
                  by = c("arm", "fold"), suffixes = c("_selected", "_baseline"))
  paired$diff <- paired$recon_selected - paired$recon_baseline
  se <- function(x) stats::sd(x) / sqrt(length(x))
  arms <- sort(unique(paired$arm))
  agg <- function(x, f) as.numeric(tapply(x, paired$arm, f)[arms])
  data.frame(arm = arms,
             selected_mean = agg(paired$recon_selected, mean), selected_se = agg(paired$recon_selected, se),
             baseline_mean = agg(paired$recon_baseline, mean), baseline_se = agg(paired$recon_baseline, se),
             diff_mean = agg(paired$diff, mean), diff_se = agg(paired$diff, se),
             n_folds = as.integer(table(paired$arm)[arms]), stringsAsFactors = FALSE)
}

#' cohenKappa
#' @param a,b Equal-length label vectors.
#' @return Cohen's kappa (NA for empty input).
#' @export
cohenKappa <- function(a, b) {
  if (!length(a)) return(NA_real_)
  lv <- union(a, b)
  tab <- table(factor(a, lv), factor(b, lv))
  n <- sum(tab)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / n^2
  if (pe == 1) return(1)
  (po - pe) / (1 - pe)
}

#' mapAgreement
#'
#' Matches two arms' reference programs (Hungarian method on H cosine, threshold `tau`) and
#' compares group calls and concordance among programs that are matched and stable in both.
#' @param H_a,H_b Reference H of each arm (same genes, same k).
#' @param programs_a,programs_b Program tables with program, stable, group, conc_rho.
#' @param tau Matching threshold.
#' @return A list with pairs (one row per program of arm a) and summary.
#' @export
mapAgreement <- function(H_a, H_b, programs_a, programs_b, tau) {
  m <- matchPrograms(H_a, H_b, tau)
  m$program_a <- rownames(H_a)[m$ref]
  m$program_b <- rownames(H_b)[m$fold_prog]
  ga <- programs_a[match(m$program_a, programs_a$program), , drop = FALSE]
  gb <- programs_b[match(m$program_b, programs_b$program), , drop = FALSE]
  pairs <- data.frame(m, stable_a = ga$stable %in% TRUE, stable_b = gb$stable %in% TRUE,
                      group_a = ga$group, group_b = gb$group, conc_a = ga$conc_rho, conc_b = gb$conc_rho,
                      stringsAsFactors = FALSE)
  sub <- pairs[pairs$matched & pairs$stable_a & pairs$stable_b & !is.na(pairs$group_a) & !is.na(pairs$group_b), , drop = FALSE]
  list(pairs = pairs,
       summary = data.frame(n_matched_stable = nrow(sub),
                            agreement = if (nrow(sub)) mean(sub$group_a == sub$group_b) else NA_real_,
                            kappa = cohenKappa(sub$group_a, sub$group_b),
                            concordance_cor = if (nrow(sub) > 2) stats::cor(sub$conc_a, sub$conc_b) else NA_real_,
                            unmatched_a = paste(m$program_a[!m$matched], collapse = ";"),
                            unmatched_b = paste(setdiff(rownames(H_b), m$program_b[m$matched]), collapse = ";"),
                            stringsAsFactors = FALSE))
}
