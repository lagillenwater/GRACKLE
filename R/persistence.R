## Persistence scoring (spec Section 8.3).

#' normRankWithin
#' @param x Numeric vector.
#' @param group Grouping vector (cancer type).
#' @return `(rank - (n + 1) / 2) / n` within each group, so every group is weighted by its size.
#' @export
normRankWithin <- function(x, group) {
  stats::ave(x, group, FUN = function(v) (rank(v) - (length(v) + 1) / 2) / length(v))
}

#' concordanceStat
#' @param zt,zo Normalized within-type ranks of tumor and organoid activity.
#' @return Pearson correlation, or 0 when either side has no variance.
#' @export
concordanceStat <- function(zt, zo) {
  if (length(zt) < 3 || !isTRUE(stats::sd(zt) > 0) || !isTRUE(stats::sd(zo) > 0)) return(0)
  stats::cor(zt, zo)
}

## Levels in byte order, so the order in which groups take seeded draws does not depend on the
## collation locale. A factor keeps the level order its caller chose.
radixFactor <- function(x) if (is.factor(x)) x else factor(x, levels = unique(sort(x, method = "radix")))

#' permuteWithin
#' @param strata Grouping vector.
#' @return Permuted indices that only move positions within a stratum. A size-one stratum stays fixed.
#' @export
permuteWithin <- function(strata) {
  groups <- split(seq_along(strata), radixFactor(strata))
  out <- seq_along(strata)
  out[unlist(groups, use.names = FALSE)] <- unlist(lapply(groups, function(i) i[sample.int(length(i))]), use.names = FALSE)
  out
}

#' resampleWithin
#' @param strata Grouping vector.
#' @return Bootstrap indices drawn with replacement within each stratum.
#' @export
resampleWithin <- function(strata) {
  groups <- split(seq_along(strata), radixFactor(strata))
  unlist(lapply(groups, function(i) i[sample.int(length(i), replace = TRUE)]), use.names = FALSE)
}

#' persistenceConcordance
#'
#' Stratified Spearman concordance between held-out tumor and organoid activity. The null
#' shuffles organoids among tumors within cancer type x fold.
#' @param pairs Data frame with tumor_act, organoid_act, cancer_type, fold.
#' @param n_perm Permutations.
#' @param n_boot Bootstrap resamples.
#' @param seed Random seed.
#' @return One-row data frame: conc_rho, conc_p (one-sided), conc_ci_low, conc_ci_high,
#'   detect_limit (95th percentile of the null), n_pairs.
#' @export
persistenceConcordance <- function(pairs, n_perm = 10000, n_boot = 2000, seed = 1) {
  zt <- normRankWithin(pairs$tumor_act, pairs$cancer_type)
  zo <- normRankWithin(pairs$organoid_act, pairs$cancer_type)
  obs <- concordanceStat(zt, zo)
  strata <- interaction(radixFactor(pairs$cancer_type), pairs$fold, drop = TRUE)
  draws <- withr::with_seed(seed, list(
    null = vapply(seq_len(n_perm), function(b) concordanceStat(zt, zo[permuteWithin(strata)]), numeric(1)),
    boot = vapply(seq_len(n_boot), function(b) {
      i <- resampleWithin(strata)
      concordanceStat(normRankWithin(pairs$tumor_act[i], pairs$cancer_type[i]),
                      normRankWithin(pairs$organoid_act[i], pairs$cancer_type[i]))
    }, numeric(1))))
  data.frame(conc_rho = obs,
             conc_p = (1 + sum(draws$null >= obs)) / (1 + n_perm),
             conc_ci_low = unname(stats::quantile(draws$boot, 0.025)),
             conc_ci_high = unname(stats::quantile(draws$boot, 0.975)),
             detect_limit = unname(stats::quantile(draws$null, 0.95)),
             n_pairs = nrow(pairs))
}

#' persistenceRetention
#'
#' Median over pairs of `log2((organoid + eps) / (tumor + eps))`, where eps is the 1st
#' percentile of the program's nonzero held-out activities. The one-sided test for reduction
#' flips signs within pairs.
#' @param pairs Data frame with tumor_act, organoid_act, cancer_type, fold.
#' @param n_perm Sign-flip draws.
#' @param n_boot Bootstrap resamples.
#' @param seed Random seed.
#' @return One-row data frame: ret_median_log2, ret_p, ret_ci_low, ret_ci_high, eps.
#' @export
persistenceRetention <- function(pairs, n_perm = 10000, n_boot = 2000, seed = 1) {
  act <- c(pairs$tumor_act, pairs$organoid_act)
  nz <- act[act > 0]
  eps <- if (length(nz)) unname(stats::quantile(nz, 0.01)) else 1e-8
  lr <- log2((pairs$organoid_act + eps) / (pairs$tumor_act + eps))
  obs <- stats::median(lr)
  strata <- interaction(radixFactor(pairs$cancer_type), pairs$fold, drop = TRUE)
  draws <- withr::with_seed(seed, list(
    null = vapply(seq_len(n_perm), function(b) stats::median(lr * sample(c(-1, 1), length(lr), replace = TRUE)), numeric(1)),
    boot = vapply(seq_len(n_boot), function(b) stats::median(lr[resampleWithin(strata)]), numeric(1))))
  data.frame(ret_median_log2 = obs,
             ret_p = (1 + sum(draws$null <= obs)) / (1 + n_perm),
             ret_ci_low = unname(stats::quantile(draws$boot, 0.025)),
             ret_ci_high = unname(stats::quantile(draws$boot, 0.975)),
             eps = eps)
}

#' classifyPrograms
#'
#' Spec Section 8.3 groups. Reduced: retention FDR < alpha and median log2 ratio <= -1.
#' Concordant: concordance FDR < alpha.
#' @param conc_fdr,ret_fdr,ret_median_log2 Vectors over programs.
#' @param alpha FDR level.
#' @return Character vector of group labels.
#' @export
classifyPrograms <- function(conc_fdr, ret_fdr, ret_median_log2, alpha = 0.05) {
  reduced <- ret_fdr < alpha & ret_median_log2 <= -1
  concordant <- conc_fdr < alpha
  lookup <- c("TRUE FALSE" = "drops_out", "TRUE TRUE" = "reduced_concordant",
              "FALSE TRUE" = "persists", "FALSE FALSE" = "present_not_patient_specific")
  unname(lookup[paste(reduced %in% TRUE, concordant %in% TRUE)])
}

#' concordancePower
#'
#' Power of the concordance test at a given pair structure (spec Section 8.3). Within each
#' simulated pair, tumor and organoid activity are standard normal with correlation `rho`. The
#' test is [persistenceConcordance()] at one-sided level `alpha` (nominal, before FDR).
#' @param strata Data frame with cancer_type and fold, one row per held-out pair.
#' @param rhos True correlations.
#' @param n_sims Simulations per rho.
#' @param n_perm Permutations per test.
#' @param alpha One-sided level.
#' @param seed Random seed.
#' @return Data frame with rho, power, n_pairs.
#' @export
concordancePower <- function(strata, rhos = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), n_sims = 200, n_perm = 999,
                             alpha = 0.05, seed = 1) {
  n <- nrow(strata)
  grid <- expand.grid(sim = seq_len(n_sims), rho = rhos)
  hit <- withr::with_seed(seed, vapply(seq_len(nrow(grid)), function(i) {
    z1 <- stats::rnorm(n)
    z2 <- grid$rho[i] * z1 + sqrt(1 - grid$rho[i]^2) * stats::rnorm(n)
    persistenceConcordance(data.frame(tumor_act = z1, organoid_act = z2, cancer_type = strata$cancer_type,
                                      fold = strata$fold), n_perm = n_perm, n_boot = 2, seed = i)$conc_p < alpha
  }, logical(1)))
  data.frame(rho = rhos, power = as.numeric(tapply(hit, grid$rho, mean)[as.character(rhos)]), n_pairs = n)
}
