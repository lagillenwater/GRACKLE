## Pair alignment: how well a factorization places each organoid with its own tumor, how many of its
## programs read as named pathways, and the balance between the two.

#' pairingScore
#'
#' Where each organoid's own tumor falls in a hierarchical tree of the samples' program activity,
#' without cutting the tree. Activities are log1p-transformed and each program standardized; the
#' tree is Ward linkage on Euclidean distance. For each organoid-tumor pair, the score is the share
#' of the other tumors that join the organoid later than its own tumor does (ties count half): 1 when
#' the own tumor joins first, about 0.5 by chance.
#' @param W Samples x programs activity; rownames are profile ids.
#' @param samples Sample table with profile_id, case_id, sample_id, parent_sample_id and role.
#' @return Data frame with organoid, tumor and score, one row per pair among the rows of `W`.
#' @export
pairingScore <- function(W, samples) {
  ids <- rownames(W)
  s <- samples[match(ids, samples$profile_id), , drop = FALSE]
  pairs <- pairTable(s)
  tumors <- ids[s$role %in% "tumor"]
  if (!nrow(pairs) || length(tumors) < 2) return(data.frame(organoid = character(0), tumor = character(0), score = numeric(0)))
  Z <- log1p(W)
  Z <- scale(Z[, apply(Z, 2, stats::sd) > 0, drop = FALSE])
  D <- as.matrix(stats::cophenetic(stats::hclust(stats::dist(Z), method = "ward.D2")))
  d_own <- D[cbind(pairs$organoid, pairs$tumor)]
  Dt <- D[pairs$organoid, tumors, drop = FALSE]
  closer <- rowSums(Dt < d_own)
  ties <- rowSums(Dt == d_own) - 1   # the own tumor is not a competitor
  data.frame(organoid = pairs$organoid, tumor = pairs$tumor,
             score = 1 - (closer + 0.5 * ties) / (length(tumors) - 1), stringsAsFactors = FALSE)
}

#' enrichedShare
#'
#' Share of programs whose top genes are enriched for at least one gene set ([enrichHypergeom()]).
#' @param H Programs x genes loadings.
#' @param gene_sets Named list of gene symbol vectors.
#' @param fdr BH threshold.
#' @param frac Fraction of genes used as each program's top set.
#' @return A number between 0 and 1.
#' @export
enrichedShare <- function(H, gene_sets, fdr = 0.05, frac = 0.05) {
  e <- enrichHypergeom(H, gene_sets, frac = frac)
  best <- tapply(e$fdr, e$program, min)
  mean(best[rownames(H)] < fdr)
}

#' balanceSelect
#'
#' Picks the setting that balances pairing and pathway enrichment. A setting must pair organoids
#' better than the baseline (plain NMF); among those, settings that no other beats on both measures
#' are kept, and the one closest to the ideal corner wins, with both measures rescaled to 0-1 over
#' the eligible settings and weighted equally.
#' @param tab Data frame with setting_id, pairing, enrichment, baseline (one TRUE row) and eligible.
#' @return `tab` with floor_ok, pareto, distance and selected.
#' @export
balanceSelect <- function(tab) {
  ok <- tab$eligible & !tab$baseline & tab$pairing > tab$pairing[tab$baseline]
  rescale <- function(x) {
    r <- range(x[tab$eligible])
    if (diff(r) > 0) (x - r[1]) / diff(r) else rep(1, length(x))
  }
  p <- rescale(tab$pairing)
  e <- rescale(tab$enrichment)
  dominated <- vapply(seq_len(nrow(tab)), function(i) any(ok & p >= p[i] & e >= e[i] & (p > p[i] | e > e[i])),
                      logical(1))
  tab$floor_ok <- ok
  tab$pareto <- ok & !dominated
  tab$distance <- sqrt((1 - p)^2 + (1 - e)^2)
  tab$selected <- FALSE
  if (any(tab$pareto)) tab$selected[which(tab$pareto)[which.min(tab$distance[tab$pareto])]] <- TRUE
  tab
}
