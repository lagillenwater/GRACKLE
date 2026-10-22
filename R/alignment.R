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

#' gseaPrograms
#'
#' Preranked GSEA (fgsea, one-sided for high loadings) of every program against every gene set, using
#' each gene's loading rather than a top-gene cut-off. Seeded per program, so results do not depend on
#' `cores`.
#' @param H Programs x genes loadings.
#' @param gene_sets Named list of gene symbol vectors.
#' @param min_size,max_size Gene-set size limits (genes present in `H`).
#' @param seed Random seed.
#' @param cores Parallel workers over programs.
#' @return Data frame with program, gene_set, nes, p and fdr (BH over all programs and sets).
#' @export
gseaPrograms <- function(H, gene_sets, min_size = 10, max_size = 500, seed = 1, cores = 1) {
  run <- function(i) {
    x <- stats::setNames(H[i, ], colnames(H))
    r <- NULL
    utils::capture.output(r <- withr::with_seed(seed + i, suppressWarnings(
      fgsea::fgsea(gene_sets, x, scoreType = "pos", minSize = min_size, maxSize = max_size, nproc = 1))))
    data.frame(program = rep(rownames(H)[i], nrow(r)), gene_set = r$pathway, nes = r$NES, p = r$pval,
               stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, parallel::mclapply(seq_len(nrow(H)), run, mc.cores = cores))
  out$fdr <- stats::p.adjust(out$p, "BH")
  out
}

#' enrichedShare
#'
#' Share of programs enriched for at least one gene set: by preranked GSEA ([gseaPrograms()], the
#' default) or by over-representation of each program's top genes ([enrichHypergeom()]).
#' @param H Programs x genes loadings.
#' @param gene_sets Named list of gene symbol vectors.
#' @param fdr BH threshold.
#' @param frac Fraction of genes used as each program's top set (method "ora").
#' @param method "gsea" or "ora".
#' @param cores Parallel workers (method "gsea").
#' @return A number between 0 and 1.
#' @export
enrichedShare <- function(H, gene_sets, fdr = 0.05, frac = 0.05, method = c("gsea", "ora"), cores = 1) {
  method <- match.arg(method)
  e <- if (method == "gsea") gseaPrograms(H, gene_sets, cores = cores) else enrichHypergeom(H, gene_sets, frac = frac)
  best <- tapply(e$fdr, e$program, min, na.rm = TRUE)
  mean(best[rownames(H)] < fdr, na.rm = FALSE)
}

#' enrichmentStrength
#'
#' Median over programs of each program's best normalized enrichment score from [gseaPrograms()]. A
#' graded measure: the share of enriched programs saturates under GSEA.
#' @param g Output of [gseaPrograms()].
#' @return A number.
#' @export
enrichmentStrength <- function(g) {
  stats::median(tapply(g$nes, g$program, max, na.rm = TRUE))
}

#' controlHolds
#'
#' The persistence map's positive control (spec Section 9): among stable programs, microenvironment
#' programs must drop out or be reduced more often than tumor-intrinsic ones. A sample graph strong
#' enough to fold tumor-only programs into shared ones fails it.
#' @param programs Program table from [combineArm()] with stable, microenv and group.
#' @return TRUE or FALSE.
#' @export
controlHolds <- function(programs) {
  p <- programs[programs$stable %in% TRUE, , drop = FALSE]
  lost <- p$group %in% c("drops_out", "reduced_concordant")
  m <- p$microenv %in% TRUE
  isTRUE(mean(lost[m]) > mean(lost[!m]))   # FALSE when either kind is absent
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
