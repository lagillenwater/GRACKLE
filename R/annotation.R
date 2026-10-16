## Post-hoc annotation of programs (spec Section 9).

.nonsilent <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
                "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site")

#' topGenes
#' @param h Named loading vector.
#' @param frac Fraction kept. Default 0.05.
#' @return Names of the top `floor(frac * length(h))` genes (at least one).
#' @export
topGenes <- function(h, frac = 0.05) {
  names(sort(h, decreasing = TRUE))[seq_len(max(1, floor(frac * length(h))))]
}

#' enrichHypergeom
#'
#' Hypergeometric enrichment of each program's top genes, with the program's gene set
#' (colnames(H)) as background. Vectorized over programs and gene sets.
#' @param H Programs x genes loadings.
#' @param gene_sets Named list of gene symbol vectors.
#' @param frac Fraction of genes used as each program's top set.
#' @return Data frame with program, gene_set, overlap, set_size, n_top, p, fdr (BH over all tests).
#' @export
enrichHypergeom <- function(H, gene_sets, frac = 0.05) {
  bg <- colnames(H)
  sets <- lapply(gene_sets, intersect, bg)
  sets <- sets[lengths(sets) > 0]
  M <- vapply(sets, function(s) bg %in% s, logical(length(bg)))                 # genes x sets
  M <- matrix(M, nrow = length(bg), dimnames = list(bg, names(sets)))
  Tm <- t(apply(H, 1, function(h) bg %in% topGenes(stats::setNames(h, bg), frac)))   # programs x genes
  Tm <- matrix(Tm, nrow = nrow(H), dimnames = list(rownames(H), bg))
  overlap <- Tm %*% M
  size <- matrix(colSums(M), nrow(overlap), ncol(overlap), byrow = TRUE)
  n_top <- rowSums(Tm)
  p <- stats::phyper(overlap - 1, size, length(bg) - size, n_top, lower.tail = FALSE)
  out <- data.frame(program = rep(rownames(H), times = ncol(M)),
                    gene_set = rep(colnames(M), each = nrow(H)),
                    overlap = as.vector(overlap), set_size = as.vector(size),
                    n_top = rep(n_top, times = ncol(M)), p = as.vector(p), stringsAsFactors = FALSE)
  out$fdr <- stats::p.adjust(out$p, "BH")
  out
}

#' programAssociation
#'
#' Spearman correlation of each program's activity with a covariate. Columns with too few
#' complete values or no variance get rho = 0 and p = 1.
#' @param act Samples x programs activity.
#' @param covariate Numeric vector aligned with rows.
#' @param label Name recorded in the `covariate` column.
#' @return Data frame with program, covariate, rho, p, fdr.
#' @export
programAssociation <- function(act, covariate, label) {
  res <- vapply(seq_len(ncol(act)), function(j) {
    a <- act[, j]
    ok <- !is.na(a) & !is.na(covariate)
    if (sum(ok) < 5 || !isTRUE(stats::sd(a[ok]) > 0) || !isTRUE(stats::sd(covariate[ok]) > 0)) return(c(0, 1))
    ct <- suppressWarnings(stats::cor.test(a[ok], covariate[ok], method = "spearman", exact = FALSE))
    c(unname(ct$estimate), ct$p.value)
  }, numeric(2))
  data.frame(program = colnames(act), covariate = label, rho = res[1, ], p = res[2, ],
             fdr = stats::p.adjust(res[2, ], "BH"), stringsAsFactors = FALSE)
}

#' estimateScores
#'
#' ESTIMATE stromal and immune scores (tidyestimate), used for annotation only.
#' @param expr Genes x samples expression (log2 TPM) with HGNC symbols as rownames.
#' @return Data frame with profile_id, stromal, immune, estimate.
#' @export
estimateScores <- function(expr) {
  if (!requireNamespace("tidyestimate", quietly = TRUE)) stop("estimateScores needs the tidyestimate package")
  df <- as.data.frame(expr)
  common <- tidyestimate::filter_common_genes(df, id = "hgnc_symbol", tidy = FALSE,
                                              tell_missing = FALSE, find_alias = FALSE)
  sc <- as.data.frame(tidyestimate::estimate_score(common, is_affymetrix = FALSE))
  data.frame(profile_id = sc$sample, stromal = sc$stromal, immune = sc$immune, estimate = sc$estimate,
             stringsAsFactors = FALSE)
}

#' mutationMatrix
#' @param mut Data frame with sample_id, Hugo_Symbol, Variant_Classification.
#' @param row_ids Row names of the result (profile ids).
#' @param row_sample Sample id of each row.
#' @param genes Columns (for example OncoKB genes).
#' @param profiled_samples Sequenced sample ids.
#' @return Logical matrix: TRUE for a non-silent mutation; NA for rows that were not sequenced.
#' @export
mutationMatrix <- function(mut, row_ids, row_sample, genes, profiled_samples) {
  m <- mut[mut$Variant_Classification %in% .nonsilent & mut$Hugo_Symbol %in% genes &
             mut$sample_id %in% row_sample, , drop = FALSE]
  M <- matrix(FALSE, length(row_ids), length(genes), dimnames = list(row_ids, genes))
  M[cbind(match(m$sample_id, row_sample), match(m$Hugo_Symbol, genes))] <- TRUE
  M[!(row_sample %in% profiled_samples), ] <- NA
  M
}

#' genotypeAssociation
#'
#' Wilcoxon test of program activity in mutated vs wild-type tumors, for genes with at least
#' `min_mut` mutated and at least 2 wild-type tumors.
#' @param act Tumors x programs.
#' @param mut Tumors x genes logical (rows aligned with `act`).
#' @param min_mut Minimum number of mutated tumors. Default 10.
#' @return Data frame with program, gene, n_mut, delta (median difference), p, fdr.
#' @export
genotypeAssociation <- function(act, mut, min_mut = 10) {
  n_mut <- colSums(mut, na.rm = TRUE)
  n_wt <- colSums(!mut, na.rm = TRUE)
  mut <- mut[, n_mut >= min_mut & n_wt >= 2, drop = FALSE]
  if (!ncol(mut)) return(data.frame(program = character(0), gene = character(0), n_mut = integer(0),
                                    delta = numeric(0), p = numeric(0), fdr = numeric(0)))
  grid <- expand.grid(program = colnames(act), gene = colnames(mut), stringsAsFactors = FALSE)
  res <- mapply(function(p, g) {
    m <- mut[, g]; a <- act[, p]; ok <- !is.na(m) & !is.na(a)
    w <- suppressWarnings(stats::wilcox.test(a[ok & m], a[ok & !m], exact = FALSE))
    c(stats::median(a[ok & m]) - stats::median(a[ok & !m]), w$p.value)
  }, grid$program, grid$gene)
  data.frame(grid, n_mut = as.integer(colSums(mut, na.rm = TRUE)[grid$gene]), delta = res[1, ], p = res[2, ],
             fdr = stats::p.adjust(res[2, ], "BH"), row.names = NULL, stringsAsFactors = FALSE)
}
