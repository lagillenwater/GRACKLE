## Gene graphs for GRACKLE (spec Section 4.5): STRING arms and co-essentiality arms.

.datatable.aware <- TRUE
.string_prior <- 0.041

#' readStringLinks
#' @param path STRING `protein.links.full` file (plain or .gz).
#' @return Data frame with protein1, protein2, the channel columns and combined_score.
#' @export
readStringLinks <- function(path) {
  data.table::fread(path, sep = " ", header = TRUE, data.table = FALSE)
}

#' readStringAliases
#' @param path STRING `protein.aliases` file (plain or .gz).
#' @return Data frame with columns string_protein_id, alias, source.
#' @export
readStringAliases <- function(path) {
  x <- data.table::fread(path, sep = "\t", header = TRUE, quote = "", data.table = FALSE)
  names(x) <- c("string_protein_id", "alias", "source")
  x
}

#' stringProteinToSymbol
#'
#' Maps STRING proteins to HGNC symbols through their Ensembl gene ids, using the pinned HGNC table.
#' @param aliases Output of [readStringAliases()].
#' @param hgnc Data frame with columns `symbol` and `ensembl_gene_id`.
#' @return Unique data frame with columns string_protein_id and symbol.
#' @export
stringProteinToSymbol <- function(aliases, hgnc) {
  a <- aliases[aliases$source == "Ensembl_gene", c("string_protein_id", "alias")]
  a$symbol <- hgnc$symbol[match(a$alias, hgnc$ensembl_gene_id)]
  out <- unique(a[!is.na(a$symbol), c("string_protein_id", "symbol")])
  rownames(out) <- NULL
  out
}

priorAway <- function(score, prior = .string_prior) {
  (pmax(score, prior) - prior) / (1 - prior)
}

#' stringRecombine
#'
#' Recombines selected STRING evidence channels with STRING's procedure: remove the prior
#' from each channel, combine direct and transferred scores per channel, combine channels as
#' `1 - prod(1 - s)`, then add the prior back once.
#' @param links Data frame from [readStringLinks()] (channels on the 0-1000 scale).
#' @param channels Channel names; each `<channel>_transferred` column is included when present.
#' @param prior STRING's prior. Default 0.041.
#' @return Numeric vector in [0, 1]; 0 where no selected channel has evidence.
#' @export
stringRecombine <- function(links, channels = c("experiments", "database"), prior = .string_prior) {
  per_channel <- vapply(channels, function(ch) {
    direct <- priorAway(links[[ch]] / 1000, prior)
    tr <- paste0(ch, "_transferred")
    transferred <- if (tr %in% names(links)) priorAway(links[[tr]] / 1000, prior) else 0
    1 - (1 - direct) * (1 - transferred)
  }, numeric(nrow(links)))
  per_channel <- matrix(per_channel, nrow = nrow(links))
  combined <- (1 - matrixStats::rowProds(1 - per_channel)) * (1 - prior) + prior
  ifelse(rowSums(per_channel > 0) > 0, combined, 0)
}

#' stringCombinedAll
#'
#' Recomputes STRING's combined_score from all channels. The channel scores in the published
#' `protein.links.full` file already carry STRING's homology correction (checked against the v12
#' file: re-applying it moves scores by up to 749), so none is applied here. Used to validate
#' [stringRecombine()].
#' @param links Data frame from [readStringLinks()].
#' @param prior STRING's prior. Default 0.041.
#' @return Integer scores on the 0-1000 scale.
#' @export
stringCombinedAll <- function(links, prior = .string_prior) {
  pa <- function(col) priorAway(links[[col]] / 1000, prior)
  both <- function(a, b) 1 - (1 - pa(a)) * (1 - pa(b))
  parts <- cbind(both("neighborhood", "neighborhood_transferred"),
                 pa("fusion"),
                 pa("cooccurence"),
                 both("coexpression", "coexpression_transferred"),
                 both("experiments", "experiments_transferred"),
                 both("database", "database_transferred"),
                 both("textmining", "textmining_transferred"))
  combined <- (1 - matrixStats::rowProds(1 - parts)) * (1 - prior) + prior
  floor(combined * 1000)
}

#' stringGraph
#'
#' Gene-by-gene STRING graph. `mode = "mech"` uses the experiments and database channels
#' recombined by [stringRecombine()]; `mode = "all"` uses combined_score / 1000. Proteins
#' that map to the same gene pair are aggregated by max. Genes without edges have zero rows.
#' @param links Data frame from [readStringLinks()].
#' @param protein_map Output of [stringProteinToSymbol()].
#' @param genes Character vector of HGNC symbols giving the graph's rows and columns.
#' @param mode "mech" or "all".
#' @return Symmetric genes x genes matrix with a zero diagonal.
#' @export
stringGraph <- function(links, protein_map, genes, mode = c("mech", "all")) {
  mode <- match.arg(mode)
  map <- protein_map[protein_map$symbol %in% genes, , drop = FALSE]
  keep <- links$protein1 %in% map$string_protein_id & links$protein2 %in% map$string_protein_id
  links <- links[keep, , drop = FALSE]
  G <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  if (!nrow(links)) return(G)
  w <- if (mode == "mech") stringRecombine(links) else links$combined_score / 1000
  e <- data.table::data.table(p1 = links$protein1, p2 = links$protein2, w = w)
  e <- e[e$w > 0]
  m <- data.table::data.table(string_protein_id = map$string_protein_id, symbol = map$symbol)
  e <- merge(e, m, by.x = "p1", by.y = "string_protein_id", allow.cartesian = TRUE)
  data.table::setnames(e, "symbol", "g1")
  e <- merge(e, m, by.x = "p2", by.y = "string_protein_id", allow.cartesian = TRUE)
  data.table::setnames(e, "symbol", "g2")
  e <- e[e$g1 != e$g2]
  e <- e[, list(w = max(w)), by = c("g1", "g2")]
  G[cbind(match(e$g1, genes), match(e$g2, genes))] <- e$w
  G <- pmax(G, t(G))
  diag(G) <- 0
  G
}

#' graphSummary
#' @param G Symmetric gene graph.
#' @return One-row data frame: n_genes_with_edges, n_edges, mean_degree (over genes with edges).
#' @export
graphSummary <- function(G) {
  deg <- rowSums(G)
  data.frame(n_genes_with_edges = sum(deg > 0),
             n_edges = sum(G[upper.tri(G)] > 0),
             mean_degree = if (any(deg > 0)) mean(deg[deg > 0]) else 0)
}

utils::globalVariables("w")

## Co-essentiality arms ---------------------------------------------------------------

#' geneArm
#' @param location HGNC cytogenetic locations, for example "17q21.31".
#' @return Chromosome arm ("17q"), or NA when the location has no single arm.
#' @export
geneArm <- function(location) {
  m <- regmatches(location, regexec("^([0-9]{1,2}|X|Y)([pq])", location))
  vapply(m, function(x) if (length(x) == 3) paste0(x[2], x[3]) else NA_character_, character(1))
}

.hcmi_pattern <- "HCM-[A-Z]{4}-[0-9]{4}-[A-Z][0-9]{2}"

#' hcmiToken
#' @param x Character vector that may contain an HCMI identifier.
#' @return The case-level HCMI id (for example "HCM-BROD-0577-C15"), or NA.
#' @export
hcmiToken <- function(x) {
  out <- rep(NA_character_, length(x))
  hit <- !is.na(x) & grepl(.hcmi_pattern, x)
  out[hit] <- regmatches(x[hit], regexpr(.hcmi_pattern, x[hit]))
  out
}

#' informativeGenes
#'
#' Genes whose dependency varies across models: SD across models above the `quantile` of the
#' SDs of DepMap's non-essential control genes.
#' @param effect Genes x models matrix of gene effects (Chronos).
#' @param nonessential HGNC symbols of the non-essential controls.
#' @param quantile Control quantile used as the threshold. Default 0.95.
#' @return Character vector of informative genes.
#' @export
informativeGenes <- function(effect, nonessential, quantile = 0.95) {
  sds <- matrixStats::rowSds(effect, na.rm = TRUE)
  names(sds) <- rownames(effect)
  ctrl <- sds[intersect(nonessential, names(sds))]
  if (length(ctrl) < 10) stop("informativeGenes: fewer than 10 non-essential controls found in the effect matrix")
  thr <- stats::quantile(ctrl, quantile, na.rm = TRUE)
  names(sds)[!is.na(sds) & sds > thr]
}

#' glsWhiten
#'
#' Whitens the model dimension by `Sigma^(-1/2)`, where `Sigma` is the model-by-model
#' covariance estimated across all genes (Wainberg et al. 2021). Missing values are replaced
#' by the gene's mean. Eigenvalues are floored at `ridge * max`.
#' @param effect Genes x models matrix.
#' @param ridge Relative eigenvalue floor. Default 1e-8.
#' @return Whitened genes x models matrix with attribute `"intercept"`: the whitened
#'   all-ones vector (`colSums(Sigma^(-1/2))`), needed by [glsCorrelation()].
#' @export
glsWhiten <- function(effect, ridge = 1e-8) {
  E <- effect
  na <- is.na(E)
  if (any(na)) E[na] <- rowMeans(E, na.rm = TRUE)[row(E)[na]]
  ev <- eigen(stats::cov(E), symmetric = TRUE)
  vals <- pmax(ev$values, ridge * max(ev$values))
  W_half <- ev$vectors %*% (t(ev$vectors) / sqrt(vals))
  out <- E %*% W_half
  dimnames(out) <- dimnames(effect)
  attr(out, "intercept") <- colSums(W_half)
  out
}

#' glsCorrelation
#'
#' GLS correlation between genes: the cosine of whitened profiles after projecting out the
#' whitened intercept. This is the GLS analogue of a Pearson correlation.
#' @param Ew Output of [glsWhiten()] (optionally subset to some genes; keep the attribute).
#' @return Genes x genes similarity matrix.
#' @export
glsCorrelation <- function(Ew) {
  v <- attr(Ew, "intercept")
  R <- Ew - outer(as.vector(Ew %*% v) / sum(v^2), v)
  nr <- sqrt(rowSums(R^2))
  nr[nr == 0] <- 1
  C <- tcrossprod(R / nr)
  dimnames(C) <- list(rownames(Ew), rownames(Ew))
  C
}

#' coessGraph
#'
#' Co-essentiality gene graph (spec Section 4.5): informative genes only, GLS correlation,
#' same-arm pairs and negative values set to 0, `k_nn` nearest neighbors among `genes`.
#' @param effect Genes x models gene-effect matrix for the arm's models.
#' @param genes Final gene set giving the graph's rows and columns.
#' @param nonessential Non-essential control genes.
#' @param arm Named character vector: gene symbol to chromosome arm.
#' @param k_nn Neighbors kept per gene.
#' @return Symmetric genes x genes matrix with a zero diagonal.
#' @export
coessGraph <- function(effect, genes, nonessential, arm, k_nn = 10) {
  G <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  use <- intersect(genes, informativeGenes(effect, nonessential))
  if (length(use) < 2) return(G)
  Ew <- glsWhiten(effect)
  sub <- Ew[use, , drop = FALSE]
  attr(sub, "intercept") <- attr(Ew, "intercept")
  C <- glsCorrelation(sub)
  a <- unname(arm[use])
  same_arm <- outer(a, a, "==")
  same_arm[is.na(same_arm)] <- FALSE
  C[same_arm] <- 0
  C[is.na(C) | C < 0] <- 0
  diag(C) <- 0
  G[use, use] <- knnSparsify(C, k_nn)
  G
}

#' dropHeldOutModels
#' @param effect Genes x models matrix.
#' @param model_case Named character vector: model id (column of `effect`) to HCMI case id (NA if not HCMI).
#' @param heldout_cases Case ids held out in this fold.
#' @return `effect` without the columns of held-out cases' models.
#' @export
dropHeldOutModels <- function(effect, model_case, heldout_cases) {
  drop <- names(model_case)[!is.na(model_case) & model_case %in% heldout_cases]
  effect[, setdiff(colnames(effect), drop), drop = FALSE]
}
