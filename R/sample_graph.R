## Sample graph for GRACKLE: pairing, patient annotations and tumor purity (spec Section 6).

#' knnSparsify
#'
#' Keeps each row's `k_nn` largest positive similarities, then symmetrizes with `pmax`.
#' Ties are broken by `secondary` (smaller first, NA last) and then by row name, so the
#' result does not depend on row order.
#'
#' @param S Square similarity matrix with row and column names.
#' @param k_nn Neighbors kept per row.
#' @param secondary Optional square matrix of tie-break distances, aligned with `S`.
#' @return Symmetric matrix with a zero diagonal.
#' @export
knnSparsify <- function(S, k_nn, secondary = NULL) {
  n <- nrow(S)
  ids <- rownames(S)
  if (is.null(ids)) ids <- as.character(seq_len(n))
  id_rank <- rank(ids, ties.method = "first")
  pick <- lapply(seq_len(n), function(i) {
    s <- S[i, ]
    s[i] <- -Inf
    s[is.na(s)] <- -Inf
    tb <- if (is.null(secondary)) rep(0, n) else {
      x <- secondary[i, ]
      x[is.na(x)] <- Inf
      x
    }
    o <- order(-s, tb, id_rank)
    utils::head(o[s[o] > 0], k_nn)
  })
  out <- matrix(0, n, n, dimnames = dimnames(S))
  idx <- cbind(rep(seq_len(n), lengths(pick)), unlist(pick, use.names = FALSE))
  out[idx] <- S[idx]
  out <- pmax(out, t(out))
  diag(out) <- 0
  out
}

#' scaleMeanDegree
#'
#' Divides a similarity matrix by its mean degree over nodes with at least one edge.
#' @param S Square similarity matrix.
#' @return `S` scaled so that its mean nonzero degree is 1 (unchanged if it has no edges).
#' @export
scaleMeanDegree <- function(S) {
  deg <- rowSums(S)
  m <- mean(deg[deg > 0])
  if (is.na(m) || m == 0) return(S)
  S / m
}

#' pairTable
#'
#' Tumor-organoid pairs. An organoid pairs with the tumor whose `sample_id` equals its
#' `parent_sample_id`. If `parent_sample_id` is NA, it pairs with every tumor of its case.
#' @param samples Sample table (see [buildSampleGraph()]).
#' @return Data frame with columns `tumor` and `organoid` (profile ids).
#' @export
pairTable <- function(samples) {
  tum <- samples[samples$role %in% "tumor", c("profile_id", "case_id", "sample_id"), drop = FALSE]
  org <- samples[samples$role %in% "organoid", c("profile_id", "case_id", "parent_sample_id"), drop = FALSE]
  if (!nrow(tum) || !nrow(org)) return(data.frame(tumor = character(0), organoid = character(0)))
  e <- merge(org, tum, by = "case_id", suffixes = c("_org", "_tum"))
  e <- e[is.na(e$parent_sample_id) | e$parent_sample_id == e$sample_id, , drop = FALSE]
  data.frame(tumor = e$profile_id_tum, organoid = e$profile_id_org, stringsAsFactors = FALSE)
}

#' pairGraph
#'
#' Weight 1 between each organoid and its parent tumor(s) (see [pairTable()]).
#' @param samples Sample table.
#' @return Square matrix over `samples$profile_id`.
#' @export
pairGraph <- function(samples) {
  ids <- samples$profile_id
  S <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  pt <- pairTable(samples)
  if (!nrow(pt)) return(S)
  idx <- cbind(match(pt$organoid, ids), match(pt$tumor, ids))
  S[idx] <- 1
  S[idx[, 2:1, drop = FALSE]] <- 1
  S
}

#' annotGraph
#'
#' `1 - Gower distance` over patient annotations, sparsified to `k_nn` neighbors. Columns
#' that are constant or entirely missing are dropped; missing values are skipped by Gower.
#' Ties are broken by absolute age difference, then by profile id.
#' @param samples Sample table.
#' @param vars Annotation columns to use (those absent from `samples` are ignored).
#' @param k_nn Neighbors kept per sample.
#' @return Square matrix over `samples$profile_id`.
#' @export
annotGraph <- function(samples,
                       vars = c("cancer_type", "tumor_descriptor", "prior_treatment", "age", "sex", "ancestry"),
                       k_nn = 10) {
  ids <- samples$profile_id
  empty <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  df <- samples[, intersect(vars, names(samples)), drop = FALSE]
  df[] <- lapply(df, function(x) if (is.numeric(x)) x else factor(x))
  informative <- vapply(df, function(x) length(unique(stats::na.omit(x))) > 1, logical(1))
  df <- df[, informative, drop = FALSE]
  if (ncol(df) == 0) return(empty)
  S <- 1 - as.matrix(cluster::daisy(df, metric = "gower"))
  S[is.na(S)] <- 0
  dimnames(S) <- list(ids, ids)
  age_gap <- if ("age" %in% names(samples)) abs(outer(samples$age, samples$age, "-")) else NULL
  knnSparsify(S, k_nn, secondary = age_gap)
}

#' purityGraph
#'
#' Gaussian-kernel similarity of DNA-based purity between tumors only, sparsified to `k_nn`
#' neighbors. The bandwidth is the median distance from each tumor to its `k_nn`-th nearest
#' purity neighbor. Organoids and tumors with NA purity get no edges.
#' @param samples Sample table with `role` and `purity_dna`.
#' @param k_nn Neighbors kept per tumor.
#' @return Square matrix over `samples$profile_id`.
#' @export
purityGraph <- function(samples, k_nn = 10) {
  ids <- samples$profile_id
  S <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  t_idx <- which(samples$role %in% "tumor" & !is.na(samples$purity_dna))
  if (length(t_idx) < 2) return(S)
  p <- samples$purity_dna[t_idx]
  D <- abs(outer(p, p, "-"))
  kth <- apply(D, 1, function(r) sort(r)[min(k_nn + 1, length(r))])  # +1 skips the zero self-distance
  sigma <- stats::median(kth)
  if (sigma <= 0) sigma <- min(D[D > 0], 1)  # every purity tied
  K <- exp(-D^2 / (2 * sigma^2))
  dimnames(K) <- list(ids[t_idx], ids[t_idx])
  S[t_idx, t_idx] <- knnSparsify(K, k_nn)
  S
}

#' buildSampleGraph
#'
#' `S = w_pair * S_pair + w_annot * S_annot + w_pur * S_pur`, with each part first divided by
#' its mean degree (spec Section 6). Build it from training samples only.
#' @param samples Sample table with columns `profile_id, case_id, sample_id, parent_sample_id,
#'   role, cancer_type, tumor_descriptor, prior_treatment, age, sex, ancestry, purity_dna`.
#' @param weights Named weights for the parts `pair`, `annot`, `purity`.
#' @param k_nn Neighbors kept per sample in the annotation and purity parts.
#' @return A list with `S` (square matrix over `samples$profile_id`) and `summary`
#'   (columns part, n_edges, mean_degree, sparsity, weight).
#' @export
buildSampleGraph <- function(samples, weights = c(pair = 1, annot = 1, purity = 1), k_nn = 10) {
  parts <- list(pair = pairGraph(samples),
                annot = annotGraph(samples, k_nn = k_nn),
                purity = purityGraph(samples, k_nn = k_nn))
  scaled <- lapply(parts, scaleMeanDegree)
  S <- Reduce(`+`, Map(function(m, w) w * m, scaled, weights[names(scaled)]))
  diag(S) <- 0
  mean_deg <- function(m) {
    d <- rowSums(m)
    if (any(d > 0)) mean(d[d > 0]) else 0
  }
  summary <- data.frame(part = names(parts),
                        n_edges = vapply(parts, function(m) sum(m[upper.tri(m)] > 0), numeric(1)),
                        mean_degree = vapply(scaled, mean_deg, numeric(1)),
                        sparsity = vapply(parts, function(m) mean(m[upper.tri(m)] == 0), numeric(1)),
                        weight = unname(weights[names(parts)]),
                        row.names = NULL, stringsAsFactors = FALSE)
  list(S = S, summary = summary)
}
