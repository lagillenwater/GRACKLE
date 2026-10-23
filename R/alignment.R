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

#' standardizeActivity
#'
#' Program activity on a common scale: log1p-transformed and each program standardized over all samples,
#' or separately within each role, which removes the shift every organoid shares (such as its lost
#' microenvironment programs). Missing values (a program unmatched in a sample's fold) and programs
#' with no variance read as 0.
#' @param A Samples x programs activity (out-of-fold W).
#' @param role Optional role for each row of `A`; when given, programs are standardized within each role.
#' @return Matrix shaped like `A`, with its dimnames.
#' @export
standardizeActivity <- function(A, role = NULL) {
  L <- log1p(A)
  idx <- if (is.null(role)) list(seq_len(nrow(L))) else split(seq_len(nrow(L)), role)
  Z <- do.call(rbind, lapply(idx, function(i) scale(L[i, , drop = FALSE])))[order(unlist(idx)), , drop = FALSE]
  Z <- matrix(Z, nrow(L), dimnames = dimnames(A))
  Z[is.na(Z)] <- 0
  Z
}

#' pairedTree
#'
#' A sample tree that keeps given clusters whole and joins each tumor-organoid pair within a cluster
#' first, so every pair is drawn as adjacent leaves (tumor first). Joins are Ward linkage at every
#' level: the samples of each pair (or of a tumor or organoid linked to two others), then those units
#' within each cluster, then the clusters, each level from the centroids and sizes of the groups below
#' it. A sample whose partner sits in another cluster is a unit of its own. A join is raised to its
#' children's height where needed so the tree draws without crossing lines.
#' @param X Samples x features, with sample ids as rownames.
#' @param cluster Cluster of each row of `X`.
#' @param pairs Data frame with tumor and organoid ids.
#' @return An `hclust` tree over the rows of `X`; its joins run pairs first, then within each cluster,
#'   then between clusters, so cutting it into as many groups as there are clusters returns them.
#' @export
pairedTree <- function(X, cluster, pairs) {
  ids <- rownames(X)
  n <- nrow(X)
  cluster <- as.vector(cluster)
  pr <- pairs[pairs$tumor %in% ids & pairs$organoid %in% ids, , drop = FALSE]
  pr <- pr[cluster[match(pr$tumor, ids)] == cluster[match(pr$organoid, ids)], , drop = FALSE]
  unit <- seq_len(n)
  if (nrow(pr)) {
    comp <- igraph::components(igraph::graph_from_data_frame(pr[c("tumor", "organoid")], directed = FALSE))$membership
    unit[match(names(comp), ids)] <- n + comp
  }
  merge <- matrix(0L, 0, 2)
  height <- numeric(0)
  ## Joins groups (nodes: -row for a sample, a join's index otherwise) by Ward linkage, starting from the
  ## Ward distance between their centroids, sqrt(2 n_a n_b / (n_a + n_b)) times the Euclidean distance.
  join <- function(nodes, rows) {
    if (length(nodes) == 1) return(nodes)
    size <- lengths(rows)
    centroid <- do.call(rbind, lapply(rows, function(r) colMeans(X[r, , drop = FALSE])))
    W <- as.matrix(stats::dist(centroid)) * sqrt(2 * outer(size, size) / outer(size, size, "+"))
    h <- stats::hclust(stats::as.dist(W), method = "ward.D2", members = size)
    m <- h$merge
    m[h$merge < 0] <- nodes[-h$merge[h$merge < 0]]
    m[h$merge > 0] <- h$merge[h$merge > 0] + nrow(merge)
    merge <<- rbind(merge, m)
    height <<- c(height, h$height)
    nrow(merge)
  }
  o <- order(unit, !(ids %in% pairs$tumor))   # tumors first within a pair
  units <- split(o, unit[o])
  unit_node <- vapply(units, function(r) join(-r, as.list(r)), numeric(1))
  by_cluster <- split(seq_along(units), cluster[vapply(units, `[`, integer(1), 1)])
  cluster_rows <- lapply(by_cluster, function(u) unlist(units[u], use.names = FALSE))
  cluster_node <- vapply(by_cluster, function(u) join(unit_node[u], units[u]), numeric(1))
  join(cluster_node, cluster_rows)
  for (i in seq_len(nrow(merge))) height[i] <- max(height[i], height[merge[i, merge[i, ] > 0]])
  leaves <- function(i) if (i < 0) -i else c(leaves(merge[i, 1]), leaves(merge[i, 2]))
  storage.mode(merge) <- "integer"
  structure(list(merge = merge, height = height, order = leaves(nrow(merge)), labels = ids, method = "ward.D2",
                 call = match.call(), dist.method = "euclidean"), class = "hclust")
}

#' clusterActivity
#'
#' Mean program activity of the tumors and of the organoids in each sample cluster. Activity is
#' log1p-transformed and each program standardized over all samples, so a program every organoid loses
#' reads low in each cluster's organoids. A program transfers in a cluster when both means exceed `high`,
#' and is lost when the tumors' mean exceeds `high` and the organoids' is below zero.
#' @param A Samples x programs activity (out-of-fold W).
#' @param role "tumor" or "organoid" for each row of `A`.
#' @param cluster Sample cluster for each row of `A`.
#' @param high Threshold on the standardized mean.
#' @return Data frame with cluster, program, tumor_z, organoid_z and status ("transfers", "lost" or ""),
#'   for the clusters that hold both tumors and organoids.
#' @export
clusterActivity <- function(A, role, cluster, high = 0.5) {
  Z <- standardizeActivity(A)
  means <- function(r) {
    s <- rowsum(Z[role == r, , drop = FALSE], cluster[role == r])
    s / as.vector(table(cluster[role == r])[rownames(s)])
  }
  zt <- means("tumor")
  zo <- means("organoid")
  both <- intersect(rownames(zt), rownames(zo))
  out <- data.frame(cluster = rep(both, times = ncol(Z)), program = rep(colnames(Z), each = length(both)),
                    tumor_z = as.vector(zt[both, , drop = FALSE]), organoid_z = as.vector(zo[both, , drop = FALSE]),
                    stringsAsFactors = FALSE)
  out$status <- ifelse(out$tumor_z > high & out$organoid_z > high, "transfers",
                       ifelse(out$tumor_z > high & out$organoid_z < 0, "lost", ""))
  out
}

#' pathwayClusters
#'
#' Groups programs by how they rank gene sets and names each group. Each program ranks every set by its
#' GSEA NES; programs are joined by Ward linkage on the Euclidean distance between their standardized
#' rankings (a monotone function of 1 - Spearman correlation), and the tree is cut into `k` groups, or
#' where the mean silhouette is highest when `k` is NULL. A group is named by the set it ranks highest
#' relative to the other programs, among the sets significant in at least half its programs (or, when
#' none is, in the most of them), so a set significant in every program does not win by that alone.
#' @param g Output of [gseaPrograms()].
#' @param programs Programs to group (at least three).
#' @param fdr Significance threshold on `g$fdr`.
#' @param k Number of groups; NULL chooses it by silhouette.
#' @param k_max Largest number of groups considered when `k` is NULL.
#' @return Data frame with program, cluster (numbered along the tree), label (a gene set, NA when no
#'   program in the group has a significant set), n_sig (programs in the group where the label is
#'   significant) and n_programs, one row per program in tree order.
#' @export
pathwayClusters <- function(g, programs, fdr = 0.05, k = NULL, k_max = 8) {
  g <- g[g$program %in% programs, , drop = FALSE]
  sets <- unique(g$gene_set)
  at <- cbind(match(g$gene_set, sets), match(g$program, programs))
  N <- matrix(0, length(sets), length(programs), dimnames = list(sets, programs))
  N[at] <- g$nes
  S <- matrix(FALSE, length(sets), length(programs), dimnames = list(sets, programs))
  S[at] <- g$fdr < fdr
  R <- apply(N, 2, rank) / nrow(N)
  D <- stats::dist(t(scale(R)))
  tree <- stats::hclust(D, method = "ward.D2")
  if (is.null(k)) {
    ks <- 2:min(k_max, length(programs) - 1)
    sil <- vapply(ks, function(j) mean(cluster::silhouette(stats::cutree(tree, j), D)[, 3]), numeric(1))
    k <- ks[which.max(sil)]
  }
  cl <- stats::cutree(tree, k)
  cl <- match(cl, unique(cl[tree$order]))
  named <- vapply(seq_len(max(cl)), function(k) {
    inn <- cl == k
    d <- rowMeans(R[, inn, drop = FALSE]) - rowMeans(R[, !inn, drop = FALSE])
    n_sig <- rowSums(S[, inn, drop = FALSE])
    ok <- if (any(n_sig >= sum(inn) / 2)) n_sig >= sum(inn) / 2 else n_sig == max(n_sig) & n_sig > 0
    best <- names(d)[ok][which.max(d[ok])]
    if (length(best)) c(best, n_sig[[best]]) else c(NA_character_, "0")
  }, character(2))
  o <- tree$order
  data.frame(program = programs[o], cluster = cl[o], label = named[1, cl[o]], n_sig = as.integer(named[2, cl[o]]),
             n_programs = as.vector(table(cl))[cl[o]], stringsAsFactors = FALSE)
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
