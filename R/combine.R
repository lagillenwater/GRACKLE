## Combine the cross-fitting results of one gene-graph arm (spec Sections 7.3, 8.2, 8.3, 9).

#' settingTable
#' @param results Task results: lists with meta, fit and error, as written by 05_fit.R.
#' @param expected Data frame with
#'   setting_id, alpha_1, alpha_2 for every setting that should exist.
#' @param n_folds Number of cross-validation folds (fold 0 is the reference fit).
#' @return One row per expected setting: complete, all_converged, recon_mean, recon_se, missing_or_failed.
#' @export
settingTable <- function(results, expected, n_folds = 5) {
  meta <- do.call(rbind, lapply(results, function(r) data.frame(
    setting_id = r$meta$setting_id, fold = r$meta$fold,
    ok = is.null(r$error) && !is.null(r$fit),
    converged = !is.null(r$fit) && isTRUE(r$fit$converged),
    recon_test = if (is.null(r$fit)) NA_real_ else r$fit$recon_test, stringsAsFactors = FALSE)))
  if (is.null(meta)) meta <- data.frame(setting_id = character(0), fold = integer(0), ok = logical(0),
                                        converged = logical(0), recon_test = numeric(0))
  do.call(rbind, lapply(seq_len(nrow(expected)), function(i) {
    m <- meta[meta$setting_id == expected$setting_id[i], , drop = FALSE]
    present <- m$fold[m$ok]
    complete <- all(0:n_folds %in% present)
    cv <- m$recon_test[m$ok & m$fold > 0]
    data.frame(setting_id = expected$setting_id[i], alpha_1 = expected$alpha_1[i], alpha_2 = expected$alpha_2[i],
               complete = complete, all_converged = complete && all(m$converged[m$ok]),
               recon_mean = if (complete) mean(cv) else NA_real_,
               recon_se = if (complete) stats::sd(cv) / sqrt(length(cv)) else NA_real_,
               missing_or_failed = paste(setdiff(0:n_folds, present), collapse = ","),
               stringsAsFactors = FALSE)
  }))
}

#' matchSetting
#' @param fits The successful results of one setting (reference fold 0 plus folds 1..5).
#' @param n_perm Gene-permutation draws per fold.
#' @param seed Random seed.
#' @return A list with tau, matches (fold, ref, fold_prog, cosine, matched), reproducibility, stability.
#' @export
matchSetting <- function(fits, n_perm = 1000, seed = 1) {
  ref <- Filter(function(r) r$meta$fold == 0, fits)[[1]]$fit$H
  folds <- Filter(function(r) r$meta$fold > 0, fits)
  null <- unlist(lapply(folds, function(r) matchNull(ref, r$fit$H, n_perm, seed + r$meta$fold)))
  tau <- unname(stats::quantile(null, 0.99))
  matches <- do.call(rbind, lapply(folds, function(r) cbind(fold = r$meta$fold, matchPrograms(ref, r$fit$H, tau))))
  list(tau = tau, matches = matches, reproducibility = mean(matches$cosine), stability = programStability(matches))
}

#' assembleOutOfFold
#' @param folds Results for folds 1..5 of one setting.
#' @param matches From [matchSetting()].
#' @return Held-out profiles x reference programs (`LV<ref>`), NA where a program is unmatched in that fold.
#' @export
assembleOutOfFold <- function(folds, matches) {
  do.call(rbind, lapply(folds, function(fr) {
    m <- matches[matches$fold == fr$meta$fold, , drop = FALSE]
    m <- m[order(m$ref), , drop = FALSE]
    A <- fr$fit$W_test[, m$fold_prog, drop = FALSE]
    A[, !m$matched] <- NA
    colnames(A) <- paste0("LV", m$ref)
    A
  }))
}

#' pairActivities
#' @param oof Out-of-fold activity.
#' @param pairs From [pairTable()].
#' @param samples Sample table.
#' @param program Column of `oof`.
#' @return Data frame with tumor_act, organoid_act, cancer_type, fold (complete rows only).
#' @export
pairActivities <- function(oof, pairs, samples, program) {
  p <- pairs[pairs$tumor %in% rownames(oof) & pairs$organoid %in% rownames(oof), , drop = FALSE]
  out <- data.frame(tumor_act = oof[p$tumor, program], organoid_act = oof[p$organoid, program],
                    cancer_type = samples$cancer_type[match(p$tumor, samples$profile_id)],
                    fold = samples$fold[match(p$tumor, samples$profile_id)], stringsAsFactors = FALSE)
  out[stats::complete.cases(out), , drop = FALSE]
}

#' scoreSetting
#' @param oof,pairs,samples See [pairActivities()].
#' @param stability From [matchSetting()].
#' @param n_perm,n_boot,seed Passed to the persistence tests.
#' @return One row per stable program with the persistence columns, conc_fdr, ret_fdr and group.
#' @export
scoreSetting <- function(oof, pairs, samples, stability, n_perm = 10000, n_boot = 2000, seed = 1) {
  stable <- paste0("LV", stability$ref[stability$stable])
  if (!length(stable)) return(data.frame())
  res <- do.call(rbind, lapply(stable, function(lv) {
    pa <- pairActivities(oof, pairs, samples, lv)
    cbind(program = lv, persistenceConcordance(pa, n_perm, n_boot, seed),
          persistenceRetention(pa, n_perm, n_boot, seed), stringsAsFactors = FALSE)
  }))
  res$conc_fdr <- stats::p.adjust(res$conc_p, "BH")
  res$ret_fdr <- stats::p.adjust(res$ret_p, "BH")
  res$group <- classifyPrograms(res$conc_fdr, res$ret_fdr, res$ret_median_log2)
  res
}

#' perTypeConcordance
#' @param oof,pairs,samples See [pairActivities()].
#' @param programs Programs to report.
#' @param min_pairs Minimum number of held-out pairs for a cancer type to be shown. Default 10.
#' @return Data frame with program, cancer_type, n_pairs, rho.
#' @export
perTypeConcordance <- function(oof, pairs, samples, programs, min_pairs = 10) {
  if (!length(programs)) return(data.frame())
  long <- do.call(rbind, lapply(programs, function(lv) cbind(program = lv, pairActivities(oof, pairs, samples, lv))))
  groups <- split(long, list(long$program, long$cancer_type), drop = TRUE)
  res <- do.call(rbind, lapply(groups, function(s) data.frame(
    program = s$program[1], cancer_type = s$cancer_type[1], n_pairs = nrow(s),
    rho = concordanceStat(normRankWithin(s$tumor_act, s$cancer_type), normRankWithin(s$organoid_act, s$cancer_type)),
    stringsAsFactors = FALSE)))
  res <- res[res$n_pairs >= min_pairs, , drop = FALSE]
  rownames(res) <- NULL
  res
}

#' annotatePrograms
#' @param H_ref Reference H.
#' @param oof Out-of-fold activity.
#' @param samples Sample table.
#' @param gene_sets Named list of gene sets (empty list to skip enrichment).
#' @param mut Tumor x gene mutation matrix, or NULL.
#' @param estimate Output of [estimateScores()], or NULL.
#' @return A list with programs (program, top_enrichment, purity_rho, purity_fdr, microenv,
#'   and stromal_rho, immune_rho, top_genotype when available), enrichment and genotype.
#' @export
annotatePrograms <- function(H_ref, oof, samples, gene_sets, mut = NULL, estimate = NULL, alpha = 0.05) {
  programs <- rownames(H_ref)
  out <- data.frame(program = programs, top_enrichment = NA_character_, stringsAsFactors = FALSE)
  enrich <- if (length(gene_sets)) enrichHypergeom(H_ref, gene_sets) else NULL
  if (!is.null(enrich)) {
    sig <- enrich[enrich$fdr < alpha, , drop = FALSE]
    sig <- sig[order(sig$program, sig$p), , drop = FALSE]
    top <- sig[!duplicated(sig$program), , drop = FALSE]
    out$top_enrichment <- top$gene_set[match(programs, top$program)]
  }
  tum <- samples$profile_id[samples$role == "tumor" & samples$profile_id %in% rownames(oof)]
  A <- oof[tum, programs, drop = FALSE]
  pur <- programAssociation(A, samples$purity_dna[match(tum, samples$profile_id)], "purity_dna")
  out$purity_rho <- pur$rho
  out$purity_fdr <- pur$fdr
  out$microenv <- pur$fdr < alpha & pur$rho < 0
  if (!is.null(estimate)) {
    out$stromal_rho <- programAssociation(A, estimate$stromal[match(tum, estimate$profile_id)], "stromal")$rho
    out$immune_rho <- programAssociation(A, estimate$immune[match(tum, estimate$profile_id)], "immune")$rho
  }
  geno <- NULL
  if (!is.null(mut)) {
    rows <- intersect(tum, rownames(mut))
    geno <- genotypeAssociation(A[rows, , drop = FALSE], mut[rows, , drop = FALSE])
    gs <- geno[geno$fdr < alpha, , drop = FALSE]
    gs <- gs[order(gs$program, gs$p), , drop = FALSE]
    gt <- gs[!duplicated(gs$program), , drop = FALSE]
    i <- match(programs, gt$program)
    out$top_genotype <- ifelse(is.na(i), NA_character_, paste0(gt$gene[i], ifelse(gt$delta[i] > 0, "(+)", "(-)")))
  }
  list(programs = out, enrichment = enrich, genotype = geno)
}

#' combineArm
#'
#' Selection, matching, persistence scoring and annotation for one gene-graph arm. Settings
#' with a missing, failed or non-converged fold are listed in `selection` and never scored.
#' @param results Task results for the arm.
#' @param expected Expected settings.
#' @param samples Sample table (included samples, with fold).
#' @param pairs From [pairTable()].
#' @param gene_sets,mut,estimate Annotation inputs (see [annotatePrograms()]).
#' @param n_perm,n_boot Persistence draws.
#' @param n_perm_match Matching-null draws.
#' @param seed Random seed.
#' @param cores Cores used to score settings in parallel.
#' @param n_folds Number of folds.
#' @return A list with selection, selected, tau_selected, programs, maps_long, matches,
#'   per_type, enrichment, genotype, ref_H, oof, convergence.
#' @export
combineArm <- function(results, expected, samples, pairs, gene_sets = list(), mut = NULL, estimate = NULL,
                       n_perm = 10000, n_boot = 2000, n_perm_match = 1000, seed = 20261015, cores = 1,
                       n_folds = 5) {
  ok <- Filter(function(r) is.null(r$error) && !is.null(r$fit), results)
  tab <- settingTable(results, expected, n_folds)
  by_setting <- split(ok, vapply(ok, function(r) r$meta$setting_id, character(1)))
  ids <- tab$setting_id[tab$complete]
  matched <- stats::setNames(lapply(ids, function(id) matchSetting(by_setting[[id]], n_perm_match, seed)), ids)
  tab$reproducibility <- NA_real_
  tab$tau <- NA_real_
  tab$reproducibility[match(ids, tab$setting_id)] <- vapply(matched, function(m) m$reproducibility, numeric(1))
  tab$tau[match(ids, tab$setting_id)] <- vapply(matched, function(m) m$tau, numeric(1))
  selected <- selectSetting(tab)
  tab$selected <- tab$setting_id == selected
  scored <- parallel::mclapply(ids, function(id) {
    folds <- Filter(function(r) r$meta$fold > 0, by_setting[[id]])
    oof <- assembleOutOfFold(folds, matched[[id]]$matches)
    list(oof = oof, scores = scoreSetting(oof, pairs, samples, matched[[id]]$stability, n_perm, n_boot, seed))
  }, mc.cores = cores)
  bad <- vapply(scored, inherits, logical(1), what = "try-error")
  if (any(bad)) stop("combineArm: scoring failed for ", paste(ids[bad], collapse = ", "), ": ", as.character(scored[[which(bad)[1]]]))
  names(scored) <- ids
  ref_H <- Filter(function(r) r$meta$fold == 0, by_setting[[selected]])[[1]]$fit$H
  oof <- scored[[selected]]$oof
  ann <- annotatePrograms(ref_H, oof, samples, gene_sets, mut, estimate)
  st <- matched[[selected]]$stability
  st$program <- paste0("LV", st$ref)
  programs <- merge(st, scored[[selected]]$scores, by = "program", all.x = TRUE)
  programs <- merge(programs, ann$programs, by = "program", all.x = TRUE)
  programs <- programs[order(programs$ref), , drop = FALSE]
  maps_long <- do.call(rbind, lapply(ids, function(id) {
    s <- scored[[id]]$scores
    if (nrow(s)) cbind(setting_id = id, s, stringsAsFactors = FALSE) else NULL
  }))
  field <- function(r, f) if (is.null(r$fit)) NA else r$fit[[f]]
  convergence <- do.call(rbind, lapply(results, function(r) data.frame(
    task_id = if (is.null(r$meta$task_id)) NA_character_ else r$meta$task_id,
    setting_id = r$meta$setting_id, fold = r$meta$fold,
    n_iter = field(r, "n_iter"), converged = field(r, "converged"),
    final_total = if (is.null(r$fit)) NA_real_ else r$fit$final_loss$total,
    error = if (is.null(r$error)) NA_character_ else r$error, stringsAsFactors = FALSE)))
  list(selection = tab, selected = selected, tau_selected = matched[[selected]]$tau, programs = programs,
       maps_long = maps_long, matches = matched[[selected]]$matches,
       per_type = perTypeConcordance(oof, pairs, samples, st$program[st$stable]),
       enrichment = ann$enrichment, genotype = ann$genotype, ref_H = ref_H, oof = oof, convergence = convergence)
}
