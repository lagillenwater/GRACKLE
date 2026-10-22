## Step 08: choose the pair-alignment setting and describe it.
##   Rscript analysis/hcmi_persistence/08_align.R [--k_clusters 8] [--n_perm 2000] [--n_boot 500] [--cores 8]
## Scores every setting in every tasks/align_<arm>.csv (a setting is named <arm>_<setting_id>) on held-out
## pairing (pairingScore() on each fold's held-out samples) and on the share of reference-fit programs
## and on pathway enrichment of the reference-fit programs: enrichmentStrength(), the median best GSEA NES
## (the GSEA and top-5% enriched shares are reported beside it, as they saturate). balanceSelect() keeps
## settings that pair better than plain NMF (alpha_1 = alpha_2 = 0) and ranks them by distance to the ideal;
## walking that order, the first setting whose persistence map passes the positive control (controlHolds(),
## from a combineArm() with --n_perm_guard permutations) is selected. An arm whose fits are absent leaves
## its settings ineligible.
## For the selected setting and the baseline it runs combineArm() (persistence map) and clusters all
## samples on their out-of-fold program activity.
source("analysis/hcmi_persistence/config.R")
library(ggplot2)
args <- parseArgs(list(k_clusters = "8", n_perm = "2000", n_boot = "500", n_perm_match = "1000", n_perm_guard = "500",
                       cores = Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
k_clusters <- as.integer(args$k_clusters)
tasks <- do.call(rbind, lapply(list.files(file.path(cfg$analysis, "tasks"), pattern = "^align_.*\\.csv$", full.names = TRUE),
                                function(f) as.data.frame(fread(f))))
tasks$key <- paste(tasks$arm, tasks$setting_id, sep = "_")
out_dir <- file.path(cfg$run, "results", "align")
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
samples <- readRDS(file.path(cfg$processed, "samples_included.rds"))
pairs <- pairTable(samples)
optional <- function(path, reader) if (file.exists(path)) reader(path) else NULL
## Genotype association is not part of this step (genotypeAssociation() fails when a program has no
## non-missing activity among a gene's mutated samples).
mut <- NULL
est <- optional(file.path(cfg$processed, "estimate.csv"), function(p) as.data.frame(fread(p)))

## Reactome (MSigDB C2:CP:REACTOME through msigdbr), cached with its version.
reactome_file <- file.path(cfg$raw, "annotation", sprintf("reactome_msigdbr_%s.rds", packageVersion("msigdbr")))
if (!file.exists(reactome_file)) {
  r <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
  dir.create(dirname(reactome_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(split(r$gene_symbol, r$gs_name), reactome_file)
}
reactome <- readRDS(reactome_file)

## 1. Score every setting.
results <- lapply(split(seq_len(nrow(tasks)), tasks$key), function(i)
  lapply(file.path(cfg$run, "fits", tasks$arm[i], "align", paste0(tasks$task_id[i], ".rds")), readResult))
genes <- colnames(Filter(Negate(is.null), unlist(results, recursive = FALSE))[[1]]$fit$H)
sets <- lapply(reactome, intersect, genes)
sets <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]
cores <- as.integer(args$cores)
## A program's leading pathway: its most significant one by GSEA (NA when none passes FDR < 0.05).
leading <- function(g, programs) {
  g <- g[order(g$fdr), ]
  g <- g[!duplicated(g$program) & g$fdr < 0.05, ]
  stats::setNames(g$gene_set[match(programs, g$program)], programs)
}
score <- function(id) {
  r <- results[[id]]
  fits <- lapply(r, function(x) if (is.null(x) || !is.null(x$error)) NULL else x$fit)
  folds <- vapply(r, function(x) if (is.null(x)) NA_integer_ else as.integer(x$meta$fold), integer(1))
  done <- !vapply(fits, is.null, logical(1))
  conv <- vapply(fits, function(f) isTRUE(f$converged), logical(1))
  held <- fits[done & folds > 0]
  ps <- do.call(rbind, lapply(held, function(f) pairingScore(`rownames<-`(f$W_test, f$test_ids), samples)))
  ref <- fits[done & folds == 0]
  g <- if (length(ref)) gseaPrograms(ref[[1]]$H, sets, cores = cores) else NULL
  lead <- if (length(ref)) leading(g, rownames(ref[[1]]$H)) else NA
  t1 <- tasks[tasks$key == id, ][1, ]
  data.frame(setting_id = id, arm = t1$arm, sample_graph = t1$sample_graph, alpha_1 = t1$alpha_1, alpha_2 = t1$alpha_2,
             n_fits = sum(done), all_converged = all(conv) && sum(done) == 6,
             pairing = if (length(held)) mean(ps$score) else NA_real_,
             pairing_se = if (length(held)) stats::sd(ps$score) / sqrt(nrow(ps)) else NA_real_,
             own_tumor_first = if (length(held)) mean(ps$score == 1) else NA_real_,
             enrichment = if (length(ref)) enrichmentStrength(g) else NA_real_,
             enriched_share = if (length(ref)) mean(!is.na(lead)) else NA_real_,       # enrichedShare(method = "gsea")
             distinct_pathways = if (length(ref)) length(unique(stats::na.omit(lead))) else NA_integer_,
             enrichment_ora = if (length(ref)) enrichedShare(ref[[1]]$H, sets, method = "ora") else NA_real_,
             recon_test = mean(vapply(held, function(f) f$recon_test, numeric(1))),
             n_iter_median = stats::median(vapply(fits[done], function(f) f$n_iter, numeric(1))))
}
tab <- do.call(rbind, lapply(names(results), score))
tab$baseline <- tab$alpha_1 == 0 & tab$alpha_2 == 0
tab$eligible <- tab$all_converged & !is.na(tab$pairing) & !is.na(tab$enrichment)
tab <- balanceSelect(tab)
## Positive control, walked in order of distance to the ideal; the first setting that passes is selected
## (it is also the closest of the passing settings that no other passing setting beats on both measures).
combineSetting <- function(id, n_perm, n_boot) {
  combineArm(Filter(Negate(is.null), results[[id]]), unique(tasks[tasks$key == id, c("setting_id", "alpha_1", "alpha_2")]),
             samples, pairs, sets, mut, est, n_perm = n_perm, n_boot = n_boot,
             n_perm_match = as.integer(args$n_perm_match), seed = cfg$seed, cores = cores)
}
tab$control_holds <- NA
tab$selected <- FALSE
for (id in tab$setting_id[tab$floor_ok][order(tab$distance[tab$floor_ok])]) {
  tab$control_holds[tab$setting_id == id] <- controlHolds(combineSetting(id, as.integer(args$n_perm_guard), 200)$programs)
  message(id, ": positive control ", if (tab$control_holds[tab$setting_id == id]) "holds" else "fails")
  if (tab$control_holds[tab$setting_id == id]) { tab$selected[tab$setting_id == id] <- TRUE; break }
}
fwrite(tab, file.path(out_dir, "settings.csv"))
writeProvenance(file.path(out_dir, "settings.csv"), inputs = reactome_file)
print(tab[order(-tab$pairing), c("setting_id", "all_converged", "pairing", "own_tumor_first", "enrichment",
                                 "enriched_share", "distinct_pathways", "enrichment_ora", "floor_ok", "control_holds", "selected")], row.names = FALSE)
selected <- tab$setting_id[tab$selected]
base_id <- tab$setting_id[tab$baseline]
if (!length(selected)) message("no setting pairs organoids better than the baseline; describing the baseline only")

p <- ggplot(tab[tab$eligible, ], aes(pairing, enrichment, colour = sample_graph, shape = factor(alpha_2))) +
  geom_point(size = 3) +
  geom_point(data = tab[tab$selected | tab$baseline, ], shape = 21, size = 6, colour = "black") +
  ggrepel::geom_text_repel(data = tab[tab$selected | tab$baseline, ], aes(label = setting_id), colour = "black",
                           size = 3) +
  labs(x = "Held-out pairing score (1 = own tumor joins first; 0.5 = chance)",
       y = "Reactome enrichment (median over programs of the best GSEA NES)", colour = "Sample graph", shape = "alpha_2") +
  facet_wrap(~arm) + theme_minimal()
ggsave(file.path(out_dir, "figures", "tradeoff.png"), p, width = 11, height = 5.5, dpi = 200)

## 2. Persistence map and clusters for the selected setting and the baseline.
describe <- function(id) {
  out <- combineSetting(id, as.integer(args$n_perm), as.integer(args$n_boot))
  d <- file.path(out_dir, id)
  dir.create(d, showWarnings = FALSE)
  saveRDS(out, file.path(d, "combined.rds"))
  fwrite(out$programs, file.path(d, "programs.csv"))
  fwrite(out$enrichment, file.path(d, "enrichment.csv"))
  ## Program label: its leading Reactome pathway by GSEA.
  g <- gseaPrograms(out$ref_H, sets, cores = cores)
  fwrite(g, file.path(d, "gsea.csv"))
  stable <- out$programs$program[out$programs$stable]
  lead <- leading(g, rownames(out$ref_H))
  label <- stats::setNames(ifelse(is.na(lead), "(none)", sub("^REACTOME_", "", lead)), names(lead))
  ## Programs that transfer: concordant between a donor's tumor and organoid (persists, reduced_concordant).
  transferring <- out$programs$program[out$programs$stable & out$programs$group %in% c("persists", "reduced_concordant")]
  if (length(transferring) < 2) transferring <- stable
  ids <- rownames(out$oof)
  role <- samples$role[match(ids, samples$profile_id)]
  pr <- pairs[pairs$organoid %in% ids & pairs$tumor %in% ids, ]
  ## log1p activity, each program standardized over all samples or separately within tumors and within
  ## organoids (which removes the shift every organoid shares, such as lost microenvironment programs).
  standardize <- function(keep, within_role) {
    A <- log1p(out$oof[, keep, drop = FALSE])
    Z <- if (within_role) do.call(rbind, lapply(split(seq_len(nrow(A)), role), function(i) scale(A[i, , drop = FALSE])))[ids, , drop = FALSE]
         else scale(A)
    Z[is.na(Z)] <- 0                               # unmatched in a sample's fold, or no variance within a role
    Z
  }
  variants <- data.frame(programs = c("transferring", "transferring", "all", "all"), within_role = c(TRUE, FALSE, TRUE, FALSE))
  co <- do.call(rbind, lapply(seq_len(nrow(variants)), function(v) {
    keep <- if (variants$programs[v] == "transferring") transferring else stable
    cuts <- stats::cutree(stats::hclust(stats::dist(standardize(keep, variants$within_role[v])), "ward.D2"), k = 2:40)
    data.frame(setting_id = id, programs = variants$programs[v], within_role = variants$within_role[v], k = 2:40,
               n_programs = length(keep), pairs_together = colSums(cuts[pr$organoid, , drop = FALSE] == cuts[pr$tumor, , drop = FALSE]),
               co_clustered = colMeans(cuts[pr$organoid, , drop = FALSE] == cuts[pr$tumor, , drop = FALSE]))
  }))
  ## Primary view: transferring programs, standardized within tumors and within organoids.
  Z <- standardize(transferring, TRUE)
  tree <- stats::hclust(stats::dist(Z), method = "ward.D2")
  cl <- stats::cutree(tree, k = k_clusters)
  ## Per cluster, every stable program standardized over all samples, so programs lost in organoids show.
  long <- clusterActivity(out$oof[, stable, drop = FALSE], role, cl)
  long$group <- out$programs$group[match(long$program, out$programs$program)]
  long$reactome <- label[long$program]
  long <- long[, c("cluster", "program", "tumor_z", "organoid_z", "group", "reactome", "status")]
  pc <- data.frame(cluster = names(table(cl)), n_tumors = as.vector(table(factor(cl[role == "tumor"], names(table(cl))))),
                   n_organoids = as.vector(table(factor(cl[role == "organoid"], names(table(cl))))),
                   pairs_together = as.vector(table(factor(cl[pr$organoid][cl[pr$organoid] == cl[pr$tumor]],
                                                           names(table(cl))))))
  fwrite(long[long$status != "", ], file.path(d, sprintf("cluster_programs_k%d.csv", k_clusters)))
  fwrite(pc, file.path(d, sprintf("clusters_k%d.csv", k_clusters)))
  fwrite(data.frame(sample = ids, role = role, cluster = cl), file.path(d, sprintf("sample_clusters_k%d.csv", k_clusters)))
  writeProvenance(file.path(d, "programs.csv"), inputs = reactome_file,
                  extra = list(setting = id, n_perm = args$n_perm, n_boot = args$n_boot, k_clusters = k_clusters))
  hm <- data.frame(sample = rep(ids, times = ncol(Z)), program = rep(colnames(Z), each = nrow(Z)), z = as.vector(Z))
  hm$sample <- factor(hm$sample, levels = ids[tree$order])
  pord <- stats::hclust(stats::dist(t(Z)), method = "ward.D2")$order
  hm$program <- factor(hm$program, levels = colnames(Z)[pord],
                       labels = paste(colnames(Z), substr(label[colnames(Z)], 1, 40))[pord])
  ann <- data.frame(sample = factor(ids, levels = levels(hm$sample)), role = role, cluster = factor(cl))
  heat <- ggplot(hm, aes(sample, program, fill = pmax(pmin(z, 3), -3))) + geom_raster() +
    scale_fill_gradient2(low = "#3b6fb6", mid = "white", high = "#c0392b", name = "z") +
    theme_minimal(base_size = 7) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = sprintf("Samples (Ward tree on held-out activity of transferring programs, standardized within tumors and organoids; %d clusters)", k_clusters), y = NULL,
         title = id)
  strip <- function(fill, values = NULL) {
    g <- ggplot(ann, aes(sample, 1, fill = .data[[fill]])) + geom_raster() + theme_void() +
      theme(legend.position = "right", legend.key.size = grid::unit(8, "pt"))
    if (is.null(values)) g else g + scale_fill_manual(values = values)
  }
  ggsave(file.path(d, "heatmap.png"),
         patchwork::wrap_plots(strip("role", c(tumor = "#555555", organoid = "#e69f00")), strip("cluster"), heat,
                               ncol = 1, heights = c(0.04, 0.04, 1)),
         width = 12, height = 8, dpi = 200)
  list(id = id, co = co, programs = out$programs)
}
desc <- lapply(unique(c(selected, base_id)), describe)
co <- do.call(rbind, lapply(desc, `[[`, "co"))
fwrite(co, file.path(out_dir, "co_clustering_by_k.csv"))
co$view <- paste0(co$programs, " programs, ", ifelse(co$within_role, "standardized within tumors / organoids", "standardized over all samples"))
ggsave(file.path(out_dir, "figures", "co_clustering_by_k.png"),
       ggplot(co, aes(k, co_clustered, colour = setting_id)) + geom_line() + geom_point(size = 0.8) + facet_wrap(~view) +
         labs(x = "Number of sample clusters (tree cut)", y = "Share of pairs in the same cluster", colour = NULL) +
         theme_minimal() + theme(legend.position = "bottom"), width = 10, height = 7, dpi = 200)
groups <- do.call(rbind, lapply(desc, function(x) {
  p <- x$programs[x$programs$stable, ]
  data.frame(setting_id = x$id, as.list(table(factor(p$group, c("persists", "drops_out", "present_not_patient_specific",
                                                                "reduced_concordant")))))
}))
fwrite(groups, file.path(out_dir, "groups.csv"))
print(groups, row.names = FALSE)
message("08_align complete: ", out_dir)
