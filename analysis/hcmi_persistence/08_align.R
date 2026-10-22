## Step 08: choose the pair-alignment setting and describe it.
##   Rscript analysis/hcmi_persistence/08_align.R [--k_clusters 8] [--n_perm 2000] [--n_boot 500] [--cores 8]
## Scores every setting in tasks/align_coess_2d.csv on held-out pairing (pairingScore() on each fold's
## held-out samples) and on the share of reference-fit programs enriched for a Reactome pathway, then
## picks the balanced setting (balanceSelect(); the baseline is plain NMF, alpha_1 = alpha_2 = 0).
## For the selected setting and the baseline it runs combineArm() (persistence map) and clusters all
## samples on their out-of-fold program activity.
source("analysis/hcmi_persistence/config.R")
library(ggplot2)
args <- parseArgs(list(k_clusters = "8", n_perm = "2000", n_boot = "500", n_perm_match = "1000",
                       cores = Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
k_clusters <- as.integer(args$k_clusters)
tasks <- as.data.frame(fread(file.path(cfg$analysis, "tasks", "align_coess_2d.csv")))
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
results <- lapply(split(tasks$task_id, tasks$setting_id), function(ids)
  lapply(file.path(cfg$run, "fits", "coess_2d", "align", paste0(ids, ".rds")), readResult))
genes <- colnames(Filter(Negate(is.null), unlist(results, recursive = FALSE))[[1]]$fit$H)
sets <- lapply(reactome, intersect, genes)
sets <- sets[lengths(sets) >= 10 & lengths(sets) <= 500]
## Reported, not used for selection: how many different pathways lead the enriched programs.
distinctTop <- function(H) {
  e <- enrichHypergeom(H, sets)
  e <- e[e$fdr < 0.05, ]
  e <- e[order(e$fdr), ]
  length(unique(e$gene_set[!duplicated(e$program)]))
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
  t1 <- tasks[tasks$setting_id == id, ][1, ]
  data.frame(setting_id = id, sample_graph = t1$sample_graph, alpha_1 = t1$alpha_1, alpha_2 = t1$alpha_2,
             n_fits = sum(done), all_converged = all(conv) && sum(done) == 6,
             pairing = if (length(held)) mean(ps$score) else NA_real_,
             pairing_se = if (length(held)) stats::sd(ps$score) / sqrt(nrow(ps)) else NA_real_,
             own_tumor_first = if (length(held)) mean(ps$score == 1) else NA_real_,
             enrichment = if (length(ref)) enrichedShare(ref[[1]]$H, sets) else NA_real_,
             distinct_pathways = if (length(ref)) distinctTop(ref[[1]]$H) else NA_integer_,
             recon_test = mean(vapply(held, function(f) f$recon_test, numeric(1))),
             n_iter_median = stats::median(vapply(fits[done], function(f) f$n_iter, numeric(1))))
}
tab <- do.call(rbind, lapply(names(results), score))
tab$baseline <- tab$alpha_1 == 0 & tab$alpha_2 == 0
tab$eligible <- tab$all_converged & !is.na(tab$pairing) & !is.na(tab$enrichment)
tab <- balanceSelect(tab)
fwrite(tab, file.path(out_dir, "settings.csv"))
writeProvenance(file.path(out_dir, "settings.csv"), inputs = reactome_file)
print(tab[order(-tab$pairing), c("setting_id", "all_converged", "pairing", "own_tumor_first", "enrichment",
                                 "distinct_pathways", "floor_ok", "pareto", "selected")], row.names = FALSE)
selected <- tab$setting_id[tab$selected]
base_id <- tab$setting_id[tab$baseline]
if (!length(selected)) message("no setting pairs organoids better than the baseline; describing the baseline only")

p <- ggplot(tab[tab$eligible, ], aes(pairing, enrichment, colour = sample_graph, shape = factor(alpha_2))) +
  geom_point(size = 3) +
  geom_point(data = tab[tab$selected | tab$baseline, ], shape = 21, size = 6, colour = "black") +
  ggrepel::geom_text_repel(data = tab[tab$selected | tab$baseline, ], aes(label = setting_id), colour = "black",
                           size = 3) +
  labs(x = "Held-out pairing score (1 = own tumor joins first; 0.5 = chance)",
       y = "Share of programs enriched for a Reactome pathway", colour = "Sample graph", shape = "alpha_2") +
  theme_minimal()
ggsave(file.path(out_dir, "figures", "tradeoff.png"), p, width = 7.5, height = 5.5, dpi = 200)

## 2. Persistence map and clusters for the selected setting and the baseline.
describe <- function(id) {
  expected <- unique(tasks[tasks$setting_id == id, c("setting_id", "alpha_1", "alpha_2")])
  out <- combineArm(Filter(Negate(is.null), results[[id]]), expected, samples, pairs, sets, mut, est,
                    n_perm = as.integer(args$n_perm), n_boot = as.integer(args$n_boot),
                    n_perm_match = as.integer(args$n_perm_match), seed = cfg$seed, cores = as.integer(args$cores))
  d <- file.path(out_dir, id)
  dir.create(d, showWarnings = FALSE)
  saveRDS(out, file.path(d, "combined.rds"))
  fwrite(out$programs, file.path(d, "programs.csv"))
  fwrite(out$enrichment, file.path(d, "enrichment.csv"))
  ## Program label: its most enriched Reactome pathway.
  e <- out$enrichment[order(out$enrichment$fdr), ]
  e <- e[!duplicated(e$program), ]
  stable <- out$programs$program[out$programs$stable]
  label <- stats::setNames(ifelse(e$fdr < 0.05, sub("^REACTOME_", "", e$gene_set), "(none)"), e$program)[stable]
  A <- out$oof[, stable, drop = FALSE]
  Z <- scale(log1p(A))
  Z[is.na(Z)] <- 0                                 # a program unmatched in a sample's fold sits at the mean
  tree <- stats::hclust(stats::dist(Z), method = "ward.D2")
  ids <- rownames(Z)
  pr <- pairs[pairs$organoid %in% ids & pairs$tumor %in% ids, ]
  cuts <- stats::cutree(tree, k = 2:40)
  co <- data.frame(setting_id = id, k = 2:40,
                   co_clustered = colMeans(cuts[pr$organoid, , drop = FALSE] == cuts[pr$tumor, , drop = FALSE]))
  cl <- cuts[, as.character(k_clusters)]
  role <- samples$role[match(ids, samples$profile_id)]
  means <- function(r) rowsum(Z[role == r, , drop = FALSE], cl[role == r]) / as.vector(table(cl[role == r]))
  zt <- means("tumor")
  zo <- means("organoid")
  both <- intersect(rownames(zt), rownames(zo))
  long <- data.frame(cluster = rep(both, times = ncol(Z)), program = rep(colnames(Z), each = length(both)),
                     tumor_z = as.vector(zt[both, , drop = FALSE]), organoid_z = as.vector(zo[both, , drop = FALSE]))
  long$reactome <- label[long$program]
  long$status <- ifelse(long$tumor_z > 0.5 & long$organoid_z > 0.5, "transfers",
                        ifelse(long$tumor_z > 0.5 & long$organoid_z < 0, "lost", ""))
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
                       labels = paste(colnames(Z), substr(label, 1, 40))[pord])
  ann <- data.frame(sample = factor(ids, levels = levels(hm$sample)), role = role, cluster = factor(cl))
  heat <- ggplot(hm, aes(sample, program, fill = pmax(pmin(z, 3), -3))) + geom_raster() +
    scale_fill_gradient2(low = "#3b6fb6", mid = "white", high = "#c0392b", name = "z") +
    theme_minimal(base_size = 7) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = sprintf("Samples (Ward tree on held-out program activity; %d clusters)", k_clusters), y = NULL,
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
ggsave(file.path(out_dir, "figures", "co_clustering_by_k.png"),
       ggplot(co, aes(k, co_clustered, colour = setting_id)) + geom_line() + geom_point(size = 1) +
         labs(x = "Number of sample clusters (tree cut)", y = "Share of pairs in the same cluster", colour = NULL) +
         theme_minimal(), width = 7, height = 4.5, dpi = 200)
groups <- do.call(rbind, lapply(desc, function(x) {
  p <- x$programs[x$programs$stable, ]
  data.frame(setting_id = x$id, as.list(table(factor(p$group, c("persists", "drops_out", "present_not_patient_specific",
                                                                "reduced_concordant")))))
}))
fwrite(groups, file.path(out_dir, "groups.csv"))
print(groups, row.names = FALSE)
message("08_align complete: ", out_dir)
