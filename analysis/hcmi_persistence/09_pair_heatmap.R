## Step 09: paired tumors and organoids by sample cluster, for the selected pair-alignment setting and the
## baseline.
##   Rscript analysis/hcmi_persistence/09_pair_heatmap.R [--k_clusters 8] [--k_programs 8]
## Reads step 08's outputs under results/align (out-of-fold activity in combined.rds, the sample clusters,
## the persistence map and the GSEA table); refits nothing. The tree and clusters are step 08's: a Ward tree
## on the transferring programs standardized within tumors and within organoids (standardizeActivity()),
## cut into k clusters. Rows are the tumors and organoids in a pair, under a tree that keeps those clusters
## and joins each pair within a cluster first (pairedTree()), so a pair the clusters hold together sits on
## adjacent rows, tumor above organoid; a gap separates clusters, each labeled with the pairs it holds together.
## Columns are the stable programs that persist, cut into k_programs groups by how they rank Reactome
## pathways (pathwayClusters() on the preranked GSEA), each named by the pathway that sets it apart; within a
## group, programs follow the cluster where they peak.
source("analysis/hcmi_persistence/config.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
args <- parseArgs(list(k_clusters = "8", k_programs = "8"))
k_clusters <- as.integer(args$k_clusters)
k_programs <- as.integer(args$k_programs)
align_dir <- file.path(cfg$run, "results", "align")
fig_dir <- file.path(align_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
samples <- readRDS(file.path(cfg$processed, "samples_included.rds"))
pairs <- pairTable(samples)
tab <- as.data.frame(fread(file.path(align_dir, "settings.csv")))
ids <- unique(c(tab$setting_id[tab$selected], tab$setting_id[tab$baseline]))

## Colors from the dataviz reference palette: diverging blue-red with its neutral midpoint for activity
## (fixed limits so the two figures compare directly), neutral inks for role.
ink <- c(primary = "#0b0b0b", secondary = "#52514e", baseline = "#c3c2b7", surface = "#fcfcfb")
diverging <- c(low = "#2a78d6", mid = "#f0efec", high = "#e34948")
role_cols <- c(Tumor = ink[["secondary"]], Organoid = ink[["baseline"]])
gap <- 3   # empty rows between clusters
theme_set(theme_minimal(base_size = 8) +
            theme(panel.grid = element_blank(), text = element_text(colour = ink[["primary"]]),
                  axis.text = element_text(colour = ink[["secondary"]]),
                  strip.text = element_text(colour = ink[["primary"]], size = 7),
                  plot.background = element_rect(fill = ink[["surface"]], colour = NA)))

## Reactome names for display: prefix dropped, sentence case with acronyms and gene symbols (any word with a
## digit) kept, long names cut.
acronyms <- c("ADME", "ECM", "HCMV", "RNA", "DNA", "SLC", "GPI", "PCP", "PD", "II", "III", "SPMS", "HFTC", "AR", "NMD",
              "TCA", "ER", "MET", "O", "TCR", "CENPA", "HATS", "HDACS", "RRNA", "CSB", "ATR", "TCF", "WNT", "AP", "K")
## Plural acronyms read with a lowercase s.
plurals <- c(HATS = "HATs", HDACS = "HDACs", SPMS = "SPMs")
pathwayLabel <- function(x, n = 40) {
  w <- strsplit(sub("^REACTOME_", "", x), "_", fixed = TRUE)
  s <- vapply(w, function(t) {
    t <- ifelse(t %in% acronyms | grepl("[0-9]", t), t, tolower(t))
    t <- ifelse(t %in% names(plurals), plurals[t], t)
    t <- sub("^RRNA$", "rRNA", t)
    t[1] <- paste0(toupper(substr(t[1], 1, 1)), substring(t[1], 2))
    paste(t, collapse = " ")
  }, character(1))
  ifelse(nchar(s) > n, paste0(substr(s, 1, n - 1), "…"), s)
}

## Tree segments with leaves at x = 0 and leaf i (in label order) at y = leaf_y[i]. Merges are processed
## in order, so each child's position is known before its parent's.
treeSegments <- function(tree, leaf_y) {
  m <- tree$merge
  pos <- numeric(nrow(m))
  segs <- vector("list", nrow(m))
  for (i in seq_len(nrow(m))) {
    a <- if (m[i, 1] < 0) c(leaf_y[-m[i, 1]], 0) else c(pos[m[i, 1]], tree$height[m[i, 1]])
    b <- if (m[i, 2] < 0) c(leaf_y[-m[i, 2]], 0) else c(pos[m[i, 2]], tree$height[m[i, 2]])
    pos[i] <- (a[1] + b[1]) / 2
    segs[[i]] <- data.frame(x = c(a[2], b[2], tree$height[i]), xend = tree$height[i],
                            y = c(a[1], b[1], a[1]), yend = c(a[1], b[1], b[1]))
  }
  do.call(rbind, segs)
}

pairHeatmap <- function(id) {
  out <- readRDS(file.path(align_dir, id, "combined.rds"))
  g <- as.data.frame(fread(file.path(align_dir, id, "gsea.csv")))
  prog <- out$programs[out$programs$stable %in% TRUE, , drop = FALSE]
  sid <- rownames(out$oof)
  role <- samples$role[match(sid, samples$profile_id)]

  ## Step 08's tree and clusters, checked against the clusters it wrote.
  transferring <- prog$program[prog$group %in% c("persists", "reduced_concordant")]
  if (length(transferring) < 2) transferring <- prog$program
  Zt <- standardizeActivity(out$oof[, transferring, drop = FALSE], role)
  tree <- stats::hclust(stats::dist(Zt), "ward.D2")
  cl <- stats::cutree(tree, k = k_clusters)
  saved <- as.data.frame(fread(file.path(align_dir, id, sprintf("sample_clusters_k%d.csv", k_clusters))))
  if (!identical(unname(cl[saved$sample]), as.integer(saved$cluster))) stop(id, ": clusters differ from step 08's")

  ## Paired samples only, each pair joined first within its cluster, with a gap where the cluster changes.
  pr <- pairs[pairs$organoid %in% sid & pairs$tumor %in% sid, ]
  keep <- unique(c(pr$tumor, pr$organoid))
  ptree <- pairedTree(Zt[keep, , drop = FALSE], cl[keep], pr)
  ord <- ptree$labels[ptree$order]
  block <- cumsum(c(1, diff(cl[ord]) != 0))
  pos <- seq_along(ord) + gap * (block - 1)
  y <- stats::setNames(max(pos) + 1 - pos, ord)
  r <- role[match(ord, sid)]
  side <- data.frame(y = y, role = factor(ifelse(r == "tumor", "Tumor", "Organoid"), levels = names(role_cols)))
  together <- table(factor(cl[pr$tumor][cl[pr$tumor] == cl[pr$organoid]], levels = seq_len(k_clusters)))
  lab <- data.frame(cluster = names(together), n = as.vector(together))
  lab <- lab[lab$cluster %in% cl[ord], ]   # a cluster of unpaired tumors only has no rows
  lab$y <- as.vector(tapply(y, cl[ord], function(v) mean(range(v)))[lab$cluster])
  lab$text <- sprintf("Cluster %s\n%d %s together", lab$cluster, lab$n, ifelse(lab$n == 1, "pair", "pairs"))

  ## Programs that persist, grouped by how they rank Reactome pathways, each group named by the pathway
  ## that sets it apart; groups run along the program tree.
  prog <- prog[prog$group == "persists", , drop = FALSE]
  pc <- pathwayClusters(g, prog$program, k = k_programs)
  first <- pc[!duplicated(pc$cluster), ]
  name_text <- ifelse(is.na(first$label), "No significant pathway", pathwayLabel(first$label, n = 80))
  head_text <- vapply(seq_along(name_text), function(i)   # about 7 characters per column
    paste(strwrap(name_text[i], max(9, 7 * first$n_programs[i])), collapse = "\n"), character(1))
  prog$pcluster <- factor(pc$cluster[match(prog$program, pc$program)], levels = first$cluster, labels = head_text)
  ## Activity is standardized within tumors and within organoids over all held-out samples, as for the tree.
  ## Within a group, programs follow the cluster (top to bottom) where they peak.
  Z <- standardizeActivity(out$oof[, prog$program, drop = FALSE], role)
  cmean <- rowsum(Z[ord, , drop = FALSE], block) / as.vector(table(block))
  peak <- apply(cmean[, prog$program, drop = FALSE], 2, which.max)
  prog <- prog[order(prog$pcluster, peak, -apply(cmean[, prog$program, drop = FALSE], 2, max)), ]
  hm <- data.frame(y = rep(y, times = ncol(Z)), program = rep(colnames(Z), each = length(ord)),
                   z = as.vector(Z[ord, ]))
  hm$label <- factor(hm$program, levels = prog$program)
  hm$pcluster <- prog$pcluster[match(hm$program, prog$program)]

  ys <- scale_y_continuous(limits = c(0.5, max(y) + 0.5), expand = c(0, 0))
  bare <- theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  p_lab <- ggplot(lab, aes(x = 1, y = y, label = text)) +
    geom_text(hjust = 1, size = 2.5, lineheight = 0.9, colour = ink[["primary"]]) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + ys + coord_cartesian(clip = "off") + labs(x = NULL) + bare +
    theme(axis.text.x = element_blank())
  p_tree <- ggplot(treeSegments(ptree, y[ptree$labels])) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.2, colour = ink[["secondary"]]) +
    scale_x_reverse(expand = c(0, 0)) + ys + labs(x = NULL) + bare + theme(axis.text.x = element_blank())
  p_role <- ggplot(side, aes(x = "Tumor or organoid", y = y)) + ys + labs(x = NULL) + bare +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = ink[["primary"]])) +
    geom_tile(aes(fill = role), height = 1) +
    scale_fill_manual(values = role_cols, name = NULL)
  p_hm <- ggplot(hm, aes(label, y, fill = pmax(pmin(z, 2.5), -2.5))) +
    geom_tile(height = 1) +
    scale_fill_gradient2(low = diverging[["low"]], mid = diverging[["mid"]], high = diverging[["high"]],
                         limits = c(-2.5, 2.5), breaks = c(-2.5, 0, 2.5), labels = c("-2.5 or less", "0", "2.5 or more"),
                         name = "Program activity\n(SD within tumors\nor within organoids)") +
    facet_grid(. ~ pcluster, scales = "free_x", space = "free_x") +
    ys + labs(x = NULL) + bare +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = ink[["primary"]], size = 6.5),
          strip.clip = "off", panel.spacing.x = grid::unit(6, "pt"))

  fig <- p_lab + p_tree + p_role + p_hm +
    plot_layout(widths = c(3, 2.5, 0.6, nrow(prog)), guides = "collect") &
    theme(legend.position = "right", plot.background = element_rect(fill = ink[["surface"]], colour = NA))
  name <- file.path(fig_dir, sprintf("pair_heatmap_%s", id))
  ggsave(paste0(name, ".pdf"), fig, width = 14, height = 12, device = grDevices::cairo_pdf)
  ggsave(paste0(name, ".png"), fig, width = 14, height = 12, dpi = 200)
  message("wrote ", name, ".{pdf,png}")
}
invisible(lapply(ids, pairHeatmap))
