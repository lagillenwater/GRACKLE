## Step 09: sample tree and program heatmap for the selected pair-alignment setting and the baseline.
##   Rscript analysis/hcmi_persistence/09_pair_heatmap.R
## Reads step 08's outputs under results/align (out-of-fold activity in combined.rds, the persistence map
## and the GSEA table); refits nothing. Rows are every held-out tumor and organoid, ordered by a Ward tree
## on all stable programs standardized within tumors and within organoids (standardizeActivity(); the
## clustering step 08 compares across settings). Beside the tree: each sample's role and one color per
## tumor-organoid pair (samples linked through a shared tumor or organoid share it; unpaired tumors are
## blank), so a pair the tree places together shows as adjacent rows of one color. Columns are the stable
## programs, named by their leading Reactome pathway (preranked GSEA) and grouped by persistence call.
source("analysis/hcmi_persistence/config.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
invisible(parseArgs(list()))
align_dir <- file.path(cfg$run, "results", "align")
fig_dir <- file.path(align_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
samples <- readRDS(file.path(cfg$processed, "samples_included.rds"))
pairs <- pairTable(samples)
tab <- as.data.frame(fread(file.path(align_dir, "settings.csv")))
ids <- unique(c(tab$setting_id[tab$selected], tab$setting_id[tab$baseline]))
t1 <- tab[match(ids, tab$setting_id), ]
title <- stats::setNames(ifelse(t1$baseline, "plain NMF (baseline)",
                                sprintf("GRACKLE, %s sample graph alpha_1 = %s, %s gene graph alpha_2 = %s",
                                        t1$sample_graph, t1$alpha_1, t1$arm, t1$alpha_2)), ids)

## Colors from the dataviz reference palette: diverging blue-red with its neutral midpoint for activity,
## fixed limits so the two figures compare directly; neutral inks for role so the pair colors stand out.
ink <- c(primary = "#0b0b0b", secondary = "#52514e", baseline = "#c3c2b7", surface = "#fcfcfb")
diverging <- c(low = "#2a78d6", mid = "#f0efec", high = "#e34948")
role_cols <- c(Tumor = ink[["secondary"]], Organoid = ink[["baseline"]])
group_levels <- c(persists = "Persists", reduced_concordant = "Reduced,\nconcordant",
                  present_not_patient_specific = "Present, not\npatient-specific", drops_out = "Drops out")
theme_set(theme_minimal(base_size = 8) +
            theme(panel.grid = element_blank(), text = element_text(colour = ink[["primary"]]),
                  axis.text = element_text(colour = ink[["secondary"]]),
                  strip.text = element_text(colour = ink[["primary"]], size = 7),
                  plot.background = element_rect(fill = ink[["surface"]], colour = NA)))

## Reactome names for display: prefix dropped, sentence case with acronyms kept, long names cut.
acronyms <- c("ADME", "ECM", "HCMV", "RNA", "DNA", "SLC", "GPI", "PCP", "PD", "II", "III", "SPMS", "HFTC", "TP53", "NOTCH1",
              "AR", "KLK2", "KLK3", "PKN1", "MLL3", "MLL4", "GALNT3", "CD274", "L1", "NMD", "TCA", "ER", "MET", "O")
pathwayLabel <- function(x, n = 40) {
  w <- strsplit(sub("^REACTOME_", "", x), "_", fixed = TRUE)
  s <- vapply(w, function(t) {
    t <- ifelse(t %in% acronyms, t, tolower(t))
    t[1] <- paste0(toupper(substr(t[1], 1, 1)), substring(t[1], 2))
    paste(t, collapse = " ")
  }, character(1))
  ifelse(nchar(s) > n, paste0(substr(s, 1, n - 1), "…"), s)
}

## Tree segments with leaves at x = 0 and leaf i of tree$order at y = n + 1 - i (first leaf on top).
## Merges are processed in order, so each child's position is known before its parent's.
treeSegments <- function(tree) {
  n <- length(tree$order)
  leaf_y <- integer(n)
  leaf_y[tree$order] <- n + 1 - seq_len(n)
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

## One color per pair, assigned down the tree with hues a golden angle apart and three lightness steps,
## so neighboring pairs never share a look.
pairColors <- function(k) {
  i <- seq_len(k)
  grDevices::hcl(h = (i * 137.508) %% 360, c = 65, l = c(52, 72, 62)[(i %% 3) + 1])
}

pairHeatmap <- function(id) {
  out <- readRDS(file.path(align_dir, id, "combined.rds"))
  g <- as.data.frame(fread(file.path(align_dir, id, "gsea.csv")))
  prog <- out$programs[out$programs$stable %in% TRUE, , drop = FALSE]
  sid <- rownames(out$oof)
  role <- samples$role[match(sid, samples$profile_id)]
  Z <- standardizeActivity(out$oof[, prog$program, drop = FALSE], role)
  tree <- stats::hclust(stats::dist(Z), method = "ward.D2")
  n <- length(sid)
  y <- integer(n)
  y[tree$order] <- n + 1 - seq_len(n)

  ## Pairs: connected tumors and organoids share a color, numbered in tree order.
  pr <- pairs[pairs$organoid %in% sid & pairs$tumor %in% sid, ]
  comp <- igraph::components(igraph::graph_from_data_frame(pr[c("tumor", "organoid")], directed = FALSE))$membership
  first <- unique(comp[intersect(sid[tree$order], names(comp))])
  pair_rank <- match(comp[sid], first)
  side <- data.frame(y = y, role = factor(ifelse(role == "tumor", "Tumor", "Organoid"), levels = names(role_cols)),
                     pair = ifelse(is.na(pair_rank), NA, pairColors(length(first))[pair_rank]))

  ## Programs: leading pathway by GSEA, grouped by persistence call, ordered by a tree within each group.
  g <- g[g$program %in% prog$program & g$fdr < 0.05, ]
  g <- g[order(g$program, g$fdr, -g$nes), ]
  lead <- g[!duplicated(g$program), ]
  prog$label <- paste(prog$program, ifelse(prog$program %in% lead$program,
                                           pathwayLabel(lead$gene_set[match(prog$program, lead$program)]), "(no pathway)"))
  col_tree <- stats::hclust(stats::dist(t(Z)), method = "ward.D2")
  prog$group <- factor(prog$group, levels = names(group_levels), labels = group_levels)
  prog <- prog[order(prog$group, match(prog$program, colnames(Z)[col_tree$order])), ]
  hm <- data.frame(y = rep(y, times = ncol(Z)), program = rep(colnames(Z), each = n), z = as.vector(Z))
  hm$label <- factor(prog$label[match(hm$program, prog$program)], levels = prog$label)
  hm$group <- prog$group[match(hm$program, prog$program)]

  ys <- scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0))
  bare <- theme(axis.text.y = element_blank(), axis.title.y = element_blank())
  p_tree <- ggplot(treeSegments(tree)) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.2, colour = ink[["secondary"]]) +
    scale_x_reverse(expand = c(0, 0)) + ys + labs(x = NULL) + bare + theme(axis.text.x = element_blank())
  strip <- function(fill_col, label) {
    ggplot(side, aes(x = label, y = y)) + ys + labs(x = NULL) + bare +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = ink[["primary"]]))
  }
  p_role <- strip("role", "Tumor or organoid") +
    geom_tile(aes(fill = role), height = 1) +
    scale_fill_manual(values = role_cols, name = NULL)
  p_pair <- strip("pair", "Pair") +
    geom_tile(data = side[!is.na(side$pair), ], aes(fill = I(pair)), height = 1)
  p_hm <- ggplot(hm, aes(label, y, fill = pmax(pmin(z, 2.5), -2.5))) +
    geom_tile(height = 1) +
    scale_fill_gradient2(low = diverging[["low"]], mid = diverging[["mid"]], high = diverging[["high"]],
                         limits = c(-2.5, 2.5), breaks = c(-2.5, 0, 2.5), labels = c("-2.5 or less", "0", "2.5 or more"),
                         name = "Program activity\n(SD within tumors\nor within organoids)") +
    facet_grid(. ~ group, scales = "free_x", space = "free_x") +
    ys + labs(x = NULL) + bare +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = ink[["primary"]], size = 6.5),
          strip.clip = "off", panel.spacing.x = grid::unit(3, "pt"))

  fig <- (p_tree + p_role + p_pair + p_hm + plot_layout(widths = c(4, 0.4, 0.4, nrow(prog)), guides = "collect")) +
    plot_annotation(title = sprintf("Tumors and organoids by program activity: %s", title[[id]])) &
    theme(legend.position = "right", plot.background = element_rect(fill = ink[["surface"]], colour = NA))
  name <- file.path(fig_dir, sprintf("pair_heatmap_%s", id))
  ggsave(paste0(name, ".pdf"), fig, width = 12, height = 12, device = grDevices::cairo_pdf)
  ggsave(paste0(name, ".png"), fig, width = 12, height = 12, dpi = 200)
  message("wrote ", name, ".{pdf,png}")
}
invisible(lapply(ids, pairHeatmap))
