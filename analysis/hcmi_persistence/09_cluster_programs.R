## Step 09: heatmap of the sample clusters and the enriched programs they carry, for the selected
## pair-alignment setting and the baseline.
##   Rscript analysis/hcmi_persistence/09_cluster_programs.R [--k_clusters 8]
## Reads step 08's outputs under results/align (out-of-fold activity in combined.rds, the sample clusters,
## the persistence map and the GSEA table); refits nothing. Left: the mean activity of each cluster's
## tumors and organoids for every stable program (clusterActivity(): log1p, standardized over all
## samples), boxed where the program transfers in that cluster. Middle: microenvironment programs.
## Right: each program's significant Reactome pathways (preranked GSEA, FDR < 0.05) among the leading
## pathways of all programs, shaded by NES. Rows are grouped by persistence call.
source("analysis/hcmi_persistence/config.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
args <- parseArgs(list(k_clusters = "8"))
k_clusters <- as.integer(args$k_clusters)
align_dir <- file.path(cfg$run, "results", "align")
fig_dir <- file.path(align_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
tab <- as.data.frame(fread(file.path(align_dir, "settings.csv")))
ids <- unique(c(tab$setting_id[tab$selected], tab$setting_id[tab$baseline]))
t1 <- tab[match(ids, tab$setting_id), ]
title <- stats::setNames(ifelse(t1$baseline, "plain NMF (baseline)",
                                sprintf("GRACKLE, %s sample graph alpha_1 = %s, %s gene graph alpha_2 = %s",
                                        t1$sample_graph, t1$alpha_1, t1$arm, t1$alpha_2)), ids)

## Colors from the dataviz reference palette: diverging blue-red with its neutral midpoint for activity,
## a one-hue violet ramp for NES (kept apart from both diverging poles); ink and hairlines from its chrome.
## Both scales are fixed so the selected and baseline figures compare directly.
ink <- c(primary = "#0b0b0b", secondary = "#52514e", muted = "#898781", hairline = "#e1e0d9", surface = "#fcfcfb")
diverging <- c(low = "#2a78d6", mid = "#f0efec", high = "#e34948")
violet <- c(low = "#e4e0fa", high = "#3a2d9c")
group_levels <- c(persists = "Persists", reduced_concordant = "Reduced,\nconcordant",
                  present_not_patient_specific = "Present, not\npatient-specific", drops_out = "Drops out")
theme_set(theme_minimal(base_size = 8) +
            theme(panel.grid = element_blank(), text = element_text(colour = ink[["primary"]]),
                  axis.text = element_text(colour = ink[["secondary"]]),
                  strip.text = element_text(colour = ink[["primary"]], size = 7),
                  panel.spacing = grid::unit(2, "pt"),
                  plot.background = element_rect(fill = ink[["surface"]], colour = NA)))

## Reactome names for display: prefix dropped, sentence case with acronyms kept, long names cut.
acronyms <- c("ADME", "ECM", "HCMV", "RNA", "DNA", "SLC", "GPI", "PCP", "PD", "II", "III", "SPMS", "HFTC", "TP53", "NOTCH1",
              "AR", "KLK2", "KLK3", "PKN1", "MLL3", "MLL4", "GALNT3", "CD274", "L1", "NMD", "TCA", "ER", "MET", "O")
pathwayLabel <- function(x, n = 44) {
  w <- strsplit(sub("^REACTOME_", "", x), "_", fixed = TRUE)
  s <- vapply(w, function(t) {
    t <- ifelse(t %in% acronyms, t, tolower(t))
    t[1] <- paste0(toupper(substr(t[1], 1, 1)), substring(t[1], 2))
    paste(t, collapse = " ")
  }, character(1))
  ifelse(nchar(s) > n, paste0(substr(s, 1, n - 1), "…"), s)
}

clusterFigure <- function(id) {
  d <- file.path(align_dir, id)
  out <- readRDS(file.path(d, "combined.rds"))
  sc <- as.data.frame(fread(file.path(d, sprintf("sample_clusters_k%d.csv", k_clusters))))
  pc <- as.data.frame(fread(file.path(d, sprintf("clusters_k%d.csv", k_clusters))))
  g <- as.data.frame(fread(file.path(d, "gsea.csv")))
  prog <- out$programs[out$programs$stable %in% TRUE, , drop = FALSE]
  sc <- sc[match(rownames(out$oof), sc$sample), ]
  long <- clusterActivity(out$oof[, prog$program, drop = FALSE], sc$role, sc$cluster)
  fwrite(long, file.path(d, sprintf("cluster_activity_k%d.csv", k_clusters)))

  ## Clusters in the order of a Ward tree on their tumor and organoid profiles; programs grouped by
  ## persistence call, then by the cluster where they peak, so each cluster's programs sit together.
  wide <- cbind(tapply(long$tumor_z, long[c("cluster", "program")], identity),
                tapply(long$organoid_z, long[c("cluster", "program")], identity))
  corder <- rownames(wide)[stats::hclust(stats::dist(wide), "ward.D2")$order]
  long$cluster <- factor(long$cluster, levels = corder)
  peak <- long[order(long$program, -(long$tumor_z + long$organoid_z)), ]
  peak <- peak[!duplicated(peak$program), ]
  prog$group <- factor(prog$group, levels = names(group_levels), labels = group_levels)
  prog$peak <- match(as.character(peak$cluster[match(prog$program, peak$program)]), corder)
  prog$peak_z <- (peak$tumor_z + peak$organoid_z)[match(prog$program, peak$program)]
  prog <- prog[order(prog$group, prog$peak, -prog$peak_z), ]

  ## Each program's leading pathway (most significant by GSEA) names its row.
  g <- g[g$program %in% prog$program & g$fdr < 0.05, ]
  g <- g[order(g$program, g$fdr, -g$nes), ]
  lead <- g[!duplicated(g$program), ]
  prog$label <- paste(prog$program, ifelse(prog$program %in% lead$program,
                                           pathwayLabel(lead$gene_set[match(prog$program, lead$program)]), "(no pathway)"))
  row_levels <- rev(prog$label)
  rows <- function(x) {
    x$label <- factor(prog$label[match(x$program, prog$program)], levels = row_levels)
    x$group <- prog$group[match(x$program, prog$program)]
    x
  }

  ## Activity: one column each for a cluster's tumors and organoids.
  act <- rows(rbind(data.frame(long[c("cluster", "program", "status")], role = "T", z = long$tumor_z),
                    data.frame(long[c("cluster", "program", "status")], role = "O", z = long$organoid_z)))
  act$role <- factor(act$role, levels = c("T", "O"))
  box <- act[act$role == "T" & act$status != "", ]
  box$status <- factor(box$status, levels = c("transfers", "lost"),
                       labels = c("Transfers (tumors and organoids > 0.5)", "Lost (tumors > 0.5, organoids < 0)"))
  pc <- pc[match(corder, pc$cluster), ]
  strip <- stats::setNames(sprintf("C%s\n%d T\n%d O\n%d pairs", pc$cluster, pc$n_tumors, pc$n_organoids, pc$pairs_together),
                           pc$cluster)
  lim <- 1.5
  p_act <- ggplot(act, aes(role, label)) +
    geom_tile(aes(fill = pmax(pmin(z, lim), -lim)), colour = ink[["surface"]], linewidth = 0.25) +
    geom_tile(data = box, aes(x = 1.5, linetype = status), width = 1.9, height = 0.86, fill = NA,
              colour = ink[["primary"]], linewidth = 0.45) +
    scale_fill_gradient2(low = diverging[["low"]], mid = diverging[["mid"]], high = diverging[["high"]],
                         limits = c(-lim, lim), breaks = c(-lim, 0, lim), labels = c("-1.5 or less", "0", "1.5 or more"),
                         name = "Mean activity\n(SD, all samples)") +
    scale_linetype_manual(values = c("solid", "22"), name = NULL, drop = FALSE) +
    facet_grid(group ~ cluster, scales = "free_y", space = "free_y", labeller = labeller(cluster = strip)) +
    labs(x = "Sample cluster: tumors (T) and organoids (O)", y = NULL) +
    theme(strip.text.y = element_blank(), axis.text.y = element_text(colour = ink[["primary"]], size = 6.5),
          legend.key = element_rect(fill = ink[["surface"]], colour = NA))

  me <- rows(data.frame(program = prog$program, microenv = prog$microenv %in% TRUE, x = "Micro-\nenv."))
  p_me <- ggplot(me, aes(x, label)) +
    geom_tile(fill = NA) +
    geom_point(data = me[me$microenv, ], size = 1.6, colour = ink[["secondary"]]) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme(axis.text.y = element_blank(), strip.text.y = element_blank(),
          panel.background = element_rect(fill = NA, colour = ink[["hairline"]], linewidth = 0.3))

  ## Pathways: the leading pathway of every program, ordered down the rows.
  sets <- unique(lead$gene_set[order(match(lead$program, prog$program))])
  nes <- rows(g[g$gene_set %in% sets, c("program", "gene_set", "nes")])
  nes$gene_set <- factor(nes$gene_set, levels = sets, labels = pathwayLabel(sets, 40))
  frame <- rows(expand.grid(program = prog$program, gene_set = levels(nes$gene_set), stringsAsFactors = FALSE))
  frame$gene_set <- factor(frame$gene_set, levels = levels(nes$gene_set))
  p_nes <- ggplot(nes, aes(gene_set, label)) +
    geom_tile(data = frame, fill = ink[["surface"]], colour = ink[["hairline"]], linewidth = 0.2) +
    geom_tile(aes(fill = nes), colour = ink[["surface"]], linewidth = 0.25) +
    scale_fill_gradient(low = violet[["low"]], high = violet[["high"]], limits = c(1, 2.2), oob = scales::squish,
                        name = "GSEA NES\n(FDR < 0.05)") +
    scale_x_discrete(position = "top", drop = FALSE) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Leading Reactome pathways", y = NULL) +
    theme(axis.text.y = element_blank(), axis.text.x.top = element_text(angle = 50, hjust = 0, vjust = 0, size = 6),
          strip.text.y = element_text(angle = 0, hjust = 0))

  n_c <- nlevels(long$cluster)
  fig <- (p_act + p_me + p_nes + plot_layout(widths = c(2 * n_c, 1.2, length(sets)), guides = "collect")) +
    plot_annotation(title = sprintf("Sample clusters and the programs they carry: %s", title[[id]]),
                    subtitle = sprintf("%s. %d stable programs; %d clusters of held-out samples (T = tumors, O = organoids); %d of %d pairs share a cluster",
                                       id, nrow(prog), n_c, sum(pc$pairs_together), nrow(out$matches))) &
    theme(legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(size = 7),
          plot.background = element_rect(fill = ink[["surface"]], colour = NA))
  name <- file.path(fig_dir, sprintf("cluster_programs_%s", id))
  ggsave(paste0(name, ".pdf"), fig, width = 13, height = 9.5, device = grDevices::cairo_pdf)
  ggsave(paste0(name, ".png"), fig, width = 13, height = 9.5, dpi = 200)
  message("wrote ", name, ".{pdf,png}")
}
invisible(lapply(ids, clusterFigure))
