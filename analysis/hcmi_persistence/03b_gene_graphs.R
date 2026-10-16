## Step 03b: the four gene-graph arms (spec Section 4.5). coess_3d is rebuilt per fold.
source("analysis/hcmi_persistence/config.R")
genes <- readLines(file.path(cfg$processed, "genes.txt"))
samples <- readRDS(file.path(cfg$processed, "samples_included.rds"))
hgnc <- as.data.frame(fread(file.path(cfg$raw, "annotation", "hgnc_complete_set_2024-10-01.txt"), quote = ""))
gg_dir <- file.path(cfg$processed, "gene_graphs"); dir.create(gg_dir, recursive = TRUE, showWarnings = FALSE)
summ <- list()
record <- function(arm, fold, G, n_models = NA, n_informative = NA, decision = NA) {
  summ[[length(summ) + 1L]] <<- cbind(data.frame(arm = arm, fold = fold, n_models = n_models,
                                                 n_informative = n_informative, decision = decision), graphSummary(G))
}

## STRING arms
links <- readStringLinks(file.path(cfg$raw, "string", "9606.protein.links.full.v12.0.txt.gz"))
pmap <- stringProteinToSymbol(readStringAliases(file.path(cfg$raw, "string", "9606.protein.aliases.v12.0.txt.gz")), hgnc)
invisible(lapply(c("mech", "all"), function(mode) {
  G <- stringGraph(links, pmap, genes, mode)
  saveRDS(G, file.path(gg_dir, sprintf("string_%s.rds", mode)))
  record(paste0("string_", mode), NA, G)
}))
rm(links); invisible(gc())

## Co-essentiality arms. Arm membership is the screen's culture format (spec Section 4.5).
dm_dir <- file.path(cfg$raw, "depmap")
models24 <- as.data.frame(fread(file.path(dm_dir, "Model.csv")))
ng_models <- as.data.frame(fread(file.path(dm_dir, "nextgen_model_metadata.csv")))
ng_screens <- as.data.frame(fread(file.path(dm_dir, "nextgen_screen_metadata.csv")))
ng <- as.data.frame(fread(file.path(cfg$analysis, "depmap_nextgen_files.csv")))
ng_path <- function(role) file.path(cfg$manual, "depmap_nextgen", ng$filename[ng$role == role])
model_tables <- list(models24)
if ("models" %in% ng$role) model_tables <- c(model_tables, list(as.data.frame(fread(ng_path("models")))))
tokens <- rbindlist(lapply(model_tables, function(m) data.table(ModelID = m$ModelID, hcmi = hcmiTokenFromRow(m))))
tokens <- tokens[!is.na(tokens$hcmi)]
hcmi_of <- function(ids) tokens$hcmi[match(ids, tokens$ModelID)]
nonessential <- readNonessential(file.path(dm_dir, "AchillesNonessentialControls.csv"), hgnc)
arm <- stats::setNames(geneArm(hgnc$location), hgnc$symbol)

E_ng <- readDepmapEffect(ng_path("gene_effect"), hgnc)
cls <- screenClass(colnames(E_ng), ng_screens, ng_models)
fwrite(cls, file.path(cfg$processed, "depmap_screen_classes.csv"))
E_ng2 <- collapseScreens(E_ng[, cls$class == "2D", drop = FALSE], cls$model_id[cls$class == "2D"])
E3_dm <- collapseScreens(E_ng[, cls$class == "3D", drop = FALSE], cls$model_id[cls$class == "3D"])
## A model in no model table cannot be checked for HCMI origin, so it is dropped (ledger, Task 22).
known_models <- unique(unlist(lapply(model_tables, `[[`, "ModelID")))
id2 <- dropUnidentifiedModels(E_ng2, known_models); E_ng2 <- id2$effect
id3 <- dropUnidentifiedModels(E3_dm, known_models); E3_dm <- id3$effect
message("models in no model table, dropped: 2D ", length(id2$dropped), ", 3D ", length(id3$dropped),
        if (length(id3$dropped)) paste0(" (", paste(id3$dropped, collapse = ", "), ")"))
if (ncol(E_ng2) >= 1000) {
  E2 <- E_ng2
  source_2d <- sprintf("NextGen manuscript run, shared with 3D (%d traditional models)", ncol(E_ng2))
} else {
  E24 <- readDepmapEffect(file.path(dm_dir, "CRISPRGeneEffect.csv"), hgnc)
  E2 <- E24[, intersect(colnames(E24), models24$ModelID[models24$GrowthPattern %in% c("Adherent", "Suspension", "Mixed")]), drop = FALSE]
  source_2d <- sprintf("DepMap 24Q4, a separate run (the NextGen file has %d traditional models)", ncol(E_ng2))
}
E2 <- E2[, is.na(hcmi_of(colnames(E2))), drop = FALSE]
G2 <- coessGraph(E2, genes, nonessential, arm)
saveRDS(G2, file.path(gg_dir, "coess_2d.rds"))
record("coess_2d", NA, G2, ncol(E2), length(informativeGenes(E2, nonessential)), source_2d)

E3_sg <- NULL
sg_tab <- file.path(cfg$analysis, "sanger_files.csv")
if (file.exists(sg_tab)) {
  sg <- as.data.frame(fread(sg_tab))
  row <- sg$role == "gene_effect_3d"
  E3_sg <- readDepmapEffect(file.path(cfg$manual, "sanger", sg$filename[row]), hgnc,
                            genes_in_rows = "genes_in_rows" %in% names(sg) && isTRUE(as.logical(sg$genes_in_rows[row])))
  colnames(E3_sg) <- paste0("SANGER:", colnames(E3_sg))
}
combined_3d <- combine3dSources(E3_dm, E3_sg)
E3 <- combined_3d$effect
model_token <- c(stats::setNames(hcmi_of(colnames(E3_dm)), colnames(E3_dm)),
                 if (!is.null(E3_sg)) stats::setNames(hcmiToken(sub("^SANGER:", "", colnames(E3_sg))), colnames(E3_sg)))
model_case <- stats::setNames(samples$case_id[match(model_token, samples$case_submitter_id, incomparables = NA)], names(model_token))
model_case <- model_case[colnames(E3)]
message("3D models: ", ncol(E3), "; mapped to an included HCMI case: ", sum(!is.na(model_case)))
invisible(lapply(0:5, function(f) {
  held <- if (f == 0) character(0) else unique(samples$case_id[samples$fold == f])
  Ef <- dropHeldOutModels(E3, model_case, held)
  G3 <- coessGraph(Ef, genes, nonessential, arm)
  saveRDS(G3, file.path(gg_dir, sprintf("coess_3d_fold%d.rds", f)))
  record("coess_3d", f, G3, ncol(Ef), length(informativeGenes(Ef, nonessential)),
         sprintf("%s; %d unidentified DepMap 3D models dropped", combined_3d$decision, length(id3$dropped)))
}))
fwrite(rbindlist(summ, fill = TRUE), file.path(cfg$processed, "gene_graph_summary.csv"))
writeProvenance(file.path(cfg$processed, "gene_graph_summary.csv"),
                inputs = c(file.path(cfg$processed, "genes.txt"), ng_path("gene_effect"), file.path(dm_dir, "Model.csv")))
print(rbindlist(summ, fill = TRUE))
