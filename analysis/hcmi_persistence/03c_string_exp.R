## Step 03c: a modular STRING gene graph for the pair-alignment run (string_exp).
## The experiments channel alone (the curated-database channel imports pathways such as Reactome), with
## edges at score >= --string_min (default 0.4, STRING's "medium confidence"). Scored by GSEA, this graph
## raised pathway enrichment in a fold-0 test, where coess_2d lowered it.
source("analysis/hcmi_persistence/config.R")
args <- parseArgs(list(string_min = "0.4"))
genes <- readLines(file.path(cfg$processed, "genes.txt"))
hgnc_file <- file.path(cfg$raw, "annotation", "hgnc_complete_set_2024-10-01.txt")
string_files <- file.path(cfg$raw, "string", c("9606.protein.links.full.v12.0.txt.gz", "9606.protein.aliases.v12.0.txt.gz"))
gg_dir <- file.path(cfg$processed, "gene_graphs")
links <- readStringLinks(string_files[1])
pmap <- stringProteinToSymbol(readStringAliases(string_files[2]), as.data.frame(fread(hgnc_file, quote = "")))
G <- stringGraph(links, pmap, genes, "exp", min_score = as.numeric(args$string_min))
saveRDS(G, file.path(gg_dir, "string_exp.rds"))
writeProvenance(file.path(gg_dir, "string_exp.rds"), inputs = c(string_files, hgnc_file),
                extra = list(min_score = args$string_min))
print(cbind(arm = "string_exp", graphSummary(G)))
