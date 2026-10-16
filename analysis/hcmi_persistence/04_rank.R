## Step 04: rank k by Horn's parallel analysis (spec Section 7.1).
source("analysis/hcmi_persistence/config.R")
Y <- readRDS(file.path(cfg$processed, "Y.rds"))
pa <- parallelAnalysisRank(Y, n_perm = 100, quantile = 0.95, nv = 200, seed = cfg$seed)
saveRDS(pa$k, file.path(cfg$processed, "k.rds"))
fwrite(data.frame(component = seq_along(pa$observed), observed = pa$observed, threshold = pa$threshold),
       file.path(cfg$processed, "k_parallel_analysis.csv"))
writeProvenance(file.path(cfg$processed, "k.rds"), inputs = file.path(cfg$processed, "Y.rds"), extra = list(k = pa$k))
message("k = ", pa$k)
