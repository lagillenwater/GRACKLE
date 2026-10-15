# HCMI Tumor-to-Organoid Program Persistence Map: Design

- **Date:** 2026-10-15
- **Branch:** `organoid-transfer` (cut from `alpine` at `1395a5d`)
- **Status:** design approved in conversation, spec awaiting review
- **Context:** K99 Specific Aim 1, preliminary result 1 ("GRACKLE runs on matched pairs and evaluates for gene set enrichment of latent variables")

---

## 1. Purpose and scope

**Question.** Which GRACKLE gene programs in parent tumors persist into the patient's own 3D organoid? Which drop out? And which are present in organoids without tracking the patient?

**Product.** A per-program persistence map for HCMI matched tumor/3D-organoid pairs, built by a leakage-free cross-fitting pipeline that runs on CU Alpine, plus the GRACKLE fixes that pipeline needs.

**In scope**
- HCMI data download, sample table, preprocessing.
- Sample-graph construction from priors.
- Rank selection and lambda tuning.
- 5-fold cross-fitting with cross-fold program matching.
- Persistence scoring and post-hoc annotation.
- Figures and tables.
- Alpine orchestration.
- Unit tests and a synthetic end-to-end fixture.

**Out of scope** (later pieces, which reuse this piece's programs and folds)
- Drug response and the program-response invariance tests.
- Derivation-propensity modeling of the 168 failed derivations.
- The cell-line and proliferation arm (GDSC2/PRISM).
- The shared/source-specific factorization.
- CI, cleanup of the existing `tests/*.R` analysis scripts, and reconciling `alpine` with `main`.

**Success criteria**
1. Every output traces to a commit, a Slurm job and input checksums. The full chain, from download to figure, reruns from the branch alone.
2. The synthetic fixture places each planted program in its correct group (Section 8.4) before any Alpine submission.
3. The main figure, the per-program table and the supplementary sensitivity panels exist for the selected setting.
4. A filter-count table reconciles the HCMI pair counts used in the grant text (Section 2).

The analysis has no target outcome. Whatever fraction of programs persists is the result.

---

## 2. Background facts (verified 2026-10-15)

From ElHarouni et al. 2026, *Nature* (doi:10.1038/s41586-026-10806-y):
- **Models:** 665 models from 637 patients. 519 are 3D organoids, 37 are 3D spheroids and 109 are 2D adherent lines.
- **Matched tumors:** 421 models have matched parental tumors. This is the DNA-level set, with 97.8% genetic concordance.
- **Transcriptional concordance:** assessed on 297 tumor-model pairs (Celligner). The "92%" is 263 of the 286 models assessed by at least 4 of 6 methods.
- **Failed derivations:** tumors from 168 failed derivation attempts were profiled. In COAD, ESCA and PAAD, failed tumors showed more DNA-repair and cell-cycle states, and successful ones more differentiation.
- **Culture medium:** affects transcriptional fidelity. For example, in GBM, NS-A medium beats conditioned medium.
- **Data locations:** the GDC project `HCMI-CMDC` (open STAR counts) and the GDC publication page `HCMI-CMDC-2026` (Supplementary Table 1: model type, medium, derivation center, treatment, clinical). Mutations and clinical data are also in cBioPortal `pancan_hcmi_2025`.

From the GDC API (2026-10-15):
- 1,153 open-access STAR-count RNA-seq files.
- 248 cases have both a tumor-type and a model-type RNA profile.
- GDC sample fields lack model type, medium and center. Passage count is present for about 370 model samples.

**Grant text correction.** Research Plan A.a.i.1 says "92% transcriptional concordance across 421 matched tumor-model pairs." The 92% applies to 286 models with RNA. The 421 is the DNA-level set.

---

## 3. Decisions and rejected alternatives

| Decision | Choice | Rejected alternative and why |
|---|---|---|
| First piece of Aim 1 | HCMI paired GRACKLE with persistence scoring | Others depend on its programs and folds |
| Headline claim | Program persistence map | Method benchmark (pairing vs. no pairing): secondary to the hypothesis |
| Models | 3D organoids only | All model types: the user chose to match the proposal's organoid framing |
| Purity handling | Purity as a sample-graph prior, tumor-tumor edges only; post-hoc labeling of microenvironment programs | Purity-aware mixture model (fixed W columns, purity-scaled WH): **rejected by the user for loss of GRACKLE's explainability**. Purity inside a combined distance with organoids: pushes each tumor away from its own organoid |
| Sample-graph priors | Pair edges, cancer type, clinical context, purity (tumor-tumor) | Genotype: kept out of the fit so later genotype-program-response tests stay clean |
| Culture variables (medium, passage, center) | Excluded from the graph | They are the environments of the later invariance test; building them into the fit would bias it |
| RNA-derived labels (TMP/MOMA subtype, ESTIMATE) | Annotation only | Circular as priors |
| Training samples | 3D-organoid pairs plus all other HCMI tumors with RNA, including failed derivations | Pairs only: programs would describe only tumors that grew into organoids |
| FFPE tumor profiles | Used only when no frozen profile exists for that tumor sample (added during spec writing) | Always include: the paper reports FFPE artifacts behind some low-concordance pairs |
| Evaluation design | 5-fold cross-fitting with Hungarian program matching | Single split: too few pairs. Full-data fit plus reprojection: each scored pair's own edge shapes H |
| Rank k | Horn's parallel analysis, fixed in advance | A k sweep: any chosen list of values is arbitrary |
| Lambda tuning | Dimensionless alpha grid; one-SE rule on held-out reconstruction, then maximum reproducibility | Held-out reconstruction alone: always picks lambda = 0. Nested CV on persistence: would reward settings that make programs look shared |
| GRACKLE objective | Unchanged | (See purity row) |

---

## 4. Data

All steps run on Alpine (Section 11). Raw data goes in `/projects/$USER/grackle_data/hcmi/raw/` and intermediates in `.../processed/`.

### 4.1 Download (`01_download.R`)
- **GDC API:** every `HCMI-CMDC` file with `data_type = "Gene Expression Quantification"` and workflow `STAR - Counts`, plus case, sample and aliquot metadata. Record the GDC data release, the download date and each file's md5.
- **Supplementary Table 1** from the GDC publication page `HCMI-CMDC-2026`. If it is not hosted there, the user places the Nature supplementary file in `raw/manual/`; the script verifies it is present and records its sha256.
- **cBioPortal `pancan_hcmi_2025`:** mutation calls and sample clinical data, used for genotype annotation only.
- **STRING v12** (`9606.protein.links.v12.0.txt`, `9606.protein.aliases.v12.0.txt`), the same files as the HTP analyses.
- **MSigDB Hallmark** via `msigdbr`, with the package version pinned and the gene sets cached as an RDS file.

### 4.2 Sample table (`02_sample_table.R`)
- **One row per RNA profile.** Columns:
  - identity: `profile_id`, `aliquot_id`, `sample_id`, `case_id`, `role` (tumor / organoid)
  - tumor context: `cancer_type`, `tumor_descriptor`
  - model: `model_id`, `model_type`, `medium`, `passage`, `derivation_center`
  - pairing: `parent_sample_id`, `pair_id`
  - covariates: `purity_dna`, `prior_treatment`, `age`, `sex`, `ancestry`
  - bookkeeping: `fold`, `included`, `exclusion_reason`
- **Included organoids:** profiles whose model type is "3D organoid" in Supplementary Table 1.
- **Included tumors:** all tumor profiles with RNA. That covers primary, metastatic and recurrent tumors, human tumor original cells, and failed-derivation tumors.
- **FFPE:** a profile from an FFPE-preserved sample (`preservation_method == "FFPE"`, including slides and FFPE scrolls) is included only when no frozen profile exists for the same tumor sample. The paper attributes some low concordance to FFPE artifacts. The preservation method is kept as a column.
- **Pairing:**
  - An organoid pairs with the tumor sample it was derived from (`parent_sample_id` from Supplementary Table 1).
  - If the parent sample is not recorded, it pairs with every tumor sample of the case, and the case is flagged.
  - ATCC re-expanded models ("Expanded Next Generation Cancer Model") are excluded as duplicates.
- **Filter-count table:** the script writes one row per filter step with profiles, cases and pairs remaining (`filter_counts.csv`).
- **Committed:** `sample_table.csv` and `filter_counts.csv`, which are small and open-access.

### 4.3 Folds (`02_sample_table.R`)
- 5 folds with seed `20261015`, stratified by cancer type and grouped by case. Every case gets a fold, including unpaired tumors.
- Committed as `folds.csv`.

### 4.4 Expression (`03_expression.R`)
1. Take `tpm_unstranded` from each STAR-count file and compute log2(TPM + 1). Keep `gene_type == "protein_coding"`.
2. Select the top 5,000 genes by variance across all included samples.
3. Keep only genes that appear in the STRING-derived adjacency, following the HTP convention. The final gene count is at most 5,000 and is recorded.
4. Min-max scale each gene across all included samples (`min_max_scale`).

Gene selection and scaling use all samples but no pairing or annotation information. The methods section states this.

### 4.5 Gene graph
- STRING v12 `combined_score / 1000`, with proteins mapped to HGNC symbols through `Ensembl_HGNC_symbol` aliases and restricted to the final genes.
- The diagonal is set to 1. This is unchanged from `tests/HTP_v2.R`.

### 4.6 Purity requirement
- Graph purity must be **DNA-based**, taken from Supplementary Table 1 or the cBioPortal clinical data.
- **If no DNA-based tumor purity is available, the pipeline stops after step 02 and the user decides.** RNA-derived purity in the graph would be circular.

---

## 5. GRACKLE changes (`R/GRACKLE.R`)

The objective, the multiplicative updates, the inputs and the SVD-based initialization are unchanged. Changes:

1. **Loss bookkeeping matches the objective.**
   - Compute `||Y - WH||_F^2 + lambda_1 tr(W' L_p W) + lambda_2 tr(H L_g H')`.
   - This replaces the unsquared norm and the `reduce_sum` over all matrix entries in the current code (lines 85-96 before the loop, 123-139 inside it).
2. **`iterations` is honored** as the maximum number of iterations, replacing the hard-coded `1:292`. The analysis uses `iterations = 2000` and `diff_threshold = 1e-6` (relative change in the objective).
3. **Convergence is reported.** The return list gains `loss_trace` (a data frame of iteration, reconstruction, sample-graph term, gene-graph term and total), `n_iter` and `converged`.
4. **Precision.** A `dtype` argument defaults to `"float32"`. `"float16"` stays available. The global Keras `mixed_float16` policy side effect and the NaN/Inf float32 fallbacks are removed.
5. **Fail-fast.** The `while(!success) tryCatch(...)` retry loops are removed, so errors stop the job.
6. **`svdInit(Y, k)` helper.** It reproduces the script logic (absolute values of the SVD factors, min-max scaled per factor) so every fold initializes the same way from its own training data.

Return fields `W`, `H`, `svd_W` and `svd_H` are kept. **The HTP and poster results remain tied to `1395a5d` on `alpine`.** Reruns on this branch may stop at different iterations.

---

## 6. Sample graph (`R/sample_graph.R`)

`buildSampleGraph(samples, weights = c(pair = 1, annot = 1, purity = 1), k_nn = 10)` returns the n x n similarity `S` plus a per-part summary (edges, mean degree, sparsity).

```
S = w_pair * S_pair + w_annot * S_annot + w_pur * S_pur
```

**S_pair.**
- Weight 1 between each tumor sample and each of its derived 3D organoids, made symmetric.
- No edges between sibling organoids.

**S_annot.**
- Similarity is `1 - Gower distance` (`cluster::daisy(metric = "gower")`) over cancer type, tumor descriptor, prior treatment, age, sex and ancestry.
- Organoids inherit their parent tumor's values. Missing values are skipped within Gower.
- Each sample keeps its `k_nn` most similar neighbors. Ties are broken by absolute age difference, then by `profile_id`.
- Symmetrized with `pmax(S, t(S))`.

**S_pur.**
- Tumor-tumor edges only, with Gaussian similarity `exp(-(pi_i - pi_j)^2 / (2 sigma^2))`.
- `sigma` is the median distance from each tumor to its `k_nn`-th nearest purity neighbor.
- kNN sparsified and symmetrized as above. Tumors with missing purity get no purity edges, and organoids never have purity edges.

**Scaling.**
- Each part is divided by its mean degree, computed over nodes with at least one edge in that part. Then the weights 1 : 1 : 1 are applied.
- The weights are **fixed in advance and not tuned.**

**Per fold.** The graph is built from training-fold samples only. GRACKLE then forms `L_p = D - S` exactly as now.

---

## 7. Rank and tuning

### 7.1 Rank k (`R/rank_selection.R`, `04_rank.R`)
- **Method:** Horn's parallel analysis on the uncentered scaled matrix of all included samples.
  - Compute the singular values with `irlba` (`nv = 200`, increased if the estimate reaches the cap).
  - Make 100 null matrices by permuting each gene column independently, with seed `20261015`.
  - **k = the number of leading components whose singular value exceeds the 95th percentile of the null singular values at the same index.**
- k is computed once and used by every fold.

### 7.2 Lambda from alpha (`R/tuning.R`)
At each fit's SVD initialization `(W0, H0)`:

```
lambda_1 = alpha_1 * ||Y - W0 H0||_F^2 / tr(W0' L_p W0)
lambda_2 = alpha_2 * ||Y - W0 H0||_F^2 / tr(H0 L_g H0')
```

The grid is `alpha_1, alpha_2 in {0, 0.01, 0.03, 0.1, 0.3, 1}`, giving 36 settings.

### 7.3 Selection rule (`R/tuning.R`, applied in `06_combine.R`)
1. **Exclude non-converged settings.** Any setting with a non-converged fold is excluded and listed.
2. **Held-out reconstruction:** for each setting, the mean over folds of the squared Frobenius error of the held-out samples' NNLS projections, divided by the number of entries. SE = SD across the 5 folds / sqrt(5).
3. **Reproducibility:** for each setting, the cosine of every Hungarian match (Section 8.2) between fold programs and that setting's reference fit, averaged over reference programs and folds. Matches below `tau` are included at their actual cosine, so the measure stays continuous.
4. **Rule:** let `m*` be the setting with the lowest held-out reconstruction and `SE*` its SE. Among settings with reconstruction <= `m* + SE*`, choose the highest reproducibility. Exact ties go to the smaller `alpha_1`, then the smaller `alpha_2`.

The persistence score is never used in selection.

### 7.4 Sensitivity
- The persistence map is computed for all 36 settings (supplementary figure).
- At the selected alpha, the full pipeline also reruns at `round(k / 2)` and `2k` (supplementary figure).

---

## 8. Cross-fitting, matching and scoring

### 8.1 Cross-fitting (`05_fit.R`)
For each setting and fold f:
1. **Training set:** all included samples not in fold f. Build `S` (Section 6), set lambda (Section 7.2), and run `svdInit` on the training Y.
2. **Fit and normalize:** fit GRACKLE (`dtype = "float32"`, `iterations = 2000`, `diff_threshold = 1e-6`). Normalize each row of H to unit L2 norm.
3. **Project:** project fold f's tumors and organoids **separately** onto H with `project_W` (NNLS).

**Reference fit.** Each setting also gets one reference fit on all included samples. It is used only for program identity (names, gene loadings, enrichment, figures) and **never contributes a score**.

### 8.2 Program matching (`R/program_matching.R`)
- **Matching:**
  - Compute the k x k cosine similarity between the reference H rows and fold f's H rows.
  - Match one-to-one by maximum total similarity using `clue::solve_LSAP(sim, maximum = TRUE)`, the Hungarian algorithm.
- **Null threshold:**
  - For each fold program, permute its gene loadings 1,000 times and take the best cosine against all reference programs.
  - The threshold `tau` is the 99th percentile of these best-match values, pooled within a setting.
  - A match with cosine below `tau` counts as **unmatched**.
- **Stable program:** a reference program matched in at least 4 of 5 folds. Unstable programs go to the supplement and are excluded from the main map.

### 8.3 Persistence scores (`R/persistence.R`)
For each stable program j, the inputs are the held-out pairs (t, o) from every fold, each using the fold program matched to j.

**Concordance.**
- **Statistic:** rank tumor activities and organoid activities separately within cancer type (pooled across folds), center the ranks, and take the Pearson correlation of the pooled centered ranks. This is a stratified Spearman correlation.
- **Null:** 10,000 permutations that shuffle organoids among tumors **within cancer type x fold strata**.
- **P-value:** one-sided, `p = (1 + #{null >= obs}) / (1 + 10000)`.
- **Confidence interval:** 95% percentile, from 2,000 bootstrap resamples of pairs within strata.
- **Call:** concordant if BH FDR < 0.05 across stable programs.

**Retention.**
- **Statistic:** per pair, `log2((o + eps) / (t + eps))`, where `eps` is the 1st percentile of program j's nonzero held-out activities, pooling tumors and organoids. Retention is the median over pairs.
- **Test:** 10,000 random sign flips within pairs, with BH FDR across stable programs. CI from the same bootstrap.

**Groups** (fixed in advance and mutually exclusive; the 2-fold threshold is a declared convention). Call a program **reduced** if retention is significantly below 0 (FDR < 0.05) and the median log2 ratio is <= -1. Then:
- **Drops out:** reduced and not concordant.
- **Reduced but concordant:** reduced and concordant.
- **Persists:** not reduced and concordant.
- **Present, not patient-specific:** not reduced and not concordant.

**Detection limit.** Each program reports the 95th percentile of its concordance null. The synthetic fixture supplies the power curve.

**Per-type view.** Concordance is computed within each cancer type that has at least 10 held-out pairs. Smaller types contribute to the pooled statistic but are not shown separately.

### 8.4 Synthetic fixture expectations
The fixture has about 60 pairs, 3 cancer types, 500 genes and a random gene graph. Its planted programs must land in these groups:

| Planted program | Expected group |
|---|---|
| Two patient-specific programs shared by tumor and organoid | Persists |
| One tumor-only microenvironment program, scaled by 1 - purity | Drops out |
| One lineage-only program (constant within type, plus noise) | Present, not patient-specific |

---

## 9. Annotation (`R/annotation.R`, from the reference fit)

- **Enrichment:**
  - Take the top 5% of genes by H loading.
  - Test against MSigDB Hallmark (pinned `msigdbr`) with a hypergeometric test (`phyper`), using the final gene set as background.
  - BH across program x gene set, FDR < 0.05.
- **Microenvironment label:**
  - Spearman correlation of out-of-fold tumor-side activity with DNA purity.
  - A program is labeled microenvironment-associated if rho < 0 and BH FDR < 0.05.
  - ESTIMATE stromal and immune scores (`tidyestimate`, computed on tumor log2 TPM) are reported as supporting columns.
  - **Positive control:** microenvironment-associated programs are expected to drop out.
- **Genotype:**
  - Drivers mutated in at least 10 included tumors.
  - Wilcoxon test of tumor-side activity, mutated vs. wild type, with BH across program x gene.
  - Annotation only.

---

## 10. Outputs

On Alpine, outputs go to scratch. They are pulled to `analysis/hcmi_persistence/results/`, which is gitignored and committed with `git add -f` only after the user has reviewed them.

- `programs.csv`: one row per reference program:
  - stability: folds matched, mean matched cosine
  - concordance with CI, p and FDR
  - retention with CI and FDR, plus the group label
  - detection limit
  - purity rho/FDR and ESTIMATE correlations
  - top enrichments and genotype associations
- `selection.csv`: every setting's reconstruction (mean and SE), reproducibility, convergence and the selected flag.
- `convergence.csv`: one row per fit with `n_iter`, `converged` and the final loss components.
- `k_parallel_analysis.csv`, plus a figure of observed vs. null singular values.
- **Main figure:** retention (x) vs. concordance (y), one point per stable program. Points are colored by microenvironment vs. tumor-intrinsic and labeled with their top enrichment. A companion heatmap shows concordance per program x cancer type.
- **Supplementary figures:** the alpha grid; k/2 and 2k; unstable programs; the parallel-analysis plot.
- **Provenance:** every output has a sidecar `.provenance.json` recording the commit, Slurm job ID, setting, input checksums and `sessionInfo()`.

---

## 11. Compute on Alpine

These conventions follow fm-pm-evaluator (`scripts/alpine/ralpine`, `docs/PROCESS.md` section 2).

- **Access**
  - Claude uses only `scripts/alpine/ralpine`, copied from fm-pm-evaluator with root `/projects/lgillenwater@xsede.org/repositories/GRACKLE-organoid-transfer`. Never raw ssh.
  - The user opens the Duo session once per working session (`! ssh -fN alpine`).
- **Checkout**
  - The existing `/projects/$USER/repositories/GRACKLE` (a detached HEAD with uncommitted HTP work) is **not touched**.
  - Once the branch is pushed, the user runs this once: `git -C /projects/$USER/repositories/GRACKLE fetch origin && git -C /projects/$USER/repositories/GRACKLE worktree add ../GRACKLE-organoid-transfer organoid-transfer`.
- **Deployment** is git-only: `git push`, then `ralpine update`, then `ralpine submit`. `submit` refuses when Alpine's HEAD differs from local HEAD.
- **Environment**
  - `setup_env.sbatch` fetches a standalone micromamba binary and builds `/projects/$USER/software/envs/grackle` from the committed `environment.yml` (conda-forge: `r-base`, `r-reticulate`, `r-tensorflow`, CPU `tensorflow`, `r-igraph`, `r-nnls`, `r-cluster`, `r-clue`, `r-irlba`, `r-msigdbr`, `r-tidyestimate`, `r-devtools`, `r-testthat`, `r-jsonlite`, `r-httr2`, `r-tidyverse`).
  - It writes `environment.lock.txt` (an explicit export), which is committed.
  - `scripts/alpine/env.sh` puts the environment's `bin/` on `PATH`; `module load anaconda` fails inside batch jobs.
- **Jobs**
  - Partition `acpu`, QoS `cpu-normal`. Memory equals cores x 3,840 MB, so a fit task uses 4 cores and 15G.
  - `set -euo pipefail`. Each job prints its commit and job ID, and logs go to `logs/%x-%j.out`.
- **Job chain** (`submit_chain.sh`, with each dependency verified by `ralpine jobinfo`):
  1. `setup_env`
  2. `01_download`
  3. `02_sample_table`
  4. `03_expression`
  5. `04_rank`
  6. `fit_array`: rows of `settings.csv` = 36 settings x (5 folds + 1 reference) = 216 tasks. Each task skips if its output exists.
  7. `06_combine`
  8. `sensitivity_array`: k/2 and 2k x 6 fits = 12 tasks.
  9. `06_combine --sensitivity`
  10. `07_figures`
- **Storage**
  - Raw data and intermediates: `/projects/$USER/grackle_data/hcmi/`
  - Run outputs: `/scratch/alpine/$USER/grackle/hcmi_persistence/<run_id>/`, where `run_id` is `<UTC date YYYYMMDD>_<short commit>`, set once by `submit_chain.sh` and passed to every job
  - Logs are pulled with `ralpine pull`.
- **Standing permissions** (approved by the user, effective once this spec is approved)
  - Commit and push to `origin` on `organoid-transfer`; `ralpine update` / `switch` / `submit` / `cancel <one id>`.
  - Ask first for a force-push, `reset --hard`, branch deletion, or anything touching `main` or `alpine`.

---

## 12. Code layout

```
R/
  GRACKLE.R              # Section 5 changes, svdInit()
  sample_graph.R         # buildSampleGraph() and its parts
  rank_selection.R       # parallelAnalysisRank()
  tuning.R               # lambdaFromAlpha(), selectSetting()
  folds.R                # makeFolds()
  program_matching.R     # matchPrograms(), matchNull(), programStability()
  persistence.R          # concordance, retention, nulls, bootstrap, detection limit, groups
  annotation.R           # enrichment, purity/ESTIMATE, genotype
  Evaluation.R           # existing project_W (reused)
analysis/hcmi_persistence/
  config.R, 01_download.R ... 07_figures.R
  settings.csv, sample_table.csv, folds.csv, filter_counts.csv
  results/               # gitignored until reviewed
scripts/alpine/
  ralpine, env.sh, setup_env.sbatch, <step>.sbatch, fit_array.sbatch,
  sensitivity_array.sbatch, submit_chain.sh
environment.yml, environment.lock.txt
tests/testthat.R, tests/testthat/
```

New exported functions get roxygen documentation, and `NAMESPACE` and `man/` are regenerated with `devtools::document()`. New package dependencies (`cluster`, `clue`, `irlba`, `msigdbr`, `tidyestimate`, `jsonlite`, `httr2`) are added to `DESCRIPTION`.

---

## 13. Testing

**Unit tests** (`testthat`, run with `devtools::test()`):

| Unit | Tests |
|---|---|
| `GRACKLE` | Objective doesn't increase across iterations, with graph terms off and on; `iterations` honored; `dtype` honored; `loss_trace` columns present; `lambda_1 = lambda_2 = 0` equals plain multiplicative-update NMF from the same initialization |
| `buildSampleGraph` | Pair edges only between a tumor and its own organoids; no purity edge touches an organoid; symmetric; deterministic under ties; mean degree 1 per part before weighting; held-out IDs absent |
| `makeFolds` | No case split across folds; stratified by type; reproducible with the seed |
| `parallelAnalysisRank` | Recovers a planted rank r from low-rank data plus noise |
| `lambdaFromAlpha` | Penalty term equals `alpha` x reconstruction at the starting point |
| `selectSetting` | Applies the one-SE rule, reproducibility and tie-break correctly on a constructed table |
| `matchPrograms` | Recovers a known permutation; solves the 2x2 greedy-failure example optimally; a below-threshold match counts as unmatched |
| Persistence | Permutation p-values are roughly uniform under a no-signal null; shuffles never cross a cancer type x fold stratum; sign-flip retention test is calibrated under the null |
| `ralpine` | Text-scan boundary tests ported from fm-pm-evaluator's `test_ralpine_boundary.py` |

**Synthetic end-to-end fixture.** This is required before any Alpine submission.
- The real `05_fit.R` and `06_combine.R` code runs on the Section 8.4 fixture with a 2 x 2 alpha grid.
- The planted programs must land in their expected groups.
- It runs locally if R and TensorFlow are available; otherwise it runs as an Alpine `atesting` job.

---

## 14. Risks and fallbacks

| Risk | Response |
|---|---|
| No DNA-based tumor purity | Stop after step 02; the user decides (Section 4.6) |
| Supplementary Table 1 not on the GDC publication page | The user places the Nature supplementary file in `raw/manual/`; its checksum is recorded |
| Model type or parent sample missing for some organoids | Unknown model type: excluded, with the reason recorded. Unknown parent: pair with all of the case's tumor samples and flag the case |
| micromamba blocked on Alpine | Use a pinned Apptainer image (`rocker/r-ver` at a fixed version plus the same packages), run with `apptainer exec` |
| Large k (> 100) | Raise `irlba` `nv`; fit cost grows roughly linearly in k and stays within the 4-core, 1-hour task budget at this size |
| Non-converged fits | Listed in `convergence.csv`; settings with any non-converged fold are excluded from selection |
| Most programs unstable | Reported as a result; the k/2 sensitivity run shows whether a smaller rank is reproducible |
| Small cancer types | Pooled statistic still valid; the per-type view is limited to types with at least 10 held-out pairs |

---

## 15. Verifications at implementation time

These are checks of facts, not open design choices. The response to each outcome is fixed above.

1. Whether Supplementary Table 1 is on the GDC publication page, and whether it has model type, parent sample and DNA purity (Sections 4.1, 4.6, 14).
2. The included organoid pair count after filtering, which goes into `filter_counts.csv` and corrects the grant text.
3. Whether R plus TensorFlow runs locally for the fixture (Section 13).
4. Whether micromamba runs on Alpine (Section 14).
