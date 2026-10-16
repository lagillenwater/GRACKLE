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
- Four gene-graph arms, each from a single source: co-essentiality 2D and 3D (co-primary), STRING mechanistic and STRING all channels (reference). All four run the identical pipeline.
- A tumor-to-organoid transfer test comparing the gene priors.
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
2. The synthetic fixture places each planted program in its correct group (Section 8.5) before any Alpine submission.
3. The main figure (co-essentiality 2D and 3D maps side by side, plus their agreement), the per-program tables and the supplementary panels exist for each arm's selected setting.
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
- **Data locations:** the GDC project `HCMI-CMDC` (open STAR counts) and the GDC publication page `HCMI-CMDC-2026` (which hosts only file manifests). Mutations and clinical data are also in cBioPortal `pancan_hcmi_2025`.

From the GDC API (2026-10-15):
- 1,153 open-access STAR-count RNA-seq files.
- 248 cases have both a tumor-type and a model-type RNA profile.
- GDC sample fields lack model type, medium and center. Passage count is present for about 370 model samples.

From STRING v12 `9606.protein.links.full.v12.0.txt` (local copy, 2026-10-15):
- **Edge shares by channel:** text mining 87%, coexpression 69%, experiments 43%, database 6%, genomic context 6%.
- **54% of edges have no experiments or database evidence**; they rest only on coexpression and/or text mining.
- 19,235 of 19,622 proteins (98%) have at least one experiments or database edge.

CRISPR dependency screens on 3D models:
- **DepMap NextGen:** 147 genome-scale screens on organoid and spheroid models across 10 cancer types (Neiswender et al. 2026).
- **Sanger:** 162 screens on HCMI organoids (Herranz-Ors et al. 2026).
- 3D growth changes dependencies relative to 2D (Han et al. 2020, *Nature*).

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
| What a gene prior must be | A property of genes measured independently of expression; the gene-side mirror of the sample-graph rule | Coexpression: a property of other samples, double-counts expression, and dominated in bulk tissue by stromal/immune composition. Text mining: study bias, and it partly encodes known drug-sensitivity genes (leakage risk for later aims) |
| Gene-graph arms | One source per arm: co-ess 2D and 3D co-primary; STRING mechanistic and STRING all channels as references | STRING mechanistic + co-essentiality combined: weakens the 2D/3D contrast, reintroduces curated-pathway circularity with Hallmark enrichment, and adds an arbitrary weight |
| Genes without gene-prior evidence | No edges, so no smoothing (plain NMF for that gene) | HTP convention of dropping genes with no STRING edge: biases the gene set toward well-studied proteins |
| Other functional sources | Not used in this piece | Drug co-response (Tahoe-100M, LINCS): expression of other samples, and the data Aim 2 evaluates on. Perturb-seq (K562/RPE1): few contexts, mostly essential genes. Genetic-interaction maps: too sparse |
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
- **cBioPortal `pancan_hcmi_2025`, via its REST API**: all patient and sample clinical attributes, plus mutations from profile `pancan_hcmi_2025_mutations` over sample list `pancan_hcmi_2025_sequenced`. Its sample IDs equal GDC sample submitter IDs (e.g., `HCM-CSHL-0085-C24-01A`). Checked 2026-10-15:
  - it carries `GROWTH_PATTERN` ("3D Organoid", "2D Cell Line", "3D Spheroid"), `CULTURE_MEDIA`, `TREATMENT_STATUS`, `DIAGNOSIS_TCGA_HCMI_CODE`, `AGE`, `SEX`, `RACE_CATEGORY`, and four DNA-based purity callers (`PURITY_ABSOULTE` [sic], `PURITY_ASCAT`, `PURITY_HACHET` [sic], `PURITY_PURPLE`) plus RNA-based `PURITY_ESTIMATE`;
  - the GDC publication page `HCMI-CMDC-2026` hosts only file manifests, so the paper's Supplementary Table 1 is not needed;
  - the cBioPortal datahub S3 archive returns 403, so the API is used.
- **STRING v12** (`9606.protein.links.full.v12.0.txt`, `9606.protein.links.v12.0.txt`, `9606.protein.aliases.v12.0.txt`). The full-links file carries the per-channel scores.
- **DepMap** (checked 2026-10-15). The DepMap portal blocks scripted downloads with a bot-verification page, but DepMap mirrors releases on figshare.
  - **Scripted from figshare, with md5 checks:**
    - DepMap 24Q4 Public (figshare+, CC BY 4.0, doi 10.25452/figshare.plus.27993248): `CRISPRGeneEffect.csv`, `Model.csv`, `AchillesNonessentialControls.csv`.
    - The NextGen manuscript metadata (CC BY 4.0, doi 10.6084/m9.figshare.29472362): `model_metadata.csv` and `screen_metadata.csv`, which give each screen's culture format (`CRISPRScreenType`).
  - **Manual, by the user:** the genome-wide NextGen gene effects exist only in the portal file set "NextGen Model Manuscript 2026". The user downloads it in a browser into `raw/manual/depmap_nextgen/`; `01_download.R` verifies it and records sha256.
  - **Which matrix each arm uses:** if the NextGen gene-effect matrix also covers the traditional 2D screens (one Chronos run), it supplies both co-essentiality arms, so the 2D/3D contrast is not confounded by processing. Otherwise the 2D arm uses 24Q4.
- **Sanger HCMI organoid CRISPR screens** (Herranz-Ors et al. 2026), from their deposited gene-level scores. Downloaded by script if a direct URL exists; otherwise by the user into `raw/manual/sanger/`, the same way as DepMap.
- **OncoKB cancer gene list** (`https://www.oncokb.org/api/v1/utils/cancerGeneList.txt`, public, checksummed). It defines which genes count as drivers for the genotype annotation.
- **Gene annotation:** HGNC complete set, the pinned 2024-10-01 file `ensemblToHGNC` already uses. It supplies symbols, Ensembl and Entrez IDs, and the cytogenetic `location` used for chromosome arms.
- **MSigDB Hallmark** via `msigdbr`, with the package version pinned and the gene sets cached as an RDS file.

### 4.2 Sample table (`02_sample_table.R`)
- **One row per RNA profile.** Columns:
  - identity: `profile_id`, `aliquot_id`, `sample_id`, `case_id`, `role` (tumor / organoid)
  - tumor context: `cancer_type`, `tumor_descriptor`
  - model: `model_id`, `model_type`, `medium`, `passage`, `derivation_center`
  - pairing: `parent_sample_id`, `pair_id`
  - covariates: `purity_dna`, `prior_treatment`, `age`, `sex`, `ancestry`
  - bookkeeping: `fold`, `included`, `exclusion_reason`
- **Field sources:**
  - **Model type:** cBioPortal `GROWTH_PATTERN` of the model sample.
  - **Medium:** `CULTURE_MEDIA`.
  - **Passage:** GDC `passage_count`.
  - **Derivation center:** the center code in the HCMI ID (`HCM-BROD-...` gives `BROD`).
  - **Cancer type:** `DIAGNOSIS_TCGA_HCMI_CODE`, set per case as the most frequent code among its samples. When a case has none, the ICD-O site in its HCMI id is used (for example `site:C18`). One value per case is required by the folds.
  - **Prior treatment:** `TREATMENT_STATUS == "Treated"`.
  - **Age and sex:** `AGE`, `SEX`.
  - **Ancestry annotation:** `RACE_CATEGORY` (reported race category, not the paper's genetic ancestry estimate).
- **Included organoids:** model profiles whose `GROWTH_PATTERN` is "3D Organoid".
- **Included tumors:** all tumor profiles with RNA. That covers primary, metastatic and recurrent tumors, human tumor original cells, and failed-derivation tumors.
- **FFPE:** a profile from an FFPE-preserved sample (`preservation_method == "FFPE"`, including slides and FFPE scrolls) is included only when no frozen profile exists for the same tumor sample. "Same tumor sample" means the same case and the same `tumor_descriptor`, because slides and scrolls get their own sample IDs. The paper attributes some low concordance to FFPE artifacts. The preservation method is kept as a column.
- **Duplicate profiles:** when one sample has more than one RNA profile, the profile with the lexicographically first `profile_id` is kept and the rest are excluded ("duplicate profile of the same sample"). This keeps a patient from entering the pairs twice through re-sequencing.
- **Unpaired organoids:** 3D organoids with no included tumor partner are excluded, following the training-set decision in Section 3.
- **Pairing:**
  - An organoid pairs with the tumor sample it was derived from. That is the case's one included tumor profile whose GDC `tumor_descriptor` matches the model sample's `tumor_descriptor` (the origin of the model).
  - If no single tumor sample matches, it pairs with every tumor sample of the case, and the case is flagged.
  - ATCC re-expanded models ("Expanded Next Generation Cancer Model") are excluded as duplicates.
- **Filter-count table:** the script writes one row per filter step with profiles, cases and pairs remaining (`filter_counts.csv`).
- **Committed:** `sample_table.csv` and `filter_counts.csv`, which are small and open-access.

### 4.3 Folds (`02_sample_table.R`)
- 5 folds with seed `20261015`, stratified by cancer type and grouped by case. Every case gets a fold, including unpaired tumors.
- Committed as `folds.csv`.

### 4.4 Expression (`03_expression.R`)
1. Take `tpm_unstranded` from each STAR-count file and compute log2(TPM + 1). Keep `gene_type == "protein_coding"`.
2. Select the top 5,000 genes by variance across all included samples. **Genes are not dropped for lacking gene-graph edges** (this departs from the HTP convention); a gene with no edges in an arm is simply not smoothed in that arm.
3. Min-max scale each gene across all included samples (`min_max_scale`).

Gene selection and scaling use all samples but no pairing or annotation information. The methods section states this.

### 4.5 Gene graphs (`R/gene_graph.R`, `03b_gene_graphs.R`)

Four arms, each a single source, all over the final gene set:

| Arm | Source | Role |
|---|---|---|
| `coess_2d` | DepMap Chronos, 2D models, HCMI-derived models excluded | Co-primary |
| `coess_3d` | Organoid and spheroid screens (DepMap NextGen + Sanger HCMI), rebuilt per fold | Co-primary |
| `string_mech` | STRING experiments + database channels, recombined | Reference |
| `string_all` | STRING `combined_score` (the HTP graph) | Reference |

**Common rules.**
- Weights are nonnegative and symmetric, with a zero diagonal. Self-loops cancel in the Laplacian and don't change the objective.
- Gene identifiers are mapped to HGNC symbols with the pinned HGNC file.
- **Graph scale needs no normalization**, because `lambda_2` is set from `alpha_2` relative to `tr(H0 L_g H0')` (Section 7.2).

**`coess_2d` and `coess_3d`.**
1. **Models.**
   - 2D: traditional screens (`CRISPRScreenType == "2DS"`; in 24Q4, growth pattern Adherent, Suspension or Mixed), excluding every HCMI-derived model.
   - 3D: NextGen screens run in 3D culture, `CRISPRScreenType` 3DO (dome organoids, 94) or 3DN (CNS spheroids, 14). This uses the screen's culture format, not the model's derivation format. NextGen models screened in 2D coat culture (2DO, 2DN) are in neither arm. If DepMap NextGen and Sanger scores cannot be put on one scale (checked at implementation), use the larger source alone. For each fold, the screens of models from that fold's held-out cases are removed before the graph is built. The reference fit uses all 3D screens.
2. **Informative genes.** A gene gets edges only if the SD of its Chronos scores across the arm's models exceeds the 95th percentile of SDs among DepMap's non-essential control genes. Without this filter, the nearest-neighbor step would give noise-only genes random edges.
3. **Nuisance correction.** GLS whitening across models (Wainberg et al. 2021, *Nat Genet*): estimate the model-by-model covariance from all genes and whiten the model dimension by its inverse square root. This removes shared line-level structure such as screen quality and proliferation.
4. **Similarity.** GLS correlation between informative genes: the cosine of their whitened profiles after projecting out the whitened intercept (`1' * Sigma^(-1/2)`). This is the GLS analogue of a Pearson correlation; plain centering after whitening would leave each gene's mean effect as a shared pattern and correlate strongly essential genes with each other. Pairs on the same chromosome arm are set to 0 to remove CRISPR proximity artifacts; the arm comes from the cytogenetic `location` field of the pinned HGNC file (e.g., `17q21.31` gives `17q`). Negative values are set to 0.
5. **Sparsify.** Keep each gene's 10 nearest neighbors among the final genes, with weight equal to the correlation, symmetrized with `pmax(S, t(S))`.

**`string_mech`.**
- Recombine the `experiments`, `experiments_transferred`, `database` and `database_transferred` channels with STRING's published channel-combination procedure:
  1. For each channel, `s' = max(0, (s / 1000 - 0.041) / (1 - 0.041))`.
  2. Combine: `c = 1 - prod(1 - s')`.
  3. Add the prior back: `c * (1 - 0.041) + 0.041`.
- **Unit-tested:** recombining all channels with the same procedure reproduces `combined_score` (to within 1 on the 0-1000 scale). For this test, direct and transferred scores are combined per channel. No homology correction is applied, because the published channel scores already carry it (checked on the v12 file).
- Weight = recombined score; no threshold.

**`string_all`.**
- `combined_score / 1000` with no threshold, as in `tests/HTP_v2.R`. The diagonal is 0 rather than 1, which leaves the objective unchanged.

### 4.6 Purity requirement
- Graph purity must be **DNA-based**: `purity_dna` is the median of the available DNA callers (`PURITY_ABSOULTE`, `PURITY_ASCAT`, `PURITY_HACHET`, `PURITY_PURPLE`) for each tumor sample. cBioPortal's `CONSENSUS_PURITY` is not a purity estimate in this study: it reads about 0.1-0.35 where all four DNA callers agree near 0.8-0.9, with r = -0.31 against ABSOLUTE (checked 2026-10-15). So agreement among the four DNA callers (median pairwise r) is reported instead.
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

The persistence score is never used in selection. Selection runs separately within each gene-graph arm.

### 7.4 Sensitivity
For each arm:
- The persistence map is computed for all 36 settings (supplementary figure).
- At the selected alpha, the full pipeline also reruns at `round(k / 2)` and `2k` (supplementary figure).

---

## 8. Cross-fitting, matching and scoring

### 8.1 Cross-fitting (`05_fit.R`)
Sections 8.1-8.3 run identically within each of the four gene-graph arms (Section 4.5). For each arm, setting and fold f:
1. **Training set:** all included samples not in fold f. Build `S` (Section 6), take the arm's gene graph (for `coess_3d`, fold f's version), set lambda (Section 7.2), and run `svdInit` on the training Y.
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
- **Statistic:** rank tumor activities and organoid activities separately within cancer type (pooled across folds). Convert to normalized centered ranks, `(rank - (n + 1) / 2) / n` with n the type's pair count, so each type is weighted by its number of pairs. Take the Pearson correlation of the pooled values. This is a stratified Spearman correlation. If either side has zero variance (for example, all NNLS activities are 0), the statistic is 0 and p = 1.
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

**Detection limit and power.** Each program reports the 95th percentile of its concordance null. The power curve of the concordance test comes from a simulation at the run's real held-out pair structure (cancer type x fold counts), giving power against true concordance (`compare/power_curve.csv`). The synthetic fixture checks that planted programs are recovered.

**Per-type view.** Concordance is computed within each cancer type that has at least 10 held-out pairs. Smaller types contribute to the pooled statistic but are not shown separately.

### 8.4 Comparing gene priors (`06_combine.R --compare`)

**Transfer test.** For each arm, at two settings: its selected `(alpha_1, alpha_2)`, and the same `alpha_1` with `alpha_2 = 0` (that arm's own no-gene-prior baseline):
1. For each fold f, fit GRACKLE on fold f's training **tumors only**. The sample graph is restricted to tumors, so there are no pair edges.
2. Project fold f's held-out **organoids** with NNLS.
3. Score by mean reconstruction error over folds, with SE. The gene prior's effect in an arm is the paired per-fold difference between its two settings.

This measures how well programs learned from tumors describe organoids, and it uses no pairing information. It is reported, not used for selection.

**Map agreement between `coess_2d` and `coess_3d`.**
- Match the two arms' reference programs with the Hungarian method on H cosine, using the Section 8.2 threshold.
- For matched programs, report agreement of group calls (percent agreement and Cohen's kappa) and the correlation of concordance values.
- List unmatched programs on each side; programs found under only one prior are a result in themselves.

### 8.5 Synthetic fixture expectations
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
  - In the STRING arms, database edges come partly from the same curated pathways behind Hallmark sets, so enrichment there is partly built in by the graph. The co-essentiality arms have no such overlap. The tables flag which arm each enrichment comes from.
- **Microenvironment label:**
  - Spearman correlation of out-of-fold tumor-side activity with DNA purity.
  - A program is labeled microenvironment-associated if rho < 0 and BH FDR < 0.05.
  - ESTIMATE stromal and immune scores (`tidyestimate`, computed on tumor log2 TPM) are reported as supporting columns.
  - **Positive control:** microenvironment-associated programs are expected to drop out.
- **Genotype:**
  - Drivers are OncoKB cancer genes. A tumor counts as mutated for a gene if it has a non-silent mutation (missense, nonsense, frameshift, in-frame indel, splice site, nonstop or translation-start). Only genes mutated in at least 10 included tumors, with at least 2 wild-type tumors, are tested.
  - Wilcoxon test of tumor-side activity, mutated vs. wild type, with BH across program x gene.
  - Annotation only.

---

## 10. Outputs

On Alpine, outputs go to scratch. They are pulled to `analysis/hcmi_persistence/results/`, which is gitignored and committed with `git add -f` only after the user has reviewed them.

Per-arm files go under `results/<arm>/`. Comparison files go under `results/compare/`.

- `<arm>/programs.csv`: one row per reference program:
  - stability: folds matched, mean matched cosine
  - concordance with CI, p and FDR
  - retention with CI and FDR, plus the group label
  - detection limit
  - purity rho/FDR and ESTIMATE correlations
  - top enrichments and genotype associations
- `<arm>/selection.csv`: every setting's reconstruction (mean and SE), reproducibility, convergence and the selected flag.
- `<arm>/convergence.csv`: one row per fit with `n_iter`, `converged` and the final loss components.
- `gene_graph_summary.csv`: per arm (and per fold for `coess_3d`), the number of models used, informative genes, genes with at least one edge among the final genes, edges, and mean degree.
- `k_parallel_analysis.csv`, plus a figure of observed vs. null singular values.
- `compare/transfer.csv`: tumor-to-organoid transfer reconstruction (mean and SE) per arm at its selected setting and at its `alpha_2 = 0` baseline, plus the paired difference.
- `compare/map_agreement.csv`: matched programs between `coess_2d` and `coess_3d`, with group-call agreement, kappa and concordance correlation.
- **Main figure:**
  - Two panels, `coess_2d` and `coess_3d`, each plotting retention (x) against concordance (y), one point per stable program. Points are colored by microenvironment vs. tumor-intrinsic and labeled with their top enrichment.
  - A third panel shows agreement between the two maps.
  - Companion heatmaps show concordance per program x cancer type for each co-primary arm.
- **Supplementary figures:**
  - The STRING arms' maps.
  - The transfer test.
  - For each arm: the alpha grid, k/2 and 2k, and unstable programs.
  - The parallel-analysis plot.
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
  - `setup_env.sbatch` fetches a standalone micromamba binary and builds `/projects/$USER/software/envs/grackle` from the committed `environment.yml` (conda-forge: `r-base`, `r-reticulate`, `r-tensorflow`, `tensorflow-cpu`, `r-igraph`, `r-nnls`, `r-cluster`, `r-clue`, `r-irlba`, `r-data.table`, `r-msigdbr`, `r-glmnet`, `r-fnn`, `r-ppcor`, `r-devtools`, `r-testthat`, `r-jsonlite`, `r-httr2`, `r-tidyverse`). `tidyestimate` is not on conda-forge, so it is installed from CRAN at a pinned version (1.1.1) inside the environment.
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
  5. `03b_gene_graphs`: all four arms, with `coess_3d` built per fold.
  6. `04_rank`
  7. `fit_array`: **one array per arm.** Each array covers the arm's rows of `settings.csv`, 36 settings x (5 folds + 1 reference) = 216 tasks, throttled with `%50`. Each task skips if its output exists.
  8. `06_combine`, per arm.
  9. `sensitivity_array`: per arm, k/2 and 2k x 6 fits = 12 tasks (48 in all).
  10. `transfer_array`: 4 arms x 2 settings x 5 folds = 40 tasks.
  11. `06_combine --sensitivity`, then `06_combine --compare`.
  12. `07_figures`.
  - **Total about 952 fits.** At an estimated 1-3 minutes each on 4 cores, that is roughly 65-190 core-hours.
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
  gene_graph.R           # stringRecombine(), stringGraph(), informativeGenes(), glsWhiten(), coessGraph()
  rank_selection.R       # parallelAnalysisRank()
  tuning.R               # lambdaFromAlpha(), selectSetting()
  folds.R                # makeFolds()
  program_matching.R     # matchPrograms(), matchNull(), programStability()
  persistence.R          # concordance, retention, nulls, bootstrap, detection limit, groups, power curve
  annotation.R           # enrichment, purity/ESTIMATE, genotype
  cross_fitting.R        # fitTask(): the core of 05_fit.R
  combine.R              # combineArm(): the core of 06_combine.R
  compare.R              # transfer table, 2D/3D map agreement
  fixture.R              # synthetic fixture and its expected calls
  hcmi_data.R            # GDC STAR parsing, buildSampleTable()
  depmap_data.R          # DepMap/Sanger parsing, 3D source combination
  Evaluation.R           # existing project_W (reused)
analysis/hcmi_persistence/
  config.R, 01_download.R, 02_sample_table.R, 03_expression.R, 03b_gene_graphs.R,
  04_rank.R, 05_fit.R, 06_combine.R, 07_figures.R
  settings.csv, sample_table.csv, folds.csv, filter_counts.csv
  results/               # gitignored until reviewed
scripts/alpine/
  ralpine, env.sh, setup_env.sbatch, <step>.sbatch, fit_array.sbatch,
  sensitivity_array.sbatch, transfer_array.sbatch, submit_chain.sh
environment.yml, environment.lock.txt
tests/testthat.R, tests/testthat/
```

New exported functions get roxygen documentation, and `NAMESPACE` and `man/` are regenerated with `devtools::document()`.

**Package loadability fix**, required for any clean environment:
- `PLIER`, `sgnesR` and `NMF` move from `@import` to `Suggests`, called with `pkg::` behind `requireNamespace()` checks. `PLIER` and `sgnesR` are GitHub-only, and `NMF` pulls in Bioconductor.
- The unused `@import dynutils` is removed.
- `@useDynLib GRACKLE`, `LinkingTo: Rcpp` and `R/RcppExports.R` are removed. That file calls compiled code, and there is no `src/` directory.

New `Imports`: `tensorflow`, `nnls`, `cluster`, `clue`, `irlba`, `data.table`, `jsonlite`, `withr` (seeded RNG without touching the caller's state). New `Suggests`: `tidyestimate`, `testthat (>= 3.0.0)`.

---

## 13. Testing

**Unit tests** (`testthat`, run with `devtools::test()`):

| Unit | Tests |
|---|---|
| `GRACKLE` | Objective doesn't increase across iterations, with graph terms off and on; `iterations` honored; `dtype` honored; `loss_trace` columns present; `lambda_1 = lambda_2 = 0` equals plain multiplicative-update NMF from the same initialization |
| `buildSampleGraph` | Pair edges only between a tumor and its own organoids; no purity edge touches an organoid; symmetric; deterministic under ties; mean degree 1 per part before weighting; held-out IDs absent |
| `makeFolds` | No case split across folds; stratified by type; reproducible with the seed |
| `stringRecombine` | Recombining all channels reproduces `combined_score` to within 1; the mechanistic recombination never exceeds `combined_score` |
| `glsWhiten` | Whitened model-by-model covariance is close to the identity on synthetic data |
| `informativeGenes` | Keeps planted variable-dependency genes and drops planted noise-only genes |
| `coessGraph` | Nonnegative, symmetric, zero diagonal; no same-arm edges; edges only among informative genes; the per-fold 3D build excludes held-out cases' models |
| `parallelAnalysisRank` | Recovers a planted rank r from low-rank data plus noise |
| `lambdaFromAlpha` | Penalty term equals `alpha` x reconstruction at the starting point |
| `selectSetting` | Applies the one-SE rule, reproducibility and tie-break correctly on a constructed table |
| `matchPrograms` | Recovers a known permutation; solves the 2x2 greedy-failure example optimally; a below-threshold match counts as unmatched |
| Persistence | Permutation p-values are roughly uniform under a no-signal null; shuffles never cross a cancer type x fold stratum; sign-flip retention test is calibrated under the null |
| `ralpine` | Text-scan boundary tests ported from fm-pm-evaluator's `test_ralpine_boundary.py` |

**Synthetic end-to-end fixture.** This is required before any Alpine submission.
- The real `05_fit.R` and `06_combine.R` code runs on the Section 8.5 fixture with a 2 x 2 alpha grid and a random gene graph standing in for one arm.
- The planted programs must land in their expected groups.
- It runs locally if R and TensorFlow are available; otherwise it runs as an Alpine `atesting` job.

---

## 14. Risks and fallbacks

| Risk | Response |
|---|---|
| No DNA-based tumor purity | Stop after step 02; the user decides (Section 4.6) |
| cBioPortal lacks growth pattern or purity for many included samples | Report the coverage after step 02. Missing growth pattern excludes the model ("model type unknown"); missing purity gives that tumor no purity edges |
| Model type or parent sample missing for some organoids | Unknown model type: excluded, with the reason recorded. Unknown parent: pair with all of the case's tumor samples and flag the case |
| micromamba blocked on Alpine | Use a pinned Apptainer image (`rocker/r-ver` at a fixed version plus the same packages), run with `apptainer exec` |
| Large k (> 100) | Raise `irlba` `nv`; fit cost grows roughly linearly in k and stays within the 4-core, 1-hour task budget at this size |
| Non-converged fits | Listed in `convergence.csv`; settings with any non-converged fold are excluded from selection |
| Most programs unstable | Reported as a result; the k/2 sensitivity run shows whether a smaller rank is reproducible |
| Small cancer types | Pooled statistic still valid; the per-type view is limited to types with at least 10 held-out pairs |
| DepMap or Sanger files can't be downloaded by script | The user downloads them in a browser and runs one `rsync` to `raw/manual/`; `01_download.R` verifies them and records checksums |
| DepMap NextGen and Sanger 3D screens can't be put on one scale | `coess_3d` uses the larger source alone; recorded in `gene_graph_summary.csv` |
| Few informative genes among the final 5,000 | Reported in `gene_graph_summary.csv`. Genes without edges are unsmoothed by design, and graph sparsity is part of the arm comparison |
| HCMI models not cleanly flagged in DepMap metadata | Match HCMI model IDs (`HCM-*`) against DepMap model names and aliases; unmatched IDs are listed. The 2D exclusion and the 3D per-fold removal use this mapping |
| Slurm array or job-count limits | One array per arm, throttled with `%50` |

---

## 15. Verifications at implementation time

These are checks of facts, not open design choices. The response to each outcome is fixed above.

1. Coverage of cBioPortal's growth pattern and DNA purity fields for the included samples (Sections 4.1, 4.6, 14). Whether Supplementary Table 1 is available was checked on 2026-10-15; it is not needed.
2. The included organoid pair count after filtering, which goes into `filter_counts.csv` and corrects the grant text.
3. Whether R plus TensorFlow runs locally for the fixture (Section 13).
4. Whether micromamba runs on Alpine (Section 14).
5. Whether the pinned DepMap release includes the NextGen 3D screens, HCMI model identifiers and the non-essential control list (Sections 4.1, 4.5, 14).
6. Whether the Sanger HCMI organoid screens are downloadable as gene-level scores and can be put on the DepMap scale (Sections 4.5, 14).
