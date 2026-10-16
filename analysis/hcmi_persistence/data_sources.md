# Data sources: HCMI persistence map

Every input was fetched by `analysis/hcmi_persistence/01_download.R` (Alpine job 34345892, commit
f0afbaf, 2026-10-16) except the two DepMap Public 26Q3 files, which the DepMap portal serves only
to a browser. Step 01 verifies those and records them with the rest. The full record, with paths,
is `raw/manifest.csv` under `$GRACKLE_DATA` on Alpine; the sha256 values below are copied from it.

## Expression: GDC

- Project HCMI-CMDC, data type "Gene Expression Quantification", workflow "STAR - Counts".
- GDC Data Release 46.0 (August 10, 2026). 1,153 STAR count files, each checked against the md5
  that the GDC files endpoint reports.
- File and sample metadata: `raw/gdc/gdc_files.tsv`
  (sha256 `0b8b23c3af8ff326ee7a808791d453c5aedbf67377e936bfac8976fe8c01b47b`).

## Clinical data, purity and mutations: cBioPortal

- Study `pancan_hcmi_2025` through the public REST API, accessed 2026-10-16.
- Sample and patient clinical attributes (growth pattern, culture media, treatment status, the
  TCGA/HCMI diagnosis code, age, sex, race, and the ABSOLUTE, ASCAT, HATCHet and PURPLE purity
  calls). `CONSENSUS_PURITY` is not used: it is not a purity estimate (r = -0.31 with ABSOLUTE).
- Mutations from `POST /molecular-profiles/pancan_hcmi_2025_mutations/mutations/fetch` with
  sample list `pancan_hcmi_2025_sequenced`. They annotate programs after fitting; they are not in
  the sample graph.

| File | sha256 |
|---|---|
| `raw/cbioportal/clinical_sample_long.tsv` | `ef66bd357913b0a31fb70fb1a66ad0df1417a388a3b8011d6d4b84061a4cef70` |
| `raw/cbioportal/clinical_patient_long.tsv` | `b6b6bfbe4c031bb58f10e3e5912ac376c651efbd6c78c4a77c558bef7659a71c` |
| `raw/cbioportal/mutations.tsv` | `b73cae3d211d1acce957a3f1cc1212edd4a54221a5305b9c2911af50ef447af4` |
| `raw/cbioportal/sequenced_samples.txt` | `f5c0ef107b0f0f5068df0a9351b9066fe500a8578d576cc57adbb0b890f5fbab` |

## Gene graphs: STRING

STRING v12.0, human (9606), from stringdb-downloads.org. The full links file supplies the
per-channel scores for `string_mech` (experiments and database channels); `combined_score` is used
for `string_all`.

| File | sha256 |
|---|---|
| `9606.protein.links.full.v12.0.txt.gz` | `07f9fa42ae5006ccb0b4694c17b97dd4d346c141bf769006cb6ffc4c8e04d016` |
| `9606.protein.links.v12.0.txt.gz` | `3e22f32572211aa341d5b4bd08d30c32e693e294603202120936872f87719d4f` |
| `9606.protein.aliases.v12.0.txt.gz` | `b65f730b993ed0c1bd72edf4565d3d425db42861101b29699704810e8f125680` |

## Gene graphs: DepMap co-essentiality

| File | Source | sha256 |
|---|---|---|
| `manual/depmap_nextgen/ScreenGeneEffect.csv` | DepMap Public 26Q3, portal download | `08e15ab81d7450db70c48c7781f312f6cebc068da7077bea775569fb516a3f04` |
| `manual/depmap_nextgen/Model.csv` | DepMap Public 26Q3, portal download | `65eb1895800e71a963a3ab56e5d08c8d006cb0d9e7aa68f7e34a3810a2a3cbc7` |
| `depmap/CRISPRGeneEffect.csv` | DepMap 24Q4 Public, figshare+ doi 10.25452/figshare.plus.27993248 | `3d8f3ec6dbf2db7ff834b79b508622ec0b226f3518003fe96ecf5a4fcf167e3b` |
| `depmap/Model.csv` | DepMap 24Q4 Public, figshare+ | `b7a0c1385e6cef30132b56aff61f1261d11e3f490490b355c430d32ee0dbdcfa` |
| `depmap/AchillesNonessentialControls.csv` | DepMap 24Q4 Public, figshare+ | `c7add12af5a37912997fdab0a7aabeef8d811b1ac7ee23dfcffa3a06471ce097` |
| `depmap/nextgen_screen_metadata.csv` | NextGen manuscript, figshare doi 10.6084/m9.figshare.29472362 | `af29e349dd19676029950d6799fec4efcefcb9b314f7622b4fcc6d564803472a` |
| `depmap/nextgen_model_metadata.csv` | NextGen manuscript, figshare | `a000a5ee1ebe2b2488a146b2a413ec4b93472488e6ad31b56126f522d652557c` |

Citation for the 26Q3 files: DepMap, Broad (2026). DepMap Public 26Q3. Dataset. Their md5s are
pinned in `depmap_nextgen_files.csv`.

**Which screens go into which arm.** Arm membership follows the culture format of the screen, read
from the NextGen `screen_metadata.csv` (`ScreenType`, `PassesQC`), not the model's growth pattern:

- `2DS` (traditional 2D) goes to `coess_2d`.
- `3DO` (organoid domes) and `3DN` (CNS spheroids) go to `coess_3d`.
- Everything else is excluded: 2D screens of NextGen models (`2DN`, `2DO`), screens that fail QC,
  and screens absent from `screen_metadata.csv`.

The 26Q3 file holds 1,611 screens: 1,386 are 2D (1,157 models), 108 are 3D (108 models), and 117
are excluded (76 of them are not in the manuscript's screen list). The file has at least 1,000
traditional models, so both arms come from this one Chronos run, and 24Q4 `CRISPRGeneEffect.csv`
is not used. Screens are averaged within a model.

**HCMI-derived models.** The three HCMI-derived 2D models are removed from `coess_2d`. In
`coess_3d`, each fold's own HCMI cases are removed before the graph is built. Models are tied to an
HCMI case by the HCMI identifier found in any text column of either model table; with the 26Q3
`HCMIID` column, 53 of the 108 3D models carry one. A model that is in neither model table cannot be checked for
HCMI origin, so it is dropped. One 3D model (`ACH-003330`, esophageal) is dropped this way, which
leaves 107. `gene_graph_summary.csv` records the count.

The non-essential control genes come from 24Q4. They only set the threshold for informative genes.

**Sanger organoid screens: not used.** The Sanger organoid CRISPR screens ("CRISPR Dependencies
2026", Cell Model Passports) were left out at the user's choice on 2026-10-16. They are available
as gene-level scores (CRISPRcleanR fold changes and BAGEL2 Bayes factors, not Chronos), but only
through a browser. `coess_3d` uses the DepMap 3D screens alone.

## Annotation

| File | Source | sha256 |
|---|---|---|
| `annotation/hgnc_complete_set_2024-10-01.txt` | HGNC monthly archive, 2024-10-01 | `517a9f4e3dff5a2d5b4fe806f955a2714797a2877b43efd0a6baa0ac4f1ecb6e` |
| `annotation/oncokb_cancerGeneList.txt` | OncoKB cancer gene list, retrieved 2026-10-16 | `576d3955b55989c9ab62c21c09a2405c4d997891c8acd7e62e6ad6272bf7b5c3` |
| `annotation/hallmark.rds` | msigdbr 26.1.1, collection H | `c3410780123b35c18231b49189e22be7c2b22ef379c7899e284680cffbc88f2a` |

## Spec Section 15 checks settled here

- **4. micromamba on Alpine:** yes. `setup_env.sbatch` (job 34344116) built the pinned environment.
  `environment.lock.linux-64.txt` lists it: 366 packages, R 4.5.3, TensorFlow 2.21.0.
- **5. The DepMap release's 3D screens, HCMI identifiers and control list:** yes, with two sources.
  DepMap Public 26Q3 holds all 108 QC-passing NextGen 3D screens, and its `Model.csv` holds the
  `HCMIID` column. The non-essential controls come from 24Q4. DepMap 24Q4 alone has no 3D screens:
  its screen map is 2D only.
- **6. Sanger HCMI organoid screens:** gene-level scores exist but are browser-only and not on the
  Chronos scale. Not used, by the user's choice. `combine3dSources()` would combine them with the
  DepMap screens only if their per-gene means correlate at 0.8 or above.
