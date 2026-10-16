## Six cases. c1: one tumor + organoid. c2: primary + metastasis, organoid from the metastasis.
## c3: two tumors, organoid with unknown parent. c4, c5: unpaired tumors. c6: tumor + organoid.
## Annotations include NA age, NA sex, a constant ancestry and one tumor without purity (t3b).
toy_samples <- function() {
  data.frame(
    profile_id = c("t1", "o1", "t2a", "t2b", "o2", "t3a", "t3b", "o3", "t4", "t5", "t6", "o6"),
    case_id = c("c1", "c1", "c2", "c2", "c2", "c3", "c3", "c3", "c4", "c5", "c6", "c6"),
    sample_id = c("c1-01", "c1-85", "c2-01", "c2-06", "c2-85", "c3-01", "c3-06", "c3-85",
                  "c4-01", "c5-01", "c6-01", "c6-85"),
    parent_sample_id = c(NA, "c1-01", NA, NA, "c2-06", NA, NA, NA, NA, NA, NA, "c6-01"),
    role = c("tumor", "organoid", "tumor", "tumor", "organoid", "tumor", "tumor", "organoid",
             "tumor", "tumor", "tumor", "organoid"),
    cancer_type = c("COAD", "COAD", "PAAD", "PAAD", "PAAD", "COAD", "COAD", "COAD", "BRCA", "COAD",
                    "PAAD", "PAAD"),
    tumor_descriptor = c("primary", "primary", "primary", "metastatic", "metastatic", "primary",
                         "metastatic", NA, "primary", "primary", "primary", "primary"),
    prior_treatment = c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, NA, TRUE, FALSE, FALSE),
    age = c(60, 60, 55, 55, 55, NA, NA, NA, 70, 62, 48, 48),
    sex = c("F", "F", "M", "M", "M", "F", "F", "F", NA, "M", "F", "F"),
    ancestry = "EUR",
    purity_dna = c(0.5, NA, 0.7, 0.4, NA, 0.6, NA, NA, 0.8, 0.45, 0.65, NA),
    stringsAsFactors = FALSE)
}
