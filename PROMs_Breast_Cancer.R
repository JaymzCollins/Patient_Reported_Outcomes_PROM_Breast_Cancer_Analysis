# =============================================================================
# PROM Subgroup Analysis — Breast Cancer Pre-Treatment Quality of Life Scores
# Author: Jaymz Collins | Cardiff University | MET584 Data Project
# Supervisor: Dr Peter Giles
# R Version: 4.5.1
# =============================================================================



# ── Load Packages ----
# Checks all required packages are installed before loading
# Stops with an informative message if any are missing rather than
# failing silently mid-analysis
required_packages <- c(
  "readxl",      # Reading Excel input file
  "dplyr",       # Sequential data cleaning pipeline
  "partykit",    # ctree model-based recursive partitioning
  "psych",       # EFA, KMO, Bartlett's test, Cronbach's alpha
  "GPArotation", # Oblimin rotation for EFA
  "ggparty",     # Publication-quality tree plots
  "ggplot2"      # Plotting support
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(paste(
    "The following required packages are not installed:",
    paste(missing_packages, collapse = ", "),
    "\nInstall with: install.packages(c('",
    paste(missing_packages, collapse = "','"), "'))"
  ))
}

invisible(lapply(required_packages, library, character.only = TRUE))
message("All packages loaded successfully.")


# ── Load Data ----
# Reads the raw Excel file into R
# Raw_Data is preserved untouched throughout — all cleaning is applied
# to a separate working copy (Clean_Data) to allow the original to be
# inspected at any point without re-reading the file

if (!file.exists("Data_PROM_Baseline_updateCF.xlsx")) {
  stop(paste(
    "Data file not found.",
    "\nCheck that 'Data_PROM_Baseline_updateCF.xlsx' is in your working directory:",
    getwd()
  ))
}

# tryCatch wraps read_excel to return an informative error if the file
# is corrupt, password-protected, or in an unreadable format
Raw_Data <- tryCatch(
  read_excel("Data_PROM_Baseline_updateCF.xlsx"),
  error = function(e) stop(paste("Failed to read Excel file:", e$message))
)

if (nrow(Raw_Data) == 0) stop("Data file loaded but contains no rows.")
message(paste("Data loaded successfully:", nrow(Raw_Data), "rows,",
              ncol(Raw_Data), "columns."))

# Create working copy — Raw_Data is preserved untouched throughout
Clean_Data <- Raw_Data


# ── Data Cleaning ----
# Performs sequential cleaning of the raw dataset
# Each step is logged with before/after row counts for transparency
# Expected final sample size: 1,391 patients

message("\n--- Beginning data cleaning ---")

# Define all 23 PROM item column names for validation
cols <- c("ql","pf","rf","ef","cf","sf","fa","nv","pa","dy","sl",
          "ap","co","di","fi","brst","brbi","brbs","brfu",
          "brsee","brsef","bras","brhl")

# Verify all expected PROM columns are present before proceeding
missing_cols <- cols[!cols %in% names(Clean_Data)]
if (length(missing_cols) > 0) {
  stop(paste("Expected PROM columns missing from data:",
             paste(missing_cols, collapse = ", ")))
}


# Step 1: Remove rows with no PROM data at all
# Patients with zero PROM responses contribute no outcome data
# rowSums(!is.na()) counts non-missing PROM values per row
# Rows where this count is 0 are entirely missing and removed
# Expected: 249 rows removed, sample = 1,478
n_before   <- nrow(Clean_Data)
Clean_Data <- Clean_Data[rowSums(!is.na(Clean_Data[, cols])) > 0, ]
message(paste("Step 1 — Removed", n_before - nrow(Clean_Data),
              "rows with no PROM data. Remaining:", nrow(Clean_Data)))


# Step 2: Remove high-missingness PROM columns
# brhl (hair loss): 65% missing
# brsee (sexual enjoyment): 38% missing
# Missingness in these items is not at random; excluding columns rather
# than rows preserves the full remaining sample
Clean_Data <- Clean_Data[, !(names(Clean_Data) %in% c("brhl", "brsee"))]
message("Step 2 — Removed brhl (65% missing) and brsee (38% missing).")


# Step 3: Remove rows with missing values on remaining BR23 items
# ctree applies listwise deletion at the modelling stage — any patient
# missing a model variable is excluded entirely. Removing these rows
# here makes the exclusion transparent and controllable
# Expected: 4 rows removed, sample = 1,474
check_cols    <- c("brbi","brbs","brfu","brsef","bras")
missing_check <- check_cols[!check_cols %in% names(Clean_Data)]
if (length(missing_check) > 0) {
  stop(paste("Expected columns for missingness check not found:",
             paste(missing_check, collapse = ", ")))
}

n_before   <- nrow(Clean_Data)
Clean_Data <- Clean_Data[complete.cases(Clean_Data[, check_cols]), ]
message(paste("Step 3 — Removed", n_before - nrow(Clean_Data),
              "rows with missing PROM scores. Remaining:", nrow(Clean_Data)))


# Step 4: Reverse symptom scales
# EORTC symptom scales are scored so higher = worse outcome
# Reversing (100 - score) aligns all items so higher = better outcome
# This is required before constructing composite scores
c30_symptoms  <- c("fa","nv","pa","dy","sl","ap","co","di","fi")
br23_symptoms <- c("bras","brbs","brst")
all_symptoms  <- c(c30_symptoms, br23_symptoms)

missing_symp <- all_symptoms[!all_symptoms %in% names(Clean_Data)]
if (length(missing_symp) > 0) {
  stop(paste("Expected symptom columns not found:",
             paste(missing_symp, collapse = ", ")))
}

# Create reversed columns with _rev suffix
for (col in all_symptoms) {
  Clean_Data[[paste0(col, "_rev")]] <- 100 - Clean_Data[[col]]
}

# Remove original unreversed columns to avoid confusion
Clean_Data <- Clean_Data[, !(names(Clean_Data) %in% all_symptoms)]
message("Step 4 — Symptom scales reversed. Higher scores now consistently indicate better outcomes.")


# Step 5: Convert categorical covariates to factors
# CRITICAL: Factor conversion must occur BEFORE covariate cleaning (Step 6)
# Reason: values outside the defined levels become NA upon factor conversion
# These new NAs are then removed cleanly in Step 6
# If Step 6 ran first, out-of-range codes would not be caught

# marital_status: 0 = no relationship, 1 = married, 2 = divorced, 3 = widowed
Clean_Data$marital_status <- factor(Clean_Data$marital_status,
                                    levels = c(0,1,2,3),
                                    labels = c("No relationship","Married","Divorced","Widowed"))

# smokingstatus: 0 = never smoked, 1 = current smoker, 2 = ex-smoker
Clean_Data$smokingstatus <- factor(Clean_Data$smokingstatus,
                                   levels = c(0,1,2),
                                   labels = c("No","Yes","Ex-smoker"))

# alcohol: 0 = none, 1 = occasional, 2 = weekly
Clean_Data$alcohol <- factor(Clean_Data$alcohol,
                             levels = c(0,1,2),
                             labels = c("None","Occasionally","Weekly"))

# diagnosis: 1 = breast cancer, 2 = DCIS, 3 = fibroadenoma, 4 = other benign
Clean_Data$diagnosis <- factor(Clean_Data$diagnosis,
                               levels = c(1,2,3,4),
                               labels = c("Breast cancer","DCIS","Fibroadenoma","Other"))

# education: 1 = low, 2 = middle, 3 = high
Clean_Data$education <- factor(Clean_Data$education,
                               levels = c(1,2,3),
                               labels = c("Low","Middle","High"))

# Binary covariates: 0 = No, 1 = Yes
Clean_Data$menopause_yn         <- factor(Clean_Data$menopause_yn,         levels = c(0,1), labels = c("No","Yes"))
Clean_Data$pre_op               <- factor(Clean_Data$pre_op,               levels = c(0,1), labels = c("No","Yes"))
Clean_Data$breastcancer_first   <- factor(Clean_Data$breastcancer_first,   levels = c(0,1), labels = c("No","Yes"))
Clean_Data$cancer_kind_family_1 <- factor(Clean_Data$cancer_kind_family_1, levels = c(0,1), labels = c("No","Yes"))
Clean_Data$comorb_depression    <- factor(Clean_Data$comorb_depression,    levels = c(0,1), labels = c("No","Yes"))

message("Step 5 — Categorical covariates converted to factors.")


# Step 6: Remove rows with missing values on core covariates
# Only seven covariates are included in the cleaning step
# Remaining covariate missingness (e.g. pregnancy_number) is handled
# by ctree's listwise deletion at the modelling stage
# Expected: 83 rows removed, final sample = 1,391
covariates_to_clean <- c("age","marital_status","alcohol","bmi",
                         "smokingstatus","diagnosis","education")

missing_cov <- covariates_to_clean[!covariates_to_clean %in% names(Clean_Data)]
if (length(missing_cov) > 0) {
  stop(paste("Expected covariate columns not found:",
             paste(missing_cov, collapse = ", ")))
}

n_before   <- nrow(Clean_Data)
Clean_Data <- Clean_Data %>%
  dplyr::filter(complete.cases(across(all_of(covariates_to_clean))))
message(paste("Step 6 — Removed", n_before - nrow(Clean_Data),
              "rows with missing covariate data. Final sample:", nrow(Clean_Data)))

# Final sample size checks
if (nrow(Clean_Data) == 0)  stop("Cleaning removed all rows — check input data.")
if (nrow(Clean_Data) < 100) warning(paste("Final sample is very small:", nrow(Clean_Data), "rows."))

cat("Final sample size:", nrow(Clean_Data), "\n")  # Expected: 1391


# ── Exploratory Factor Analysis (EFA) ----

# Step 1: Select PROM items for EFA
prom_items <- c(
  "pf", "rf", "dy_rev", "pa_rev", "fa_rev",
  "bras_rev", "ql", "brbs_rev", "di_rev",
  "co_rev", "nv_rev", "brst_rev", "ap_rev",
  "ef", "brfu", "sf", "cf", "brbi",
  "fi_rev", "brsef", "sl_rev"
)

# Verify all PROM items are present after cleaning
missing_items <- prom_items[!prom_items %in% names(Clean_Data)]
if (length(missing_items) > 0) {
  stop(paste("PROM items not found in cleaned data:",
             paste(missing_items, collapse = ", ")))
}

df_prom <- Clean_Data[, prom_items]

# Check for residual missingness — should be zero after cleaning
n_missing <- sum(is.na(df_prom))
if (n_missing > 0) {
  warning(paste(n_missing, "missing values found in PROM items. Check cleaning pipeline."))
}


# Step 2: Check factorability
# KMO (Kaiser-Meyer-Olkin): values > 0.90 indicate exceptional suitability
# Bartlett's test: significant p confirms correlation structure is suitable
kmo_result <- KMO(df_prom)
message(paste("KMO sampling adequacy:", round(kmo_result$MSA, 3)))
# Expected: KMO = 0.94

if (kmo_result$MSA < 0.60) {
  warning(paste("KMO of", round(kmo_result$MSA, 3),
                "is below the acceptable threshold of 0.60.",
                "Factor analysis may not be appropriate."))
}

cortest.bartlett(df_prom)


# Step 3: Parallel analysis
# Compares real data eigenvalues against eigenvalues from randomly
# simulated datasets of the same dimensions
# fm = "ml": maximum likelihood — consistent with subsequent EFA
# fa = "fa": factor analysis rather than PCA
# n.iter = 100: 100 simulated datasets for stable comparison
# set.seed() must be called immediately before fa.parallel() to ensure
# the same 100 simulated datasets are generated each run
set.seed(42)
parallel <- fa.parallel(
  df_prom,
  fm     = "ml",   # Maximum likelihood estimation
  fa     = "fa",   # Factor analysis (not PCA)
  n.iter = 100,    # Number of simulated datasets for comparison
  SMC    = TRUE,
  main   = "Parallel Analysis — Number of Factors"
)

cat("\nParallel analysis suggests:", parallel$nfact, "factors\n")
# Expected: 7 factors suggested, 3 retained on clinical and statistical grounds


# Step 4: Run EFA
# oblimin rotation: allows factors to correlate — appropriate because
# physical and psychological health domains are related constructs
# fm = "ml": maximum likelihood estimation
efa_result <- tryCatch(
  fa(
    df_prom,
    nfactors = 7,          # Suggested by parallel analysis
    rotate   = "oblimin",  # Oblique rotation — factors allowed to correlate
    fm       = "ml",       # Maximum likelihood estimation
    SMC      = TRUE
  ),
  error = function(e) stop(paste("EFA failed:", e$message))
)


# Step 5: Inspect factor loadings
# cutoff = 0.30: items loading below this on all factors are not assigned
print(efa_result$loadings, cutoff = 0.30, sort = TRUE)

cat("\nVariance explained:\n")
print(efa_result$Vaccounted)

# Factor intercorrelations — non-zero values confirm oblimin was appropriate
cat("\nFactor correlations:\n")
print(efa_result$Phi)


# Step 6: Create composite scores

# Physical functioning composite (ML1)
# Items: pf (physical functioning), dy_rev (dyspnoea), brsef (sexual functioning)
ml1_items   <- c("pf", "dy_rev", "brsef")
missing_ml1 <- ml1_items[!ml1_items %in% names(Clean_Data)]
if (length(missing_ml1) > 0) {
  stop(paste("ML1 composite items not found:", paste(missing_ml1, collapse = ", ")))
}
Clean_Data$ml1_physical <- rowMeans(Clean_Data[, ml1_items], na.rm = TRUE)

# Psychosocial wellbeing composite (ML6+ML3)
# ML6 and ML3 merged on clinical and statistical grounds:
#   Clinical: emotional functioning (ef) sits naturally within the
#             psychosocial domain alongside social functioning,
#             body image and future perspective
# Items: sf (social functioning), brbi (body image),
#        brfu (future perspective), ef (emotional functioning)
ml6_items   <- c("sf", "brbi", "brfu", "ef")
missing_ml6 <- ml6_items[!ml6_items %in% names(Clean_Data)]
if (length(missing_ml6) > 0) {
  stop(paste("ML6 composite items not found:", paste(missing_ml6, collapse = ", ")))
}
Clean_Data$ml6_psychosoc <- rowMeans(Clean_Data[, ml6_items], na.rm = TRUE)


# Step 7: Check composite reliability
cat("\n--- ML1 Physical Functioning ---\n")
alpha(Clean_Data[, ml1_items])
# Expected: alpha = 0.55

cat("\n--- ML6+ML3 Psychosocial Wellbeing ---\n")
alpha(Clean_Data[, ml6_items])
# Expected: alpha = 0.78

# Confirm composites were created successfully
if (all(is.na(Clean_Data$ml1_physical)))  stop("ML1 composite contains only missing values.")
if (all(is.na(Clean_Data$ml6_psychosoc))) stop("ML6 composite contains only missing values.")


# ── Model-Based Recursive Partitioning ----

# Step 1: Train/test split
# 70/30 split: conventional proportion for moderately sized datasets
#   Training set (70%): used to fit the ctree models
#   Test set (30%):     held out for unbiased performance evaluation
# set.seed() ensures the same split is produced on every run
set.seed(42)
train_idx <- sample(seq_len(nrow(Clean_Data)), size = floor(0.7 * nrow(Clean_Data)))
df_train  <- Clean_Data[train_idx, ]   # Expected: 973 patients
df_test   <- Clean_Data[-train_idx, ]  # Expected: 418 patients

cat("Training set:", nrow(df_train), "\n")  # Expected: 973
cat("Test set:",     nrow(df_test),  "\n")  # Expected: 418

if (nrow(df_train) == 0) stop("Training set is empty after split.")
if (nrow(df_test)  == 0) stop("Test set is empty after split.")


# Step 2: Fit ML1 Physical Functioning tree
# Covariates selected on clinical grounds — factors plausibly affecting
# baseline physical capacity prior to treatment
# minbucket = 50: minimum patients in any terminal node — ensures
#   clinically interpretable subgroup sizes and stable estimates
# maxdepth = 4: maximum splits from root to terminal node — prevents
#   overfitting while allowing meaningful subgroup structure
tree_ml1 <- tryCatch(
  ctree(
    ml1_physical ~ age + bmi + pre_op + diagnosis,
    data    = df_train,
    control = ctree_control(minbucket = 50, maxdepth = 4)
  ),
  error = function(e) stop(paste("ML1 ctree failed:", e$message))
)


# Step 3: Evaluate ML1 on test set
# RMSE: average distance between predicted and actual scores (0-100 scale)
# R²: proportion of outcome variance explained by the model
pred_ml1 <- tryCatch(
  predict(tree_ml1, newdata = df_test),
  error = function(e) stop(paste("ML1 prediction failed:", e$message))
)

rmse_ml1 <- sqrt(mean((df_test$ml1_physical - pred_ml1)^2, na.rm = TRUE))
r2_ml1   <- 1 - sum((df_test$ml1_physical - pred_ml1)^2, na.rm = TRUE) /
  sum((df_test$ml1_physical - mean(df_test$ml1_physical, na.rm = TRUE))^2, na.rm = TRUE)

cat("\n--- ML1 Physical Functioning ---\n")
cat("RMSE:", round(rmse_ml1, 3), "\n")  # Expected: 16.343
cat("R²:  ", round(r2_ml1,   3), "\n")  # Expected:  0.152


# Step 4: Fit ML6+ML3 Psychosocial Wellbeing tree
# Covariates selected on clinical grounds — broader set reflecting the
# multiple dimensions known to influence psychosocial wellbeing at diagnosis
tree_ml6 <- tryCatch(
  ctree(
    ml6_psychosoc ~ diagnosis + marital_status + education + age +
      menopause_yn + breastcancer_first + cancer_kind_family_1 +
      comorb_depression + pregnancy_number + pre_op,
    data    = df_train,
    control = ctree_control(minbucket = 50, maxdepth = 4)
  ),
  error = function(e) stop(paste("ML6 ctree failed:", e$message))
)


# Step 5: Evaluate ML6 on test set
pred_ml6 <- tryCatch(
  predict(tree_ml6, newdata = df_test),
  error = function(e) stop(paste("ML6 prediction failed:", e$message))
)

rmse_ml6 <- sqrt(mean((df_test$ml6_psychosoc - pred_ml6)^2, na.rm = TRUE))
r2_ml6   <- 1 - sum((df_test$ml6_psychosoc - pred_ml6)^2, na.rm = TRUE) /
  sum((df_test$ml6_psychosoc - mean(df_test$ml6_psychosoc, na.rm = TRUE))^2, na.rm = TRUE)

cat("\n--- ML6+ML3 Psychosocial Wellbeing ---\n")
cat("RMSE:", round(rmse_ml6, 3), "\n")  # Expected: 19.229
cat("R²:  ", round(r2_ml6,   3), "\n")  # Expected:  0.108


# Step 6: Node summary
# Mean composite score for each terminal node in the training set
# These describe the quality of life profile of each identified subgroup
cat("\nML1 Node means (training set):\n")
print(aggregate(ml1_physical ~ predict(tree_ml1, newdata = df_train, type = "node"),
                data = df_train, FUN = mean))
# Expected nodes: 3, 5, 6, 9, 10, 11
# Expected means: 77.99, 74.11, 69.17, 66.81, 57.54, 50.17

cat("\nML6 Node means (training set):\n")
print(aggregate(ml6_psychosoc ~ predict(tree_ml6, newdata = df_train, type = "node"),
                data = df_train, FUN = mean))
# Expected nodes: 3, 5, 6, 8, 9
# Expected means: 52.73, 65.89, 59.22, 75.78, 67.83


# ── Model Validation ----

# Step 1: 10-fold cross-validation
# The dataset is divided into 10 folds; the model is trained on 9 folds
# and tested on the held-out fold, cycling through all folds
# This checks whether results are stable across different data partitions
# rather than an artefact of the particular 70/30 split chosen
# Implemented manually using partykit — caret's ctree method requires the
# legacy party package rather than partykit and is therefore incompatible
set.seed(42)
n     <- nrow(Clean_Data)
folds <- sample(rep(1:10, length.out = n))  # Assign each row to a fold

# Storage for per-fold results
cv_results <- data.frame(
  fold     = 1:10,
  rmse_ml1 = NA,
  rmse_ml6 = NA,
  rsq_ml1  = NA,
  rsq_ml6  = NA
)


# Step 2: Run 10-fold CV
for (i in 1:10) {
  
  cat("Fold", i, "of 10...\n")
  
  fold_train <- Clean_Data[folds != i, ]  # Train on all folds except i
  fold_test  <- Clean_Data[folds == i, ]  # Test on fold i
  
  # Guard against empty folds
  if (nrow(fold_train) == 0 | nrow(fold_test) == 0) {
    warning(paste("Fold", i, "produced an empty train or test set. Skipping."))
    next
  }
  
  # Fit ML1 tree — tryCatch allows remaining folds to continue if one fails
  t1 <- tryCatch(
    ctree(
      ml1_physical ~ age + bmi + pre_op + diagnosis,
      data    = fold_train,
      control = ctree_control(minbucket = 50, maxdepth = 4)
    ),
    error = function(e) {
      warning(paste("ML1 fold", i, "failed:", e$message))
      NULL  # Return NULL so downstream code can check and skip
    }
  )
  
  # Fit ML6 tree
  t6 <- tryCatch(
    ctree(
      ml6_psychosoc ~ diagnosis + marital_status + education + age +
        menopause_yn + breastcancer_first + cancer_kind_family_1 +
        comorb_depression + pregnancy_number + pre_op,
      data    = fold_train,
      control = ctree_control(minbucket = 50, maxdepth = 4)
    ),
    error = function(e) {
      warning(paste("ML6 fold", i, "failed:", e$message))
      NULL
    }
  )
  
  # Calculate RMSE and R² for folds where trees fitted successfully
  if (!is.null(t1)) {
    p1 <- predict(t1, newdata = fold_test)
    cv_results$rmse_ml1[i] <- sqrt(mean((fold_test$ml1_physical - p1)^2, na.rm = TRUE))
    cv_results$rsq_ml1[i]  <- 1 - sum((fold_test$ml1_physical - p1)^2, na.rm = TRUE) /
      sum((fold_test$ml1_physical - mean(fold_test$ml1_physical, na.rm = TRUE))^2, na.rm = TRUE)
  }
  
  if (!is.null(t6)) {
    p6 <- predict(t6, newdata = fold_test)
    cv_results$rmse_ml6[i] <- sqrt(mean((fold_test$ml6_psychosoc - p6)^2, na.rm = TRUE))
    cv_results$rsq_ml6[i]  <- 1 - sum((fold_test$ml6_psychosoc - p6)^2, na.rm = TRUE) /
      sum((fold_test$ml6_psychosoc - mean(fold_test$ml6_psychosoc, na.rm = TRUE))^2, na.rm = TRUE)
  }
}


# Step 3: CV summary
# SD of R² indicates fold-level stability — high SD suggests results
# may be sensitive to which patients were in the training set
cat("\n--- ML1 Physical Functioning ---\n")
cat("Mean RMSE:", round(mean(cv_results$rmse_ml1, na.rm = TRUE), 3), "\n")  # Expected: 16.337
cat("SD RMSE:  ", round(sd(cv_results$rmse_ml1,   na.rm = TRUE), 3), "\n")  # Expected:  1.128
cat("Mean R²:  ", round(mean(cv_results$rsq_ml1,  na.rm = TRUE), 3), "\n")  # Expected:  0.119
cat("SD R²:    ", round(sd(cv_results$rsq_ml1,    na.rm = TRUE), 3), "\n")  # Expected:  0.085

cat("\n--- ML6+ML3 Psychosocial Wellbeing ---\n")
cat("Mean RMSE:", round(mean(cv_results$rmse_ml6, na.rm = TRUE), 3), "\n")  # Expected: 20.076
cat("SD RMSE:  ", round(sd(cv_results$rmse_ml6,   na.rm = TRUE), 3), "\n")  # Expected:  1.161
cat("Mean R²:  ", round(mean(cv_results$rsq_ml6,  na.rm = TRUE), 3), "\n")  # Expected:  0.105
cat("SD R²:    ", round(sd(cv_results$rsq_ml6,    na.rm = TRUE), 3), "\n")  # Expected:  0.041


# Step 4: Null model comparison
# The null model predicts the training set mean for every patient regardless
# of their characteristics — representing zero analytical effort
# If the ctree models cannot outperform the null, they add no predictive value
null_mean_ml1 <- mean(df_train$ml1_physical,  na.rm = TRUE)
null_mean_ml6 <- mean(df_train$ml6_psychosoc, na.rm = TRUE)

null_rmse_ml1 <- sqrt(mean((df_test$ml1_physical  - null_mean_ml1)^2, na.rm = TRUE))
null_rmse_ml6 <- sqrt(mean((df_test$ml6_psychosoc - null_mean_ml6)^2, na.rm = TRUE))

cat("\n--- Null Model Comparison ---\n")
cat("ML1 Null RMSE:", round(null_rmse_ml1, 3), "\n")  # Expected: 17.753
cat("ML6 Null RMSE:", round(null_rmse_ml6, 3), "\n")  # Expected: 20.458
cat("ML1 improvement:", round((null_rmse_ml1 - rmse_ml1) / null_rmse_ml1 * 100, 1), "%\n")  # Expected: 7.9%
cat("ML6 improvement:", round((null_rmse_ml6 - rmse_ml6) / null_rmse_ml6 * 100, 1), "%\n")  # Expected: 6.0%

# Warn if either model fails to outperform the null
if (rmse_ml1 >= null_rmse_ml1) {
  warning("ML1 model does not outperform the null model. Review covariate selection.")
}
if (rmse_ml6 >= null_rmse_ml6) {
  warning("ML6 model does not outperform the null model. Review covariate selection.")
}


# Step 5: MCID analysis
# The EORTC Minimum Clinically Important Difference (MCID) is 10 points
# on the 0-100 scale — the smallest change a patient would themselves
# perceive as meaningful (Osoba et al., 1998)
# Node mean range is compared against this threshold to confirm subgroup
# differences are not merely statistical artefacts
ml1_nodes <- c(77.99, 74.11, 69.17, 66.81, 57.54, 50.17)
ml6_nodes <- c(52.73, 65.89, 59.22, 75.78, 67.83)

ml1_range <- max(ml1_nodes) - min(ml1_nodes)
ml6_range <- max(ml6_nodes) - min(ml6_nodes)

cat("\n--- MCID Analysis ---\n")
cat("ML1 node mean range:", round(ml1_range, 2), "\n")              # Expected: 27.82
cat("ML1 MCID multiplier:", round(ml1_range / 10, 1), "x MCID\n")  # Expected:  2.8x
cat("ML6 node mean range:", round(ml6_range, 2), "\n")              # Expected: 23.05
cat("ML6 MCID multiplier:", round(ml6_range / 10, 1), "x MCID\n")  # Expected:  2.3x

if (ml1_range < 10) warning("ML1 node mean range does not exceed the MCID threshold.")
if (ml6_range < 10) warning("ML6 node mean range does not exceed the MCID threshold.")