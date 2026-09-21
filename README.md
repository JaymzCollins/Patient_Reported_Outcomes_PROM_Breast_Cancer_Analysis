# PROM Subgroup Analysis — Breast Cancer Pre-Treatment Quality of Life

Identifying clinically meaningful patient subgroups in pre-treatment breast cancer 
quality-of-life scores, using exploratory factor analysis and model-based recursive 
partitioning (conditional inference trees).

**Author:** Jaymz Collins 
**Supervisor:** Dr Peter Giles  
**R Version:** 4.5.1

## Background

Patient-Reported Outcome Measures (PROMs) — the EORTC QLQ-C30 and QLQ-BR23 
questionnaires — capture quality of life across physical, functional, symptom, and 
psychosocial domains in breast cancer patients prior to treatment. This project asks 
whether patients cluster into interpretable subgroups with meaningfully different 
baseline quality-of-life profiles, and which patient characteristics predict 
subgroup membership.

## Data

Baseline PROM scores and patient covariates (age, diagnosis, marital status, 
menopausal status, comorbidities, etc.) for breast cancer patients. 
**The dataset is not included in this repository due to patient confidentiality** — 
`Data_PROM_Baseline_updateCF.xlsx` is expected in the working directory when running 
the script, but is not published here or elsewhere.

## Method

1. **Data cleaning** — sequential, logged cleaning pipeline: remove patients with no 
   PROM data, drop high-missingness items (hair loss, sexual enjoyment), remove rows 
   with missing scores/covariates, reverse-score symptom scales so higher is always 
   better, convert categorical covariates to factors. Final sample: **1,391 patients** 
   (from 1,727).
2. **Exploratory Factor Analysis (EFA)** — KMO and Bartlett's test confirm 
   factorability (KMO = 0.94); parallel analysis suggests 7 factors, of which 3 were 
   retained on clinical and statistical grounds. Maximum likelihood extraction with 
   oblimin (oblique) rotation, since physical and psychological domains are expected 
   to correlate.
3. **Composite scores** — two composites built from EFA-informed item groupings:
   - **ML1 (Physical functioning)**: physical functioning, dyspnoea (reversed), sexual functioning (α = 0.55)
   - **ML6+ML3 (Psychosocial wellbeing)**: social functioning, body image, future perspective, emotional functioning (α = 0.78)
4. **Model-based recursive partitioning** — conditional inference trees (`ctree`, 
   `partykit`) predict each composite from patient covariates, on a 70/30 train/test 
   split (minbucket = 50, maxdepth = 4).
5. **Validation** — 10-fold cross-validation for stability, comparison against a null 
   (mean-only) model, and subgroup differences checked against the EORTC Minimum 
   Clinically Important Difference (MCID = 10 points).

## Results

| Model | Test RMSE | Test R² | vs. null model | CV mean RMSE (SD) |
|---|---|---|---|---|
| ML1 — Physical functioning | 16.343 | 0.152 | 7.9% improvement | 16.337 (1.128) |
| ML6+ML3 — Psychosocial wellbeing | 19.229 | 0.108 | 6.0% improvement | 20.076 (1.161) |

Both models outperformed the null model and identified subgroups with mean 
differences well exceeding the MCID threshold:
- **ML1** node means ranged 50.17–77.99 (2.8× MCID)
- **ML6** node means ranged 52.73–75.78 (2.3× MCID)

This indicates the patient subgroups identified by covariates (diagnosis, age, 
marital status, comorbidities, etc.) differ by clinically — not just statistically — 
meaningful amounts in baseline quality of life.

## Running it

```r
# Required packages
install.packages(c("readxl", "dplyr", "partykit", "psych", "GPArotation", "ggparty", "ggplot2"))
```

Place `Data_PROM_Baseline_updateCF.xlsx` in your working directory, then run 
`PROMs_Breast_Cancer.R` top to bottom. The script validates package 
installation, data presence, and expected columns at each stage, and stops with an 
informative message if anything is missing.

## Tools

R, dplyr, psych (EFA, KMO, Cronbach's alpha), GPArotation, partykit (ctree), ggparty, ggplot2
