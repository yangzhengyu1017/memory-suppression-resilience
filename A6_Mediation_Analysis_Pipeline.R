# ==============================================================================
# Script Name: A6_Mediation_Analysis_Pipeline.R
# Project:     Resilience fMRI PKU
# 
# Description: 
#   Comprehensive pipeline for mediation analyses used in the study.
#   Part 1: Evaluates competing psychometric mediation models (CTQ, TACQ, Resi).
#   Part 2: Batch processing for high-dimensional neurobehavioral ROI data, 
#           incorporating covariates (Age, Gender) standard to neuroimaging GLMs.
#
# Methodology:
#   - Non-parametric Bootstrapping (1000 simulations) via the 'mediation' package.
#   - Tidyverse workflow for robust batch processing and result extraction.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Environment Setup & Library Initialization
# ------------------------------------------------------------------------------
# Clear workspace and set seed for reproducible bootstrapping
rm(list=ls())
set.seed(12345) 

# Initialize required packages (using pacman for efficient package management)
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, mediation, performance, broom)

# ==============================================================================
# 1. Data Loading & Preprocessing
# ==============================================================================

tryCatch({
  clinical_data <- read_csv("data_all.csv", show_col_types = FALSE) %>% drop_na()
  neuro_data <- read_csv("ROI_behavior_NT_T.csv", show_col_types = FALSE)
  message("Data loaded successfully.")
}, error = function(e) {
  stop("Failed to load data. Please ensure the CSV files exist in the directory.\n", e)
})

# ==============================================================================
# PART 1: Clinical Competing Models (Without Covariates)
# Model A: CTQ (IV) -> TACQ (Mediator) -> Resi (DV)
# Model B: CTQ (IV) -> Resi (Mediator) -> TACQ (DV)
# ==============================================================================
message("\n>>> Executing Part 1: Clinical Competing Models...")

competing_mediators <- c("TACQ", "Resi")
clinical_results_list <- list()

for (med_var in competing_mediators) {
  # Dynamically assign the dependent variable based on the current mediator
  dv_var <- ifelse(med_var == "TACQ", "Resi", "TACQ") 
  
  # Construct base linear models (Path a and Path b/c')
  formula_M <- as.formula(paste(med_var, "~ CTQ"))
  formula_Y <- as.formula(paste(dv_var, "~", med_var, "+ CTQ"))
  
  model.M <- lm(formula_M, data = clinical_data)
  model.Y <- lm(formula_Y, data = clinical_data)
  
  # Execute non-parametric bootstrap mediation analysis
  # Parameters: treat (Independent Var), mediator (Mediating Var), boot (Enable bootstrapping)
  med_out <- mediate(model.M, model.Y, treat = 'CTQ', mediator = med_var, 
                     boot = TRUE, sims = 1000)
  
  # Extract core coefficients and store in the list
  res_sum <- summary(med_out)
  clinical_results_list[[med_var]] <- tibble(
    Model_Path    = paste("CTQ ->", med_var, "->", dv_var),
    Mediator      = med_var,
    ACME_Estimate = res_sum$d0,        # Average Causal Mediation Effect (Indirect)
    ACME_CI_lower = res_sum$d0.ci[1],
    ACME_CI_upper = res_sum$d0.ci[2],
    ACME_P_value  = res_sum$d0.p,
    ADE_Estimate  = res_sum$z0,        # Average Direct Effect
    ADE_P_value   = res_sum$z0.p,
    Prop_Mediated = res_sum$n0         # Proportion Mediated
  )
}

# Bind and display results
clinical_results <- bind_rows(clinical_results_list)
print(clinical_results %>% mutate_if(is.numeric, round, 4))


# ==============================================================================
# PART 2: Neurobehavioral Batch Mediation (With Covariates)
# Model: inhibitionRate (IV) -> Brain ROI (Mediator) -> CD_RISC_totalScore (DV)
# Covariates: Age, Gender
# ==============================================================================
message("\n>>> Executing Part 2: Neurobehavioral ROI Mediation (Batch Process)...")

# Data cleaning: Exclude specific cohorts (e.g., Group 3) and remove missing cases
neuro_clean <- neuro_data %>% 
  filter(Group != 3) %>% 
  drop_na()

# Dynamically extract all Brain ROI column names (assuming they start with "ROI_")
roi_biomarkers <- neuro_clean %>% select(starts_with("ROI_")) %>% colnames()

neuro_results_list <- list()

for (biomarker in roi_biomarkers) {
  
  # Path a: IV -> Mediator (controlling for Age, Gender)
  form_M <- as.formula(paste(biomarker, "~ inhibitionRate + Gender + Age"))
  model.M <- lm(form_M, data = neuro_clean)
  
  # Path b & c': Mediator & IV -> DV (controlling for Age, Gender)
  form_Y <- as.formula(paste("CD_RISC_totalScore ~", biomarker, "+ inhibitionRate + Gender + Age"))
  model.Y <- lm(form_Y, data = neuro_clean)
  
  # Bootstrap Mediation Analysis
  med_out <- mediate(model.M, model.Y, treat = 'inhibitionRate', mediator = biomarker, 
                     boot = TRUE, sims = 1000)
  
  res_sum <- summary(med_out)
  
  neuro_results_list[[biomarker]] <- tibble(
    ROI_Mediator  = biomarker,
    ACME_Estimate = res_sum$d0,
    ACME_CI_lower = res_sum$d0.ci[1],
    ACME_CI_upper = res_sum$d0.ci[2],
    ACME_P_value  = res_sum$d0.p,
    ADE_Estimate  = res_sum$z0,
    ADE_P_value   = res_sum$z0.p
  )
}

# Consolidate results
neuro_results <- bind_rows(neuro_results_list)

# Filter for ROIs with significant indirect effects (using raw p-value for demo purposes)
sig_neuro_results <- neuro_results %>% filter(ACME_P_value < 0.05)

message("\n--- Significant ROI Mediators (Uncorrected p < 0.05) ---")
print(sig_neuro_results %>% mutate_if(is.numeric, round, 4))

message("\n>>> All Mediation Operations Completed Successfully.")
