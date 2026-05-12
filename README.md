# memory-suppression-resilience
Analysis protocols investigating how memory suppression (Think/No-Think paradigm) fosters psychological resilience in young adults with childhood trauma.

# The Protective Role of Memory Suppression in Psychological Resilience among Young Adults with Childhood Trauma

## Overview
This repository contains the analysis scripts and protocol guidelines used in our study investigating the relationship between childhood trauma, psychological resilience, and memory suppression. 

The study integrates a large-scale behavioral survey with a functional MRI experiment involving a Think/No-Think (TNT) task. The fMRI cohort consists of three distinct groups: a high-resilience trauma group (HR), a low-resilience trauma group (LR), and a trauma-naïve healthy control group (HC).

## Repository Contents

### 1. fMRI First-Level Analysis (MATLAB / SPM12)
* **`A1_FirstLevel_Analysis_TNT.m`**: Performs first-level General Linear Model (GLM) specification, estimation, and contrast generation for the main Think/No-Think task.
* **`A4_FirstLevel_Perception_Ident.m`**: Executes first-level analysis for the post-TNT Perception Identification task, modeling specific condition contrasts.

### 2. fMRI Second-Level Analysis (MATLAB / SPM12)
* **`A2_SecondLevel_NT_T.m`**: Conducts a second-level single-sample T-test aggregating the No-Think vs. Think (NT-T) contrast across all three groups while covarying for age and gender.
* **`A3_group_ana_NT_T_ANOVA.m`**: Runs a second-level One-Way ANOVA to compare the NT-T contrast between the HR, LR, and HC groups.
* **`A5_group_ana_Perception_Ident_ANOVA_F_Test.m`**: Automates second-level ANOVA F-tests and post-hoc T-tests across 16 first-level contrasts for the Perception Identification task.

### 3. Functional Connectivity Protocol
* **`CONN_Analysis_Protocol.pdf`**: A step-by-step description of the graphical user interface operations performed in the CONN toolbox (version 22.v2407). It details the setup, denoising, and application of the generalized psychophysiological interaction (gPPI) model to assess seed-to-voxel connectivity.

### 4. Statistical & Mediation Analysis (R)
* **`A6_Mediation_Analysis_Pipeline.R`**: A comprehensive R pipeline for running non-parametric bootstrap mediation analyses. It includes:
    * Evaluation of competing clinical psychometric models (e.g., Childhood Trauma Questionnaire, Thought Control Ability, and Resilience).
    * Batch processing for high-dimensional neurobehavioral ROI data (incorporating Age and Gender as covariates).

## Dependencies & Requirements
* **MATLAB** with **SPM12** (Statistical Parametric Mapping).
* **DPABI** Toolbox (for fMRI preprocessing).
* **CONN Toolbox** (version 22.v2407).
* **R** (version 4.0 or higher) with the following packages: `tidyverse`, `mediation`, `lavaan`, `performance`, and `broom`.
