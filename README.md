# Survival Analysis of Breast Cancer (METABRIC)

This repository contains an **end-to-end survival analysis project** on breast cancer patients using the **METABRIC dataset**. The project focuses on understanding factors associated with **Overall Survival (OS)** and comparing traditional statistical models with modern machine-learning–based survival models. The analysis is implemented in **R**, and a detailed written report is provided explaining the data, methods, results, and interpretations.

## Project Overview

Breast cancer survival varies widely across patients depending on clinical and biological characteristics. The goals of this project are to:

- Explore relationships between patient characteristics and **overall survival**
- Compare survival across different patient groups using **Kaplan–Meier analysis**
- Build and interpret **Cox proportional hazards models**
- Apply **Random Survival Forests (RSF)** for improved prediction
- Evaluate and compare model performance using **Brier scores**
- Demonstrate **survival prediction** and **conditional survival estimation**

## Dataset

- **Source:** METABRIC (Molecular Taxonomy of Breast Cancer International Consortium) https://www.cbioportal.org/study/clinicalData?id=brca_metabric
- **Sample size:** ~2,500 breast cancer patients
- **Data types:**
  - Clinical variables (age, tumor stage, surgery type, treatments, receptor status)
  - Molecular subtypes (ER, PR, HER2, PAM50, Integrative Clusters)
  - Survival outcomes

### Survival Outcome Used
- **Overall Survival (OS):**
  - Time from diagnosis to death or last follow-up
  - Right-censored observations handled appropriately

## Exploratory Data Analysis

Key exploratory steps include:

- Descriptive statistics of age, tumor size, and receptor status
- Spearman correlation analysis to detect multicollinearity
- Kaplan–Meier survival curves for:
  - Tumor stage
  - Surgery type
  - ER and HER2 status
- Log-rank tests to assess survival differences between groups

## Methodology

### 1. Data Preprocessing
- Removed missing survival time/status
- Imputed missing values:
  - Median (numeric)
  - Mode (categorical)
  - **MICE** for variables with >5% missingness

### 2. Cox Proportional Hazards Models
- **Linear Cox model** with stepwise AIC variable selection
- Model diagnostics:
  - Martingale residuals (linearity)
  - Schoenfeld residuals (proportional hazards)
- **Non-linear Cox model**:
  - Natural cubic splines for age
  - Log-transformed tumor size

### 3. Random Survival Forests (RSF)
- Default RSF model
- Tuned RSF model:
  - Improved node size, mtry, and splitting
- Variable importance analysis

## Model Evaluation

Models are compared using the **Integrated Brier Score (IBS)**:

| Model | IBS (lower is better) |
|------|-----------------------|
| Kaplan–Meier (Null) | 0.174 |
| Linear Cox | 0.152 |
| Non-linear Cox | 0.149 |
| RSF (Default) | 0.148 |
| **RSF (Tuned)** | **0.147** |

**Random Survival Forest** achieved the best predictive performance.

## Survival Prediction

The project demonstrates:

- **12-year survival prediction** for representative patient profiles
- **Conditional survival** (e.g., survival beyond 12 years given survival up to 7 years)
- Interpretation of survival differences by:
  - Age
  - Molecular subtype
  - Tumor stage
  - Chemotherapy status

## Extensions: Genomic Data

The report also discusses extending the analysis by incorporating **high-dimensional gene expression data**, using:

- Variance-based gene filtering
- Regularized Cox models (Lasso / Elastic Net)
