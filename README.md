# Peripheral Metabolic Profiling of Maternal Environmental Exposure

[![R Version](https://img.shields.io/badge/R-4.4.0-blue.svg)](https://cran.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete, fully reproducible R scripts and raw data files for the manuscript: **"Peripheral Metabolic Profiling of Maternal Environmental Exposure: Translational Insights for Early Screening and Precision Nutrition"** (Submitted to *Heliyon*).

Our study performs a secondary, machine-learning-driven data mining on an open-source murine metabolomic dataset (Project PR002895 / Study ST004597). By employing a rigorous triple-algorithm framework, we identified early-warning systemic biomarkers (KYNA) and quantified a Blood-Brain Barrier (BBB) amino acid transport competition index for precision nutritional nursing.

## 📂 Repository Contents

All files are placed in the root directory for easy access and execution.

### 1. Raw Data Files
*   `MSdata_ST004597_1.txt`: Raw positive ion mode (`POS`) mass spectrometry expression matrix.
*   `MSdata_ST004597_2.txt`: Raw negative ion mode (`NEG`) mass spectrometry expression matrix.
> **Note:** Clinical metadata is not uploaded manually. Our preprocessing script is designed to dynamically fetch the latest clinical cohort information directly from the NIH NMDR cloud API via `jsonlite`.

### 2. Reproducible R Scripts
*   `01_Data_Preprocessing_and_QC.R`: Handles data fusion, dynamic metadata matching, strict quality control (filtering missing values > 30%), TIC normalization, KNN imputation, and Log2 transformation. Generates the fundamental `Table_S1_Cleaned_Metabolomics_Matrix.csv`.
*   `02_Global_Profiling_and_Volcano.R`: Performs differential abundance analysis and generates the publication-ready Volcano plot and `Table_S2` & `Table_S3`.
*   `03_MachineLearning_and_ClinicalIndex.R`: Executes the triple-algorithm cross-validation framework (Random Forest, LASSO, and PLS-DA with LOOCV) to identify robust consensus biomarkers. Calculates and visualizes the BBB Transport Competition Index (BCAA/Tryptophan).
*   `04_GutBrain_Cometabolic_Network.R`: Constructs the targeted exploratory gut-brain axis co-metabolic network based on strict mathematical thresholds (Spearman |r| >= 0.75, p < 0.05).

## 🚀 How to Reproduce the Analysis

To fully reproduce the findings, Supplementary Tables, and high-resolution PDF figures:

1. **Clone this repository** to your local machine:
   ```bash
   git clone https://github.com/bingmoon/Maternal_TCS_Metabolomics.git
   ```
2. **Set your Working Directory** in R/RStudio to the cloned folder. 
3. **Run the scripts sequentially**:
   * **You MUST run `01_Data_Preprocessing_and_QC.R` first**. It generates the cleaned matrix required for all subsequent steps.
   * Then, run `02`, `03`, and `04` to generate the corresponding statistical results and figures. All outputs (PDFs and CSVs) will be saved directly in the same folder.

## 🛠️ Dependencies

All computational analyses were performed using **R software (version 4.4.0)**. Required packages:
*   **Data Manipulation:** `dplyr`, `stringr`, `jsonlite`, `reshape2`, `impute` (Bioconductor)
*   **Machine Learning & Statistics:** `randomForest`, `glmnet`, `caret`, `pls`, `Hmisc`
*   **Visualization & Network:** `ggplot2`, `ggrepel`, `ggpubr`, `ggVennDiagram`, `igraph`, `ggraph`, `tidygraph`

## ✉️ Contact & Authorship

**Xiaohuan Pei** (Corresponding Author)  
Department of Cardiovascular Medicine, The Fourth People's Hospital of Chengdu, Chengdu, Sichuan, China.  
Email: 421107707@qq.com  

## 📜 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
