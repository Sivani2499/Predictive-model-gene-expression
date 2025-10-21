# Identification of Potential Therapeutic Targets in Idiopathic Pulmonary Fibrosis (IPF)

## Project Overview

Idiopathic Pulmonary Fibrosis (IPF) is a chronic, progressive, and fatal lung disease characterized by the scarring of lung tissue. Despite advancements, effective therapeutic targets remain elusive. This project aims to leverage advanced machine learning techniques on gene expression data to identify robust gene signatures that can accurately distinguish IPF patients from healthy individuals, thereby pinpointing potential therapeutic targets.

This rigorous analytical pipeline focuses on preventing data leakage and ensuring the reliability of findings through cross-validation, stability selection, and a dedicated holdout test set.

## Dataset

This analysis utilizes a gene expression dataset (RNA-seq / microarray counts) and associated clinical metadata for IPF patients and healthy controls.

-   **`Raw_count_IPF.csv`**: Contains raw gene expression counts.
-   **`metadata_modified.csv`**: Contains sample-specific metadata, including `Disease_status` (Normal/Disease).

## Analytical Workflow

The pipeline consists of several key steps designed for robustness and generalizability:

1.  **Data Loading and Preprocessing**: Initial setup, including assigning gene names as row identifiers and aligning metadata with expression data.
2.  **Normalization**: Variance Stabilizing Transformation (VST) using `DESeq2` to ensure robust gene expression values.
3.  **Data Splitting**: The dataset is initially partitioned into an independent **training set (70%)** and a **holdout test set (30%)**. All subsequent feature selection and model training are performed *only* on the training set. The test set is used just once, at the very end, to evaluate the final model's performance on unseen data.
4.  **Handling Class Imbalance**: I utilized Synthetic Minority Over-sampling Technique (SMOTE) within repeated cross-validation folds to mitigate bias towards the majority class.
<img width="607" height="896" alt="class_distribution" src="https://github.com/user-attachments/assets/d5e1a616-c3bf-4c3f-afce-bf64550709bf" />

5.  **Feature Selection (Stability Selection with LASSO)**: To identify a robust and minimal set of genes, I employed LASSO regression combined with stability selection. This process involves repeatedly training LASSO on bootstrapped samples of the *training data* and selecting genes that are consistently identified as important (selected in at least 60% of iterations). This ensures that selected features are not just random noise.
6.  **Model Training**: I trained two powerful classification models – **Random Forest** and **LASSO (Generalized Linear Model with Elastic Net regularization)** – on the training data using only the stable genes identified in the previous step. Hyperparameter tuning is performed via repeated cross-validation.
7.  **Model Evaluation**: Final models are evaluated on the completely unseen holdout test set using metrics like Area Under the Receiver Operating Characteristic curve (AUC) and confusion matrices.
8.  **Overfitting Check & Permutation Testing**: I compared training and test set performance and conducted permutation testing to statistically assess if the observed performance is significantly better than random chance.

## Key Findings & Results

This analysis yielded extremely high performance across both Random Forest and LASSO models, demonstrating an **Area Under the Curve (AUC) of 1.0** on both the training and independent holdout test sets. This indicates a perfect ability of the models to distinguish between IPF patients and healthy controls in this dataset.

### The Phenomenon of AUC = 1.0: A Deeper Look

An AUC of 1.0 in real-world biological data is rare and often prompts scrutiny for data leakage or artifacts. However, after carefully implementing a robust, leakage-proof workflow (including rigorous data splitting and internal cross-validation for feature selection and model training), the results indicate a fundamental property of this specific dataset:

**The identified 'stable genes' exhibit absolute, non-overlapping expression differences between the Disease and Normal groups.** This means that for each of these genes, there is a clear threshold where all healthy samples fall on one side (e.g., lower expression) and all IPF samples fall on the other side (e.g., higher expression). This perfect separation by a few key genes allows the machine learning models to achieve flawless classification.

While this is an exceptionally strong signal, it warrants further investigation:
*   **Biological Validation**: Are these genes known to be absolutely defining markers of IPF with no overlap in expression in healthy individuals in other studies?
*   **External Dataset Validation**: The ultimate test would be to apply this gene signature and model to a completely independent IPF cohort from a different source. If the performance holds, these genes represent highly potent diagnostic and potentially therapeutic targets.

---
