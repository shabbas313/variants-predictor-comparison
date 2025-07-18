---

## Manuscript Figures

This directory contains the final figures presented in the manuscript. Below is a detailed description of each figure.

### Figure 1: Study Design Workflow
*A schematic overview of the study design for the comprehensive evaluation of in-silico pathogenicity predictors. The workflow illustrates the curation of six distinct benchmark datasets from four primary sources (ClinVar, gnomAD, SGE study, and CFTR2) and the subsequent stages of model selection, performance evaluation, and comparative analysis.*

### Figure 2: Performance Evaluation on ClinVar Datasets
*Performance evaluation of 29 models on the ClinVar benchmark datasets.*
-   **(A)** Distribution of pathogenic (21,718) and benign (30,173) variants in the primary ClinVar dataset, with a bar chart showing the percentage of missing predictions for each model.
-   **(B)** AU-ROC on the primary dataset, sorted by performance.
-   **(C)** AU-ROC on the balanced subset.
-   **(D)** Distribution and missing prediction rates for the independent ClinVar dataset.
-   **(E)** AU-ROC on the independent ClinVar dataset, demonstrating model generalization.
-   **(F)** Performance of six representative models stratified by variant allele frequency, showing the impact of rarity on predictive accuracy.

### Figure 3: Performance by Model Category and Diagnostic Utility
*Comparative performance by model category and analysis of clinical utility trade-offs.*
-   **(A)** A box and strip plot illustrating the distribution of AU-ROC scores for each model across the ClinVar datasets. Boxes are colored by model category, and individual points are colored by dataset, showing both overall performance and consistency.
-   **(B-D)** Scatter plots of sensitivity versus specificity for all 29 models on the primary, balanced, and independent ClinVar datasets, respectively.

### Figure 4: Disparity Analysis: Clinical vs. Functional Prediction
*Comparative performance of models on *BRCA1/2* clinical versus *BRCA2* functional (SGE) datasets.*
-   **(A)** AU-ROC and **(B)** Matthews Correlation Coefficient (MCC) for all models on two distinct benchmark sets: variants from *BRCA1* and *BRCA2* with clinical annotations in gnomAD, and variants from *BRCA2* with functional classifications from a Saturation Genome Editing (SGE) study. Models are sorted by their performance on the SGE dataset. The drastic reduction in performance highlights the challenge of predicting functional impact compared to clinical pathogenicity. Dashed lines indicate the performance of a random classifier (AU-ROC = 0.5; MCC = 0).

---
