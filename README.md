# Performance Comparison of *in-silico* models for Variant Pathogenicity Prediction
## Abstract 
This repository contains results, figures and necessary codes for our manuscript, "Evaluating the Performance of Variant Effect Predictors: A Comprehensive Analysis Utilizing Clinical and High-throughput Functional Data". In this study, we performed a rigorous, multi-faceted performance validation of 29 widely used *in-silico* pathogenicity predictors. Our findings identify the most reliable predictors for clinical application but also uncover a significant disparity in model performance between clinical and functional benchmarks, providing a vital guide for the informed application of these powerful models in advanced genomic medicine. 

**Link to Published Paper:** [To be added upon publication- In-progress]

## Benchmark Datasets 
Due to their large file size, the benchmark datasets analyzed in this study are not hosted on this GitHub repository. The full datasets are available from the author upon reasonable request.  
Upon receipt, will be provided with the datasets in .csv format. In all files, the ground truth labels are located in a column named Label, where 1 represents Pathogenic and 0 represents Benign. 
To request access, please contact [Syed Hassan Abbas] at [abbas.hassan@stu.ecnu.edu.cn].

- **ClinVar Datasets:** Derived from the ClinVar database (accessed [December 2024]). Includes a primary imbalanced set, a balanced subset, and a temporal validation set containing variants submitted between Jan-May 2025 to ensure non-circularity.
- **BRCA1 and BRACA2 Clinical Dataset:** Curated from gnomAD v4.1.0, containing missense variants in *BRCA1* and *BRCA2* with high-confidence clinical classifications.
- **BRCA2 SGE Functional Dataset:** A large-scale experimental dataset derived from a saturation genome editing study of *BRCA2* [1].
- **CFTR2 Disease-Specific Dataset:** A high-confidence set of CF-causing and non-CF-causing variants from the CFTR2 database.

## Citation 
1. Sahu S, Galloux M, Southon E, Caylor D, Sullivan T, Arnaudi M, et al. Saturation genome editing-based clinical classification of BRCA2 variants. Nature. 2025;638(8050):538-45.

## Contact 
For questions about the paper or data, please contact [Syed Hassan Abbas] at [abbas.hassan@stu.ecnu.edu.cn] or open an issue in this repository.
