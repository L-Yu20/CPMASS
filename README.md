# Project Overview

This repository includes R and Rcpp code for change point detection experiments based on the MASS framework.

## File Descriptions

### Parameter Selection
- `alpha_choice.R`, `alpha_choiceim.R`, `alpha_choiceout.R`  
  Scripts for selecting the tuning parameter \(\alpha_n\) under different experimental settings.

### Experiment 1: When G follows a Gaussian distribution
- `gus.R`, `gusim.R`, `gusimout.R`  
  Simulation code for Experiment 1.

### Experiment 2: When G follows a t-distribution
- `tf.R`, `tfim.R`, `tfimout.R`  
  Simulation code for Experiment 2.

### Experiment 3: When G follows a Uniform distribution
- `unif.R`, `unifim.R`, `unifimout.R`  
  Simulation code for Experiment 3.

### Real Data Applications
- `macroeconomic.R`, `medical_data.R`  
  Analysis scripts for two real-world datasets: U.S. macroeconomic indicators and genomic medical data.

### Rcpp Source Files: Loss Function Computation
- `F4DDS_binary.cpp`, `F4DDS_pairwise.cpp`, `dF4DDS_binary.cpp`, `dF4DDS_pairwise.cpp`  
  Rcpp implementations for computing the loss function and its derivatives, intended to be called from R.

### Auxiliary Optimization Functions
- `F.R`, `dF.R`, `ybi.R`, `ypa.R`  
  R functions for gradient-based optimization.

### Evaluation Metric
- `randi.R`  
  Computes the Rand Index to assess clustering accuracy.

---

## Datasets

### Macroeconomic Data
- Files: `1.csv`, `4.csv`, `5.csv`, `6.csv` (located in the `macro` folder)  
- Description: U.S. macroeconomic indicators used in the real data experiments.  
- Source: Downloaded from the [Federal Reserve Bank of St. Louis (FRED)](https://fred.stlouisfed.org/)

### Medical Data
- Dataset: `ACGH` (Array Comparative Genomic Hybridization)  
- Description: Genomic data used to detect DNA copy number variations.  
- Source: Available in the R package [`ecp`](https://cran.r-project.org/package=ecp)
