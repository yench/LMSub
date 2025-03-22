# LMSub
R code for "Dynamic prediction by landmarking with data from cohort subsampling designs" by Yen Chang, Anastasia Ivanova, Jason P. Fine, Demetrius Albanes, and Yei Eun Shin.

## Simulated data sets
Examples of data used in the simulation study (Section 3). Each subject has one row for each landmark time he/she is still at risk.
- fuldat.csv: Full cohort data for training. This dataset is used to run full cohort and IPW (with subsetting) landmark analyses. 
- nccdat.csv: NCC data for training. This dataset is used to run conditional likelihood analyses. 
- valdat.csv: validation data.

## Landmarking analyses
- LMSub_analysis.R: R code for dynamic prediction by landmarking with NCC and full cohort data.
- LMSub_functions.R: R functions for influence function (IF)-based variance estimation. Sourced by LMSub_analysis.R.
  Parts of the code in the above two files are adapted from the R code accompanying Wenhao Li, Liang Li, and Brad C. Astor (2023) "A comparison of two approaches to dynamic prediction: Joint modeling and landmark  
modeling" published by Statistics in Medicine (https://github.com/liwh0904/Compare_JM_and_LM).
- Rho.cpp: C++ function to compute joint inclusion probabilities needed in IF-based variance estimation.
- dLambda0.cpp: C++ function to compute estimated increments in cumulative baseline hazards and their IFs.
