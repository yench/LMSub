# LMSub
R code for "Dynamic prediction by landmarking with data from cohort subsampling designs" by Yen Chang, Anastasia Ivanova, Jason P. Fine, Demetrius Albanes, and Yei Eun Shin.

## Simulated data sets
Examples of data used in the simulation study (Section 3). Each subject has one row for each landmark time he/she is still at risk.
- fuldat.csv: Full cohort data for training. This dataset is used to run full cohort and IPW (with subsetting) landmark analyses. 
- nccdat.csv: NCC data for training. This dataset is used to run conditional analyses. 
- valdat.csv: validation data.

## Landmarking analyses
- LMSub_analysis.R: code for dynamic prediction by landmark analyses with NCC and full cohort data.
- LMSub_functions.R: R functions for influence-function based variance estimation. Sourced by LMSub_analysis.R.
