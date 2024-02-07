# simbasilica

R package to generate synthetic data from the Basilica generative model.

The synthetic validation is performed on 270 cohorts of patients:
  * Number of patients: 150, 500, 1000
  * Number of groups: 1, 3, 6
  * 30 datasets for each configuration
  

Each dataset has been fitted with Basilica:
  * Datasets with sample sizes 150 and 500 -> CPU
  * Datasets with sample size 1000 -> GPU 

