Files used for analysis in the SPEAR manuscript. Organized into simulations, TCGA-BC, and COVID-19.

The COVID-19 analysis can be run locally, as the data size is not too large. The TCGA-BC analysis was done using a computing cluster via jobs, with each job running one model given the size of the dataset. As such, a job file is included for reference in the TCGA-BC analysis and not in the COVID-19 analysis.

Each file begins with a python path requirement (for MOFA+), save path (to save the results), and load data path (to point to the data).
