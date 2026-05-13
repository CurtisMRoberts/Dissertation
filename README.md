#Dissertation code

This repository contains R code used for my dissertation on handling missing values in cross-sectional and longitudinal data analysis.

The repository is organised into three folders:

- 'Chapter2/': R code used to generate the figures in Chapter 2. 
- 'Simulation-Study/': R code for the Chapter 5 simulation study. 
- 'Real-data Application/': R code for the Chapter 5 real-data application using the PRIME colorectal cancer trial. A 


```text
Dissertation/
├── Chapter2/
│   └── Chapter 2 figure code.R
├── Simulation-Study/
│   └── Chapter 5 Simulation Study Code UPDATED.R
└── Real-data Application/
    ├── Real-data application on PRIME dataset.R
    ├── Create Github Demo Subset from PRIME dataset.R
    └── prime_full_demo_50.csv
```


#Data availability

The original PRIME clinical trial dataset is not included in this repository because it is restricted-access clinical trial data.

The file `prime_full_demo_50.csv` is a privacy-preserved demonstration subset. Subject identifiers have been replaced and observed tumour-burden values have been perturbed. This file is provided only to illustrate the data structure and code workflow; hence, it will not reproduce the numerical results reported in the dissertation.
