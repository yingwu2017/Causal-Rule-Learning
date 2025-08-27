This repository contains the code and resources for the simulation
studies presented in the paper **“A New Causal Rule Learning Approach to
Interpretable Estimation of Heterogeneous Treatment Effect”**.

The project implements a three-stage framework (Rule Discovery, Rule
Selection, and Rule Analysis) for estimating Heterogeneous Treatment
Effects (HTE) using an interpretable, rule-based approach.

Please refer to the [CRL paper on
arxiv](https://arxiv.org/pdf/2310.06746) for details of this project.

## Project Structure

    Causal_Rule_Learning/          # Project Root (RStudio Project)
    │
    ├── Causal_Rule_Learning.Rproj # RStudio Project file
    │
    ├── functions/                 # Directory for core R functions
    │   ├── CRL_functions.R       # Main functions for the CRL workflow
    │   └── rule.process.R        # Utilities for rule processing & analysis
    │
    ├── simulation/               # Directory for all simulation studies
    │   │
    │   ├── study one/            # Simulation Study 1
    │   │   ├── CRL/              # Directory for scripts of CRL workflow
    │   │   │   ├── CRL_study_one.R   # Main script to run study
    │   │   │   └── CRL_evaluation.R  # Script for evaluating CRL results
    │   │   ├── model.result      # Folder containing saved CRL outputs
    │   │   └── data/             # Folder containing simulated data for study one
    │   │       └── data_generation.R # Script for generating data for study one
    │   │
    │   └── study two/            # Simulation Study 2
    │       ├── data/             # Folder containing simulated data for study two
    │       │   └── data_generation.R # Script for generating data for study two
    │       └── study_two.R       # Main script to run study two
    │
    └── .Rproj.user/              # RStudio IDE settings (auto-generated)

## Prerequisites & Installation

1.  **R**: Ensure you have R (version 4.1) or later installed.

2.  **RStudio**: (Recommended) The project is configured as an RStudio
    project.

3.  **R Packages**: This project depends on the following key R
    packages. Please install them before proceeding:

    `install.packages(c("grf", "glmnet", "dplyr", "ggplot2"))`

## Usage

To reproduce the results for a specific study, open the
`Causal_Rule_Learning.Rproj` file in RStudio in which you can run the
codes without changing working directory.

You have to first generate the simulated data running
**data\_generation.R** script and the data will be generated and saved
in the *data/*.

1.  For study one, run the **CRL\_study\_one.R** to implement CRL on the
    simulated data and the results will be saved in *model.result/ .*
    Run **CRL\_evaluation.R** to evaluate model results and obtain the
    performance metric values.

2.  For study two, you can implement and evaluate CRL by running
    **study\_two.R**.
