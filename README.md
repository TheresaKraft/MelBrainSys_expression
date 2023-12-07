# MelBrainSys_expression
Script for patient-matched melanoma metastases expression analysis


# Kraft et. al 2024, Workflow

**Personalized identification and characterization of genome-wide gene expression differences between patient-matched intracranial and extracranial melanoma metastasis pairs**

This is the pipeline that we used for the analysis of patient-matched RNA-seq epression data of melanoma metastases.

# Instructions

This repository has five parts: 

* **MainScript.R**: The main script to use is MainScript.R . It includes all analysis and can be run in R. 
* **HelperFunctions.R**: All functions that are used in the main script are stored in the separate file, HelperFunctions.R. Please load all these helper functions to your R workspace before applying the pipeline.
* **Folder: Annotations**: Contains all functional annotations that are used by the pipeline to analyze the expression data.
* **Folder: Figures**: This folder contains all scripts to create the main figures that are part of the manuscript. A separate folder is used for each figure. Each folder contains an .RData file that comprises all data required to create the specific figure. These R.Data files are created in the pipeline (see MainScript.R). Each folder additionally contains a separate R-script that creates the figure in PNG format.
* **Folder: Programs**: This folder contains the Java JAR file implemetation of the HMM by Seifert et. al, 2014 that was used to make the predictions for the expression states of individual genes (see MainScript.R).


# License

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
