
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="data/figure/Title.png" width="100%" style="margin-top: 20px;" />

<!-- badges: start -->

[![R
\>3.5](https://img.shields.io/badge/R-%3E3.5-success.svg)](https://www.r-project.org/)
<a href='#devtools'>![installed with
devtools](https://img.shields.io/badge/installed%20with-devtools-blueviolet.svg)</a>
[![WebServer](https://img.shields.io/badge/Web_Server-MMEASE-blue)](http://idrblab.org/mmease2025/)
<!-- badges: end -->

# How to Use `mmease`?

## Contents

- [Overview](#overview)
- [Installation](#installation)
- [Example Data](#example-data)
- [Usage and Examples](#usage-and-examples)

## Overview

Metabolomics closest to the phenotype is shifting to the single-cell
level (SCM), which is a powerful tool for studying cellular
heterogeneity by providing insight into the differences be-tween
individual cells. Because analytical workflow of single-cell
metabolomics (SCM) in-cludes various processes, a standardized and
transparent pipeline is essential to connect the raw data to biological
interpretation especially for non-bioinformatic researchers. However, a
uni-fied tool to provide the entire workflow of SCM is still
lacking.<br><br>In this study, MMEASE was updated to its 2.0 version by
advancing from bulk metabolomics to the single-cell metabolom-ics, which
realize a comprehensive analytical workflow of SCM for the first time.
Specifically, MMEASE 2.0 is a unique tool in SCM by (a) providing the
most comprehensive workflow of data processing, (b) realizing
systematical functions of analyzing both metabolic heterogeneity and
functional heterogeneity, and (c) significantly enhancing the capability
of metabolite an-notation with biological interpretations using tandem
spectra.<br><br>In summary, MMEASE 2.0 pro-vides an indispensable online
service of whole analytical workflow in SCM. MMEASE 2.0 is freely
accessible at
<a href="https://idrblab.org/mmease2025/">http://idrblab.org/mmease2025</a>.

## Installation

Install a variety of *R* packages imported in this protocol.

Installed from ***CRAN*** (can also from other
repositories):dummies,e1071,adabag,C50,pROC,kknn,MASS,multiROC,caret,mlbench,randomForest,varSelRF,magrittr.

Installation commands:

``` r
CRAN_packages <- c("dummies","e1071","adabag","C50","pROC","kknn","MASS","multiROC","caret","mlbench","randomForest","varSelRF","magrittr")
install.packages(CRAN_packages, dependencies = TRUE)
```

Installed from ***Bioconductor*** (can also from other repositories):
ropls,mixOmics,multtest,AUC,metabolomics.

Installation commands:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

Bioconductor_packages <- c("multtest","mixOmics","ropls","AUC","metabolomics")

BiocManager::install(Bioconductor_packages, ask = FALSE)
```

Install the `mmease` package by running the following command in
RStudio:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
library(devtools)

devtools::install_github("mmease2025/mmease")
```

## Example Data

You can use the example data provided in `mmease` to try it out.

**Example single-cell metabolomic data**
([Download](https://iboalab.cn/mmease2025/MTBLS78_co_culture.csv)): As
shown in the sample data, cell name, cell classes, cell types, and batch
information are required in the first four columns of the input file. In
the following columns, the raw peak intensities across all cells are
further provided. Unique metabolite IDs or peaks are listed in the first
row of the csv file.

## Usage and Examples

Data Filtering

``` r
data <- read.csv("/yourpath/SCMdata.csv")
filted_data <- filtering(data,Percent of Missing Values)
#Percent of Missing Values must be between 0 and 1
```

Data Imputation

``` r
imputed_data <- imputation(filted_data, method = "1/5 of minimum positive value")
imputed_data <- imputation(filted_data, method = "KNN")
```

Data Transformation

``` r
transformed_data <- transformation(imputed_data, method = "G-log")
transformed_data <- transformation(imputed_data, method = "Log2")
transformed_data <- transformation(imputed_data, method = "Log10")
```

Data normalization

``` r
normalized_data <- normalization(transformed_data, method = "Auto Scaling")
normalized_data <- normalization(transformed_data, method = "Mean")
normalized_data <- normalization(transformed_data, method = "Median")
normalized_data <- normalization(transformed_data, method = "MSTUS")
normalized_data <- normalization(transformed_data, method = "SIS")
```

Batch correction

``` r
corrected_data <- batch_correction(normalized_data, method = "ComBat")
corrected_data <- batch_correction(normalized_data, method = "Limma")
```

DifferentialAnalysis

``` r
#Functional analysis
functional_table <- functional_differential(corrected_data,method = "t_test")
functional_table <- functional_differential(corrected_data,method = "ANOVA")
functional_table <- functional_differential(corrected_data,method = "FC")
functional_table <- functional_differential(corrected_data,method = "PLS-DA")
functional_table <- functional_differential(corrected_data,method = "OPLS-DA")
functional_table <- functional_differential(corrected_data,method = "RF_RFE")
functional_table <- functional_differential(corrected_data,method = "Kruskal_Wallis")
functional_table <- functional_differential(corrected_data,method = "svmrfeFeatureRanking")


#Metabolic analysis
metabolic_table <- metabolic_differential(corrected_data,method = "t_test")
metabolic_table <- metabolic_differential(corrected_data,method = "ANOVA")
metabolic_table <- metabolic_differential(corrected_data,method = "FC")
metabolic_table <- metabolic_differential(corrected_data,method = "PLS-DA")
metabolic_table <- metabolic_differential(corrected_data,method = "OPLS-DA")
metabolic_table <- metabolic_differential(corrected_data,method = "RF_RFE")
metabolic_table <- metabolic_differential(corrected_data,method = "Kruskal_Wallis")
metabolic_table <- metabolic_differential(corrected_data,method = "svmrfeFeatureRanking")
```

Classification model

``` r
#Functional analysis
functional_class <- functional_classification(corrected_data,method = "AdaBoost")
functional_class <- functional_classification(corrected_data,method = "Bagging")
functional_class <- functional_classification(corrected_data,method = "Decision Trees")
functional_class <- functional_classification(corrected_data,method = "K-Nearest Neighbor")
functional_class <- functional_classification(corrected_data,method = "Linear Discriminat Analysis")
functional_class <- functional_classification(corrected_data,method = "Naive Bayes")
functional_class <- functional_classification(corrected_data,method = "Partial Least Squares")
functional_class <- functional_classification(corrected_data,method = "Random Forest")
functional_class <- functional_classification(corrected_data,method = "Support Vector Machine")


#Metabolic analysis
metabolic_class <- metabolic_classification(corrected_data,method = "AdaBoost")
metabolic_class <- metabolic_classification(corrected_data,method = "Bagging")
metabolic_class <- metabolic_classification(corrected_data,method = "Decision Trees")
metabolic_class <- metabolic_classification(corrected_data,method = "K-Nearest Neighbor")
metabolic_class <- metabolic_classification(corrected_data,method = "Linear Discriminat Analysis")
metabolic_class <- metabolic_classification(corrected_data,method = "Naive Bayes")
metabolic_class <- metabolic_classification(corrected_data,method = "Partial Least Squares")
metabolic_class <- metabolic_classification(corrected_data,method = "Random Forest")
metabolic_class <- metabolic_classification(corrected_data,method = "Support Vector Machine")
```

Plot model

``` r

plots <- classification_plots(functional_class)
```
