# Time Variation in the Sticky Information Phillips Curve Research

## Overview
This repository contains the files and data related to the research paper titled "Time Variation in the Sticky Information Phillips Curve." The paper explores the dynamics of inflation using a Sticky Information Phillips Curve (SIPC).

## Research Paper
1. **Time Variation in the Sticky Information Phillips Curve.pdf**: This is the main research paper in PDF format. It contains detailed findings, analysis, and conclusions regarding the time variation in the SIPC.

## R Scripts for Estimation
1. **Estimating Fixed and Time Varying Parameters.R**: This R script utilizes the Metropolis-Hastings Sampler to estimate the parameters of both fixed and time-varying models. It employs a joint likelihood approach to account for an endogenous variable.
2. **analyze_data.R**: This R script file contains the code used for the analysis of the data. It includes functions and routines for processing and analyzing the estimations produced by `Estimating Fixed and Time Varying Parameters.R`.
3. **Bayesian_Model_Comparison.R**: This R script performs Bayesian Model Comparison to find the marginal likelihood of each model and calculate the relative probability of each model.

## Data
1. **MyData.Rda**: This R data file contains the processed data used in the analysis. It includes the variables necessary for evaluating the Phillips Curve dynamics under the sticky information framework.
2. **Raw data**: This folder contains preprocessed data that can be found online.

## Output
1. **TVP_Data.Rda**: This R data file contains the output generated from running the `Estimating Fixed and Time Varying Parameters.R` script for Time-Varying Parameters. It includes samples and the acceptance rate from estimation.
2. **Fixed_Data.Rda**: This R data file contains the output generated from running the `Estimating Fixed and Time Varying Parameters.R` script for Fixed Parameters. It includes samples and the acceptance rate from estimation.
3. **attentionbounds.png**: This image file is a visual representation of the time-varying parameter with 95% confidence bounds.

