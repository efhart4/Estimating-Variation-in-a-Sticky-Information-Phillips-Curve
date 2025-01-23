# Time Variation in the Sticky Information Phillips Curve

## Overview
In this research, I investigate the empirical evidence of structural parameter variation within the sticky information model, exploring both the statistical significance and economic implications of changing attention patterns among firms. Using survey data from the survey of professional forecasters, collected by the Philadelphia Federal Reserve, I estimate models capturing firms’ attention to macroeconomic conditions, uncovering a systematic decline in attention over time. This decline is linked to a reduction in the expected volatility of prices. Through Bayesian model comparison, I assess the performance of fixed and time-varying attention models, finding that while a fixed parameter model is favored, the time-varying model still holds relevance.


## Research Paper
1. [**Time Variation in the Sticky Information Phillips Curve**](https://github.com/efhart4/Estimating-Variation-in-a-Sticky-Information-Phillips-Curve/blob/main/Time%20Variation%20in%20the%20Sticky%20Information%20Phillips%20Curve.pdf)

## Data Processing
1. https://html-preview.github.io/?url=https://github.com/efhart4/Estimating-Variation-in-a-Sticky-Information-Phillips-Curve/blob/main/Data/Data-Preparation.html

## R Scripts for Estimation and Interpretation
1. **Metropolis Hasting Sampler.R**: This R script utilizes the Metropolis-Hastings Sampler to estimate the parameters of both fixed and time-varying models. It employs a joint likelihood approach to account for an endogenous variable.
2. **Analyze Samples.R**: This R script contains the code to analyze samples produced by `Metropolis Hasting Sampler.R`.
3. **Bayesian Model Comparison.R**: This R script performs Bayesian Model Comparison of the fixed and time-varying models.

## Keywords
Nonlinear Estimation, Bayesian Statistics, Aggregate Price Model, R
