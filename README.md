# metaBLUEX <img src="https://raw.githubusercontent.com/SungManhin/metaBLUEX/refs/heads/main/man/figures/logo.png" alt="metaBLUEX" height="150px" align="right" />


## Introduction

In meta-analysis for continuous outcomes, the classical methods are designed to handle the studies that report the sample mean and standard deviation. However, some researchers choose to report other summary statistics, e.g. quartiles, mean, min/max. 

To effectively utilize these summary statistics in meta-analysis, numerous methods have been developed for estimating the mean and standard deviation. Among these, the BLUE estimators based on generalized least squares have garnered significant attention due to their favorable theoretical properties. 

This package provides implementations of BLUE estimators for the mean and standard deviation across different scenarios of summary statistics and distributional assumptions. The methods included are as follows: under normality assumption, approximation formulae for BLUE estimators proposed by [Wan et al. (2014)](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-14-135), [Luo et al. (2018)](https://journals.sagepub.com/doi/10.1177/0962280216669183)  and [Shi et al. (2020)](https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1429), approximation for the expectation vector and variance-covariance matrix in BLUE estimators proposed by [Yang et al. (2022)](https://www.tandfonline.com/doi/full/10.1080/02664763.2021.1967890) and [Balakrishnan et al. (2022)](https://journals.sagepub.com/doi/10.1177/09622802221111546); under other location-scale family distributions: approximations for BLUE estimators specifically for the logistic and Laplace distributions (developed by the authors of this package).

> The logo consists of a blue background, representing the **BLUE** method, and forest plot, symbolizing **meta-analysis**.

## Installation 

You can install this package by running the following code in your `R` console:

```R
# R package `crayon` is needed
# install.packages("crayon")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("SungManhin/metaBLUEX")
```

## Methods


## Examples

```R
library(metaBLUEX)

### normal in S3 ###
set.seed(1)
n = 100
# mu = 1, sigma = 1
r = rnorm(n, 1, 1)
fivenumber = fivenum(r)

# Wan, Luo and Shi's method
metablue(fivenumber, n, "S3", "normal")
# 1.0981794 0.9141252

# Yang's method
metablue(fivenumber, n, "S3", "normal", "yang")
# 1.0857491 0.9065274

# Balakrishnan's method
metablue(fivenumber, n, "S3", "normal", "bala")
# 1.092199 0.933297
```

```R
library(VGAM)

### laplace in S3 ###
set.seed(1)
n = 100
# mu = 1, sigma = 2
# Laplace dist. with parameters (0, 1/sqrt(2)) has mean 0 and sigma 1
l = 1 + 2 * rlaplace(n, 0, 1/sqrt(2))
fivenumber_l = fivenum(l)

# Yang's method
metablue(fivenumber_l, n, "S3", "laplace", "yang")
# 0.9505243 1.7241492

# Our approximation method
metablue(fivenumber_l, n, "S3", "laplace", "approx")
#  0.9887838 1.7481349
```

## References

1. Wan, X., Wang, W., Liu, J., & Tong, T. (2014). Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. *BMC medical research methodology, 14*, 1-13.

2. Luo, D., Wan, X., Liu, J., & Tong, T. (2018). Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. *Statistical methods in medical research, 27*(6), 1785-1805.

3. Shi, J., Luo, D., Weng, H., Zeng, X. T., Lin, L., Chu, H., & Tong, T. (2020). Optimally estimating the sample standard deviation from the five‐number summary. *Research synthesis methods, 11*(5), 641-654.

4. Balakrishnan, N., Rychtář, J., Taylor, D., & Walter, S. D. (2022). Unified approach to optimal estimation of mean and standard deviation from sample summaries. *Statistical Methods in Medical Research, 31*(11), 2087-2103.

5. Yang, X., Hutson, A. D., & Wang, D. (2022). A generalized BLUE approach for combining location and scale information in a meta-analysis. *Journal of Applied Statistics, 49*(15), 3846-3867.