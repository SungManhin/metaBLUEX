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
install.packages("crayon")
if(!require(devtools)) install.packages("devtools")
devtools::install_github("SungManhin/metaBLUEX")
```

## Methods

### Wan, Luo and Shi's method

$$
\hat{\mu} = w_{1,\text{opt}} \left(\frac{a+b}{2}\right) + w_{2,\text{opt}} \left(\frac{q_1+q_3}{2}\right) +\left(1 - w_{1,\text{opt}} - w_{2,\text{opt}}\right) m
$$

$$
\hat{\sigma} = w_{3,\text{opt}} \left(\frac{b-a}{2 \mathrm{E}(Z_{(n)})}\right) + 
\left(1 - w_{3,\text{opt}}\right) \left(\frac{q_3 - q_1}{2 \mathrm{E}(Z_{(3Q+1)})}\right)
$$

| Weights / Scenarios | $\widetilde{w}_{1,\text{opt}}$ | $\widetilde{w}_{2,\text{opt}}$ | $\widetilde{w}_{3,\text{opt}}$ |
| :-----: | :------------------------------: | :-------------------------------: | :------------------------------: |
| $S_1$ | $$\frac{4}{4 + n^{0.75}}$$   | $$0$$                         | $$1$$                        |
| $S_2$ | $$0$$                        | $$0.7 + \frac{0.39}{n}$$      | $$0$$                        |
| $S_3$ | $$\frac{2.2}{2.2 + n^{0.75}}$$ | $$0.7 - \frac{0.72}{n^{0.55}}$$ | $$\frac{1}{1 + 0.07 n^{0.6}}$$ |

$2 \mathrm{E}(Z_{(n)})=2\Phi^{-1}[(n-0.375)/(n+0.25)]$, $2 \mathrm{E}(Z_{(n)})=2\Phi^{-1}[(0.75n-0.125)/(n+0.25)]$.

### Yang's method

Let $Z_i\sim\mathcal{N}(0,1)$, $i=1,\dots,n$, and $Z_{(1)},\dots,Z_{(n)}$ be the order statistics of $Z_{1},\dots,Z_{n}$ with sample size $n$.

$$
\mathrm{E}\left(Y_{(r)}\right)= Q_r+\frac{p_r q_r}{2(n+2)} Q_r^{\prime \prime}+\frac{p_r q_r}{(n+2)^2}\left[\frac{1}{3}\left(q_r-p_r\right) Q_r^{\prime \prime \prime}+\frac{1}{8} p_r q_r Q_r^{\prime \prime \prime \prime}\right],
$$

$$
\mathrm{Var}\left(Y_{(r)}\right)= \frac{p_r q_r}{n+2} Q_r^{\prime 2}+\frac{p_r q_r}{(n+2)^2}\left[2\left(q_r-p_r\right) Q_r^{\prime} Q_r^{\prime \prime}+p_r q_r\left(Q_r^{\prime} Q_r^{\prime \prime \prime}+\frac{1}{2} Q_r^{\prime \prime 2}\right)\right],
$$

$$
\mathrm{Cov}\left(Y_{(r)}, Y_{(s)}\right)= \frac{p_r q_s}{n+2} Q_r^{\prime} Q_s^{\prime}+\frac{p_r q_s}{(n+2)^2}\left[\left(q_r-p_r\right) Q_r^{\prime \prime} Q_s^{\prime}+\left(q_s-p_s\right) Q_r^{\prime} Q_s^{\prime \prime}\right. 
\left.+\frac{1}{2} p_r q_r Q_r^{\prime \prime \prime} Q_s^{\prime}+\frac{1}{2} p_s q_s Q_s^{\prime \prime \prime} Q_r^{\prime}+\frac{1}{2} p_r q_s Q_r^{\prime \prime} Q_s^{\prime \prime}\right],
$$

where $Q_r = Q(p_r)$, $p_r = r /(n+1)$ , $q_r = 1 − p_r$, and $Q^{′}_r$, $Q^{′′}_r$, $Q^{′′′}_r$, $Q^{′′′′}_r $ are the first four derivatives of the quantile function $Q(u)$ at $u = p_r$. For the most commonly used standard normal  distribution, in particular, $Q^′_r = \sqrt{2\pi} e^{\gamma^2}$ , $Q^{′'}_r = −2\sqrt{2\pi} \gamma e^{2\gamma^2}$ , $Q^{'''}_r = 2\sqrt{2}\pi^{3/2}e^{3\gamma^2} (1 +  4\gamma^2)$, $Q^{''''}_r = −4\sqrt{2}\pi^2\gamma e^{4\gamma^ 2} (7 + 12\gamma^2)$, where $\gamma=-\Phi^{-1}(\frac{r}{n+1})/\sqrt{2}$​.

For the Laplace distribution $f_{\text{Lap}}(x)=\frac{1}{2b}\exp\left(-|x-\mu|/b\right)$, its quantile function is $F^{-1}_{\text{Lap}}(p)=\mu-b\ \mathrm{sgn}(p-0.5)\ln(1-2|p-0.5|)$; for the Logistic distribution $f_{\text{Logit}}(x)=\frac{\exp(-(x-\mu)/s)}{s(1+\exp(-(x-\mu)/s))^2}$, its quantile function is $F^{-1}_{\text{Logit}}(p)=\mu+s\ln\left({p}/{(1-p)}\right)$.

### Balakrishnan's method

Let $\boldsymbol{B}=\{\beta_{i,j:n}\}$ be the variance-covariance matrix of  $Z_{(1)},\dots,Z_{(n)}$.

$$
\beta_{r, s: n}^{1 / 4} \approx a_{\frac{r}{n}, \frac{s}{n}} \ln (\ln n)+b_{\frac{r}{n}, \frac{s}{n}},
$$

$$
a_{q_1, q_2} \approx \frac{1}{8} \sqrt{1-4\left(q_1-0.5\right)\left(q_2-0.5\right)+2\left|q_1-q_2\right|^{0.5}}-0.51, 
$$

$$
b_{q_1, q_2} \approx \frac{1}{2.8} \sqrt{1-4\left(q_1-0.5\right)\left(q_2-0.5\right)+2\left|q_1-q_2\right|^{0.5}}+1.3.
$$

For $\beta_{1,1:n}=\beta_{n,n:n}$, use $\beta^{1/4}_{1,1:n}=\beta^{1/4}_{n,n:n}=−0.1565 \ln (\ln n) + 0.8949$ instead due to a different order of min/max order statistics's variaince.

Expectation of $Z_{(1)},\dots,Z_{(n)}$  is derived similar to Yang's method.

### Our approximation for Laplace and Logistic

**Laplace**

$$
\hat{\mu}=w^{\prime}_{1,\mathrm{opt}}\left(\frac{a+b}{2}\right)+w^{\prime}_{2,\mathrm{opt}}\left(\frac{q_{1}+q_{3}}{2}\right)+(1-w^{\prime}_{1,\mathrm{opt}}-w^{\prime}_{2,\mathrm{opt}})m
$$

$$
\hat{\sigma}=w^{\prime}_{3,\mathrm{opt}}\left(\frac{b-a}{2/(-0.008+1.462(\log(x))^{-0.988})}\right)+(1-w^{\prime}_{3,\mathrm{opt}})\left(\frac{q_3-q_1}{2/(2.042+2.919x^{-1.209})}\right)
$$

| Weights / Scenarios | $\widetilde{w}_{1,\text{opt}}^{\prime}$  | $\widetilde{w}_{2,\text{opt}}^{\prime}$ |           $\widetilde{w}_{3,\text{opt}}^{\prime}$           |
| :-----------------: | :--------------------------------------: | :-------------------------------------: | :---------------------------------------------------------: |
|        $S_1$        |  $$\frac{1}{-0.469 + 0.599 x^{1.606}}$$  |                  $$0$$                  |                            $$1$$                            |
|        $S_2$        |                  $$0$$                   | $$\frac{1}{-0.045 + 0.871 x^{0.555}}$$  |                            $$0$$                            |
|        $S_3$        | $$\frac{1.122}{\exp(-1.137 - 0.476 x)}$$ | $$\frac{1}{0.297 + 0.759  x^{0.577}}$$  | $$\frac{9.306}{10.870 + 0.065x + 0.581 (\log(x))^{1.359}}$$ |

**Logistic**

$$
\hat{\mu}=w^{\prime\prime}_{1,\mathrm{opt}}\left(\frac{a+b}{2}\right)+w^{\prime\prime}_{2,\mathrm{opt}}\left(\frac{q_{1}+q_{3}}{2}\right)+(1-w^{\prime\prime}_{1,\mathrm{opt}}-w^{\prime\prime}_{2,\mathrm{opt}})m
$$

$$
\hat{\sigma}=w^{\prime\prime}_{3,\mathrm{opt}}\left(\frac{b-a}{2/(0.183+1.198x^{-0.426})}\right)+(1-w^{\prime\prime}_{3,\mathrm{opt}})\left(\frac{q_3-q_1}{2/(1.653+3,185x^{-1.132})}\right)
$$

| Weights / Scenarios | $\widetilde{w}_{1,\text{opt}}^{\prime\prime}$ | $\widetilde{w}_{2,\text{opt}}^{\prime\prime}$ |    $\widetilde{w}_{3,\text{opt}}^{\prime\prime}$     |
| :-----------------: | :-------------------------------------------: | :-------------------------------------------: | :--------------------------------------------------: |
|        $S_1$        |   $$\frac{65.452}{50.318+27.929x^{0.993}}$$   |                     $$0$$                     |                        $$1$$                         |
|        $S_2$        |                     $$0$$                     |      $$\frac{1}{1.667-1.057x^{-0.997}}$$      |                        $$0$$                         |
|        $S_3$        | $$\frac{109.931}{-17.958+117.682x^{0.992}}$$  |          $$\frac{0.605x}{0.991+x}$$           | $$\frac{448.738}{404.359x^{0.143}+16.884x^{0.803}}$$ |

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