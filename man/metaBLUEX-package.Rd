\name{metaBLUEX-package}
\alias{metaBLUEX-package}
\alias{metaBLUEX}
\docType{package}
\title{
\packageTitle{metaBLUEX}
}
\description{
\packageDescription{metaBLUEX}
}
\author{
\packageAuthor{metaBLUEX}

Maintainer: \packageMaintainer{metaBLUEX}
}
\references{
Wan, X., Wang, W., Liu, J., & Tong, T. (2014). Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. \emph{BMC medical research methodology, 14}, 1-13.

Luo, D., Wan, X., Liu, J., & Tong, T. (2018). Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. \emph{Statistical methods in medical research, 27}(6), 1785-1805.

Shi, J., Luo, D., Weng, H., Zeng, X. T., Lin, L., Chu, H., & Tong, T. (2020). Optimally estimating the sample standard deviation from the five‐number summary. \emph{Research synthesis methods, 11}(5), 641-654.

Balakrishnan, N., Rychtář, J., Taylor, D., & Walter, S. D. (2022). Unified approach to optimal estimation of mean and standard deviation from sample summaries. \emph{Statistical Methods in Medical Research, 31}(11), 2087-2103.

Yang, X., Hutson, A. D., & Wang, D. (2022). A generalized BLUE approach for combining location and scale information in a meta-analysis. \emph{Journal of Applied Statistics, 49}(15), 3846-3867.
}

\keyword{Best linear unbiased estimator (BLUE), Five-number summary, Meta-analysis}
\examples{
library(metaBLUEX)
set.seed(1)
n = 100
# mu = 1, sigma = 1
r = rnorm(n, 1, 1)
fivenumber = fivenum(r)

metablue(fivenumber, n, "S3", "normal")
# 1.0981794 0.9141252

metabluex(fivenumber, n, "S3")
}
