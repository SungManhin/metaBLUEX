% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metaBLUEX.R
\name{metablue}
\alias{metablue}
\title{BLUE Estimators for Mean and Standard Deviation}
\usage{
metablue(
  summary = NULL,
  n = NULL,
  scenario = "S1",
  dist = "normal",
  method = "wls"
)
}
\arguments{
\item{summary}{a vector of ordered summary statistics in ascending order.}

\item{n}{the sample size.}

\item{scenario}{a character string indicating the scenario in which summary statistics is reported. The options for the \code{scenario} argument are as follows.
\tabular{ll}{
\code{"S1"} \tab for the median, minimum and maximum values. \cr
\code{"S2"} \tab for the median, first and third quartiles. \cr
\code{"S3"} \tab for the median, first and third quartiles, and minimum and maximum values. \cr
\code{"tertiles"} \tab for tertiles. \cr
\code{"quintiles"} \tab for quintiles. \cr
\code{"deciles"} \tab for deciles. \cr}}

\item{dist}{a character string indicating which distribution of data is assumed. The options for the \code{dist} argument are as follows.
\tabular{ll}{
\code{"normal"} \tab for the normal distribution. \cr
\code{"laplace"} \tab for the Laplace distribution. \cr
\code{"logistic"} \tab for the logistic distribution. \cr}}

\item{method}{a character string indicating which method to use. The options for the \code{method} argument are as follows.
\tabular{ll}{
\code{"wls"} \tab an abbreviation for approximation formulae of BLUE estimators proposed by Wan et al. (2014), Luo et al. (2018) and Shi et al. (2020). Only applicable when \code{dist} is set to \code{"normal"} and \code{scenario} is set to \code{"S1"}, \code{"S2"} or \code{"S3"}. \cr
\code{"bala"} \tab Balakrishnan et al. (2022)'s approximation method. Only applicable when \code{dist} is set to \code{"normal"}. \cr
\code{"yang"} \tab Yang et al. (2022)'s approximation method. \cr
\code{"approx"} \tab proposed by the authors of this package to handle the Laplace and logistic distribution. \cr}}
}
\value{
A vector of length two, the first element is \eqn{\hat\mu} and the second one is \eqn{\hat\sigma}.
}
\description{
A comprehensive function for estimating the mean and standard deviation from summary statistics, incorporating both published and novel BLUE methods.
}
\details{
\code{metablue} implements BLUE methods for estimating the mean and standard deviation from summary statistics,
including methods proposed by the authors.
Below are important notes regarding the function's arguments.
\itemize{
   \item \code{method = "wls"} can only be used under normality assumption and three basic scenarios.
   \item \code{method = "bala"} can only be used under normality assumption. Balakrishnan et al. (2022)'s approximation method is proposed for three basic scenarios only, and the authors of this package extend their method to tertiles, quintiles and deciles. Additional message will be printed when using the method in extended scenarios.
   \item \code{method = "yang"} can be used for all six scenarios under normality assumption. The authors of this package extend their method to the Laplace and logistic distribution. Additional message will be printed when using the method in extended distributions.
   \item \code{method = "approx"} is proposed for the Laplace and logistic distribution and three basic scenarios, inspired by Wan et al. (2014), Luo et al. (2018) and Shi et al. (2020)'s works, check \url{https://github.com/sungmanhin/metaBLUEX} for approximation details.
}
}
\examples{
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

### laplace in S3 ###
library(VGAM)
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

}
\references{
Wan, X., Wang, W., Liu, J., & Tong, T. (2014). Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. \emph{BMC medical research methodology, 14}, 1-13.

Luo, D., Wan, X., Liu, J., & Tong, T. (2018). Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. \emph{Statistical methods in medical research, 27}(6), 1785-1805.

Shi, J., Luo, D., Weng, H., Zeng, X. T., Lin, L., Chu, H., & Tong, T. (2020). Optimally estimating the sample standard deviation from the five‐number summary. \emph{Research synthesis methods, 11}(5), 641-654.

Balakrishnan, N., Rychtář, J., Taylor, D., & Walter, S. D. (2022). Unified approach to optimal estimation of mean and standard deviation from sample summaries. \emph{Statistical Methods in Medical Research, 31}(11), 2087-2103.

Yang, X., Hutson, A. D., & Wang, D. (2022). A generalized BLUE approach for combining location and scale information in a meta-analysis. \emph{Journal of Applied Statistics, 49}(15), 3846-3867.
}
\author{
Sung Manhin \email{songwenxuan@ruc.edu.cn}
}
