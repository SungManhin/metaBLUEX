#' @title BLUE Estimators for Mean and Standard Deviation
#' @description A comprehensive function for estimating the mean and standard deviation from summary statistics, incorporating both published and novel BLUE methods.
#' 
#' @param summary a vector of ordered summary statistics in ascending order.
#' @param n the sample size.
#' @param scenario a character string indicating which scenario of summary statistics is reported. The options for the \code{scenario} argument are as follows. 
#' \tabular{ll}{
#' \code{"S1"} \tab for the median, minimum and maximum values. \cr
#' \code{"S2"} \tab for the median, first and third quartiles. \cr
#' \code{"S3"} \tab for the median, first and third quartiles, and minimum and maximum values. \cr
#' \code{"tertiles"} \tab for tertiles. \cr
#' \code{"quintiles"} \tab for quintiles. \cr
#' \code{"deciles"} \tab for deciles. \cr}
#' @param dist a character string indicating which distribution of data is assumed. The options for the \code{dist} argument are as follows.
#' \tabular{ll}{
#' \code{"normal"} \tab for the normal distribution. \cr
#' \code{"laplace"} \tab for the Laplace distribution. \cr
#' \code{"logistic"} \tab for the logistic distribution. \cr}
#' @param method a character string indicating which method to use. The options for the \code{method} argument are as follows.
#' \tabular{ll}{
#' \code{"wls"} \tab an abbreviation for approximation formulae of BLUE estimators proposed by Wan et al. (2014), Luo et al. (2018) and Shi et al. (2020). Only applicable when \code{dist} is set to \code{"normal"} and \code{scenario} is set to \code{"S1"}, \code{"S2"} or \code{"S3"}. \cr
#' \code{"bala"} \tab Balakrishnan et al. (2022)'s approximation method. Only applicable when \code{dist} is set to \code{"normal"}. \cr
#' \code{"yang"} \tab Yang et al. (2022)'s approximation method. \cr
#' \code{"approx"} \tab proposed by the authors of this package to handle the Laplace and logistic distribution. \cr}
#' 
#' @returns A vector of length two, the first element is \eqn{\hat\mu} and the second one is \eqn{\hat\sigma}.
#' 
#' @details
#' \code{metablue} implements BLUE methods for estimating the mean and standard deviation from summary statistics, 
#' including methods proposed by the authors. 
#' Below are important notes regarding the function's arguments.
#' \itemize{
#'    \item \code{method = "wls"} can only be used under normality assumption and three basic scenarios. 
#'    \item \code{method = "bala"} can only be used under normality assumption. Balakrishnan et al. (2022)'s approximation method is proposed for three basic scenarios only, and the authors of this package extend their method to tertiles, quintiles and deciles. Additional message will be printed when using the method in extended scenarios. 
#'    \item \code{method = "yang"} can be used for all six scenarios under normality assumption. The authors of this package extend their method to the Laplace and logistic distribution. Additional message will be printed when using the method in extended distributions.
#'    \item \code{method = "approx"} is proposed for the Laplace and logistic distribution and three basic scenarios, inspired by Wan et al. (2014), Luo et al. (2018) and Shi et al. (2020)'s works.
#' }
#' 
#' @examples
#' set.seed(1)
#' n=100
#' r=rnorm(n, 1, 1)
#' fivenumber=fivenum(r)
#' metablue(fivenumber, n, "S3")
#' 
#' @author Sung Manhin \email{songwenxuan@ruc.edu.cn}
#' 
#' @references Luo D, Wan X, Liu J, and Tong T. (2016). Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. \emph{Statistical Methods in Medical Research}, arXiv:1505.05687.
#' @references
#' @references
#' @references
#' @references
#' 
#' @export
metablue = function(summary = NULL,
                    n = NULL,
                    scenario = "S1",
                    dist = "normal",
                    method = "wls") {
  if (is.null(summary)) {
    stop("No summary statistics are given.")
  }
  
  if (is.null(n)) {
    stop("No sample size is given.")
  }
  
  if (!(is.numeric(n) && n > 0 && n %% 1 == 0)) {
    stop("Sample size is not a positive integer.")
  }
  
  if (!scenario %in% c("S1", "S2", "S3", "tertiles", "quintiles", "deciles")) {
    stop(
      "Paramter 'scenario' must be set to 'S1', 'S2', 'S3', 'tertiles', 'quintiles' or 'deciles'. For more information, use ?metablue."
    )
  }
  
  if (length(scenario) != 1) {
    stop("Please handle one scenario at a time.")
  }
  
  if (!dist %in% c("normal", "laplace", "logistic")) {
    stop(
      "Paramter 'dist' must be set to 'normal', 'laplace' or 'logistic'. For more information, use ?metablue."
    )
  }
  
  if (length(dist) != 1) {
    stop("Please use one distributional assumption at a time.")
  }
  
  if (!method %in% c("wls", "bala", "yang", "approx")) {
    stop("Paramter 'method' must be set to 'wls', 'bala', 'yang' or 'approx'.")
  }
  
  if (length(dist) != 1) {
    stop("Please choose one method at a time.")
  }
  
  if (method == "wls") {
    if ((dist != "normal") || (!scenario %in% c("S1", "S2", "S3"))) {
      stop(
        "Wan et al. (2014), Luo et al. (2018) and Shi et al. (2020)'s methods are only applicable under normality assumption and three classic scenarios."
      )
    }
    if (scenario == "S1") {
      if (length(summary) != 3) {
        stop("Incorrect length of summary statistics in scenario 1.")
      }
      return(as.numeric(S1_wls(summary, n)))
    }
    else if (scenario == "S2") {
      if (length(summary) != 3) {
        stop("Incorrect length of summary statistics in scenario 2.")
      }
      return(as.numeric(S2_wls(summary, n)))
    }
    else if (scenario == "S3") {
      if (length(summary) != 5) {
        stop("Incorrect length of summary statistics in scenario 3.")
      }
      return(as.numeric(S3_wls(summary, n)))
    }
  }
  
  if (method == "bala") {
    if (dist != "normal") {
      stop(
        "Balakrishnan et al. (2022)'s methods are only applicable under normality assumption."
      )
    }
    if (scenario == "S1") {
      if (length(summary) != 3) {
        stop("Incorrect length of summary statistics in scenario 1.")
      }
      return(as.numeric(S1_bala(summary, n)))
    }
    else if (scenario == "S2") {
      if (length(summary) != 3) {
        stop("Incorrect length of summary statistics in scenario 2.")
      }
      return(as.numeric(S2_bala(summary, n)))
    }
    else if (scenario == "S3") {
      if (length(summary) != 5) {
        stop("Incorrect length of summary statistics in scenario 3.")
      }
      return(as.numeric(S3_bala(summary, n)))
    }
    else if (scenario == "tertiles") {
      if (length(summary) != 2) {
        stop("Incorrect length of tertiles.")
      }
      #message(
      #  "This is an extension of Balakrishnan et al. (2022)'s original methods on other summary statistics."
      #)
      return(as.numeric(tertiles_bala(summary, n)))
    }
    else if (scenario == "quintiles") {
      if (length(summary) != 4) {
        stop("Incorrect length of quintiles.")
      }
      #message(
      #  "This is an extension of Balakrishnan et al. (2022)'s original methods on other summary statistics."
      #)
      return(as.numeric(quintiles_bala(summary, n)))
    }
    else if (scenario == "deciles") {
      if (length(summary) != 9) {
        stop("Incorrect length of deciles.")
      }
      #message(
      #  "This is an extension of Balakrishnan et al. (2022)'s original methods on other summary statistics."
      #)
      
      return(as.numeric(deciles_bala(summary, n)))
    }
    
  }
  
  if (method == "yang") {
    if (dist == "normal") {
      if (scenario == "S1") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 1.")
        }
        return(as.numeric(S1_yang(summary, n)))
      }
      else if (scenario == "S2") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 2.")
        }
        return(as.numeric(S2_yang(summary, n)))
      }
      else if (scenario == "S3") {
        if (length(summary) != 5) {
          stop("Incorrect length of summary statistics in scenario 3.")
        }
        return(as.numeric(S3_yang(summary, n)))
      }
      else if (scenario == "tertiles") {
        if (length(summary) != 2) {
          stop("Incorrect length of tertiles.")
        }
        return(as.numeric(tertiles_yang(summary, n)))
      }
      else if (scenario == "quintiles") {
        if (length(summary) != 4) {
          stop("Incorrect length of quintiles.")
        }
        return(as.numeric(quintiles_yang(summary, n)))
      }
      else if (scenario == "deciles") {
        if (length(summary) != 9) {
          stop("Incorrect length of deciles.")
        }
        return(as.numeric(deciles_yang(summary, n)))
      }
    }
    else if (dist == "laplace") {
      #message(
      #  "This is an extension of Yang et al. (2022)'s original methods on the Laplace distribution."
      #)
      if (scenario == "S1") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 1.")
        }
        return(as.numeric(S1_yang_lap(summary, n)))
      }
      else if (scenario == "S2") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 2.")
        }
        return(as.numeric(S2_yang_lap(summary, n)))
      }
      else if (scenario == "S3") {
        if (length(summary) != 5) {
          stop("Incorrect length of summary statistics in scenario 3.")
        }
        return(as.numeric(S3_yang_lap(summary, n)))
      }
      else if (scenario == "tertiles") {
        if (length(summary) != 2) {
          stop("Incorrect length of tertiles.")
        }
        return(as.numeric(tertiles_yang_lap(summary, n)))
      }
      else if (scenario == "quintiles") {
        if (length(summary) != 4) {
          stop("Incorrect length of quintiles.")
        }
        return(as.numeric(quintiles_yang_lap(summary, n)))
      }
      else if (scenario == "deciles") {
        if (length(summary) != 9) {
          stop("Incorrect length of deciles.")
        }
        return(as.numeric(deciles_yang_lap(summary, n)))
      }
    }
    else if (dist == "logistic") {
      #message(
      #  "This is an extension of Yang et al. (2022)'s original methods on the logistic distribution."
      #)
      if (scenario == "S1") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 1.")
        }
        return(as.numeric(S1_yang_logit(summary, n)))
      }
      else if (scenario == "S2") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 2.")
        }
        return(as.numeric(S2_yang_logit(summary, n)))
      }
      else if (scenario == "S3") {
        if (length(summary) != 5) {
          stop("Incorrect length of summary statistics in scenario 3.")
        }
        return(as.numeric(S3_yang_logit(summary, n)))
      }
      else if (scenario == "tertiles") {
        if (length(summary) != 2) {
          stop("Incorrect length of tertiles.")
        }
        return(as.numeric(tertiles_yang_logit(summary, n)))
      }
      else if (scenario == "quintiles") {
        if (length(summary) != 4) {
          stop("Incorrect length of quintiles.")
        }
        return(as.numeric(quintiles_yang_logit(summary, n)))
      }
      else if (scenario == "deciles") {
        if (length(summary) != 9) {
          stop("Incorrect length of deciles.")
        }
        return(as.numeric(deciles_yang_logit(summary, n)))
      }
    }
  }
  
  if (method == "approx") {
    if (!dist %in% c("laplace", "logistic")) {
      stop("Proposed for the Laplace distribution and logistic distribution only.")
    }
    if (!scenario %in% c("S1", "S2", "S3")) {
      stop("Proposed for three classic scenarios only.")
    }
    if (dist == "laplace") {
      if (scenario == "S1") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 1.")
        }
        return(as.numeric(S1_lap(summary, n)))
      }
      else if (scenario == "S2") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 2.")
        }
        return(as.numeric(S2_lap(summary, n)))
      }
      else if (scenario == "S3") {
        if (length(summary) != 5) {
          stop("Incorrect length of summary statistics in scenario 3.")
        }
        return(as.numeric(S3_lap(summary, n)))
      }
    }
    else if (dist == "logistic") {
      if (scenario == "S1") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 1.")
        }
        return(as.numeric(S1_logit(summary, n)))
      }
      else if (scenario == "S2") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 2.")
        }
        return(as.numeric(S2_logit(summary, n)))
      }
      else if (scenario == "S3") {
        if (length(summary) != 5) {
          stop("Incorrect length of summary statistics in scenario 3.")
        }
        return(as.numeric(S3_logit(summary, n)))
      }
    }
    
  }
  
  
  
  
}
