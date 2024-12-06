#' @title BLUE Estimators for Mean and Standard Deviation
#' @description A comprehensive function for estimating the mean and standard deviation from summary statistics, incorporating both published and novel BLUE methods.
#'
#' @param summary a vector of ordered summary statistics in ascending order.
#' @param n the sample size.
#' @param scenario a character string indicating the scenario in which summary statistics is reported. The options for the \code{scenario} argument are as follows.
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
#'    \item \code{method = "approx"} is proposed for the Laplace and logistic distribution and three basic scenarios, inspired by Wan et al. (2014), Luo et al. (2018) and Shi et al. (2020)'s works, check \url{https://github.com/sungmanhin/metaBLUEX} for approximation details.
#' }
#'
#' @examples
#' ### normal in S3 ###
#' set.seed(1)
#' n = 100
#' # mu = 1, sigma = 1
#' r = rnorm(n, 1, 1)
#' fivenumber = fivenum(r)
#'
#' # Wan, Luo and Shi's method
#' metablue(fivenumber, n, "S3", "normal")
#' # 1.0981794 0.9141252
#'
#' # Yang's method
#' metablue(fivenumber, n, "S3", "normal", "yang")
#' # 1.0857491 0.9065274
#'
#' # Balakrishnan's method
#' metablue(fivenumber, n, "S3", "normal", "bala")
#' # 1.092199 0.933297
#'
#' ### laplace in S3 ###
#' library(VGAM)
#' set.seed(1)
#' n = 100
#' # mu = 1, sigma = 2
#' # Laplace dist. with parameters (0, 1/sqrt(2)) has mean 0 and sigma 1
#' l = 1 + 2 * rlaplace(n, 0, 1/sqrt(2))
#' fivenumber_l = fivenum(l)
#'
#' # Yang's method
#' metablue(fivenumber_l, n, "S3", "laplace", "yang")
#' # 0.9505243 1.7241492
#'
#' # Our approximation method
#' metablue(fivenumber_l, n, "S3", "laplace", "approx")
#' #  0.9887838 1.7481349
#'
#' @author Sung Manhin \email{songwenxuan@ruc.edu.cn}
#'
#' @references Wan, X., Wang, W., Liu, J., & Tong, T. (2014). Estimating the sample mean and standard deviation from the sample size, median, range and/or interquartile range. \emph{BMC medical research methodology, 14}, 1-13.
#' @references Luo, D., Wan, X., Liu, J., & Tong, T. (2018). Optimally estimating the sample mean from the sample size, median, mid-range, and/or mid-quartile range. \emph{Statistical methods in medical research, 27}(6), 1785-1805.
#' @references Shi, J., Luo, D., Weng, H., Zeng, X. T., Lin, L., Chu, H., & Tong, T. (2020). Optimally estimating the sample standard deviation from the five‐number summary. \emph{Research synthesis methods, 11}(5), 641-654.
#' @references Balakrishnan, N., Rychtář, J., Taylor, D., & Walter, S. D. (2022). Unified approach to optimal estimation of mean and standard deviation from sample summaries. \emph{Statistical Methods in Medical Research, 31}(11), 2087-2103.
#' @references Yang, X., Hutson, A. D., & Wang, D. (2022). A generalized BLUE approach for combining location and scale information in a meta-analysis. \emph{Journal of Applied Statistics, 49}(15), 3846-3867.
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
  
  if (!is.numeric(summary)) {
    stop("Summary statistics are not of numeric type.")
  }
  
  if (is.null(n)) {
    stop("No sample size is given.")
  }
  
  if (!(is.numeric(n) && n > 0 && n %% 1 == 0)) {
    stop("Sample size should be a positive integer.")
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
    stop(
      "Paramter 'method' must be set to 'wls', 'bala', 'yang' or 'approx'. For more information, use ?metablue."
    )
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
      return(as.numeric(.S1_wls(summary, n)))
    }
    else if (scenario == "S2") {
      if (length(summary) != 3) {
        stop("Incorrect length of summary statistics in scenario 2.")
      }
      return(as.numeric(.S2_wls(summary, n)))
    }
    else if (scenario == "S3") {
      if (length(summary) != 5) {
        stop("Incorrect length of summary statistics in scenario 3.")
      }
      return(as.numeric(.S3_wls(summary, n)))
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
      return(as.numeric(.S1_bala(summary, n)))
    }
    else if (scenario == "S2") {
      if (length(summary) != 3) {
        stop("Incorrect length of summary statistics in scenario 2.")
      }
      return(as.numeric(.S2_bala(summary, n)))
    }
    else if (scenario == "S3") {
      if (length(summary) != 5) {
        stop("Incorrect length of summary statistics in scenario 3.")
      }
      return(as.numeric(.S3_bala(summary, n)))
    }
    else if (scenario == "tertiles") {
      if (length(summary) != 2) {
        stop("Incorrect length of tertiles.")
      }
      #message(
      #  "This is an extension of Balakrishnan et al. (2022)'s original methods on other summary statistics."
      #)
      return(as.numeric(.tertiles_bala(summary, n)))
    }
    else if (scenario == "quintiles") {
      if (length(summary) != 4) {
        stop("Incorrect length of quintiles.")
      }
      #message(
      #  "This is an extension of Balakrishnan et al. (2022)'s original methods on other summary statistics."
      #)
      return(as.numeric(.quintiles_bala(summary, n)))
    }
    else if (scenario == "deciles") {
      if (length(summary) != 9) {
        stop("Incorrect length of deciles.")
      }
      #message(
      #  "This is an extension of Balakrishnan et al. (2022)'s original methods on other summary statistics."
      #)
      
      return(as.numeric(.deciles_bala(summary, n)))
    }
    
  }
  
  if (method == "yang") {
    if (dist == "normal") {
      if (scenario == "S1") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 1.")
        }
        return(as.numeric(.S1_yang(summary, n)))
      }
      else if (scenario == "S2") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 2.")
        }
        return(as.numeric(.S2_yang(summary, n)))
      }
      else if (scenario == "S3") {
        if (length(summary) != 5) {
          stop("Incorrect length of summary statistics in scenario 3.")
        }
        return(as.numeric(.S3_yang(summary, n)))
      }
      else if (scenario == "tertiles") {
        if (length(summary) != 2) {
          stop("Incorrect length of tertiles.")
        }
        return(as.numeric(.tertiles_yang(summary, n)))
      }
      else if (scenario == "quintiles") {
        if (length(summary) != 4) {
          stop("Incorrect length of quintiles.")
        }
        return(as.numeric(.quintiles_yang(summary, n)))
      }
      else if (scenario == "deciles") {
        if (length(summary) != 9) {
          stop("Incorrect length of deciles.")
        }
        return(as.numeric(.deciles_yang(summary, n)))
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
        return(as.numeric(.S1_yang_lap(summary, n)))
      }
      else if (scenario == "S2") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 2.")
        }
        return(as.numeric(.S2_yang_lap(summary, n)))
      }
      else if (scenario == "S3") {
        if (length(summary) != 5) {
          stop("Incorrect length of summary statistics in scenario 3.")
        }
        return(as.numeric(.S3_yang_lap(summary, n)))
      }
      else if (scenario == "tertiles") {
        if (length(summary) != 2) {
          stop("Incorrect length of tertiles.")
        }
        return(as.numeric(.tertiles_yang_lap(summary, n)))
      }
      else if (scenario == "quintiles") {
        if (length(summary) != 4) {
          stop("Incorrect length of quintiles.")
        }
        return(as.numeric(.quintiles_yang_lap(summary, n)))
      }
      else if (scenario == "deciles") {
        if (length(summary) != 9) {
          stop("Incorrect length of deciles.")
        }
        return(as.numeric(.deciles_yang_lap(summary, n)))
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
        return(as.numeric(.S1_yang_logit(summary, n)))
      }
      else if (scenario == "S2") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 2.")
        }
        return(as.numeric(.S2_yang_logit(summary, n)))
      }
      else if (scenario == "S3") {
        if (length(summary) != 5) {
          stop("Incorrect length of summary statistics in scenario 3.")
        }
        return(as.numeric(.S3_yang_logit(summary, n)))
      }
      else if (scenario == "tertiles") {
        if (length(summary) != 2) {
          stop("Incorrect length of tertiles.")
        }
        return(as.numeric(.tertiles_yang_logit(summary, n)))
      }
      else if (scenario == "quintiles") {
        if (length(summary) != 4) {
          stop("Incorrect length of quintiles.")
        }
        return(as.numeric(.quintiles_yang_logit(summary, n)))
      }
      else if (scenario == "deciles") {
        if (length(summary) != 9) {
          stop("Incorrect length of deciles.")
        }
        return(as.numeric(.deciles_yang_logit(summary, n)))
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
        return(as.numeric(.S1_lap(summary, n)))
      }
      else if (scenario == "S2") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 2.")
        }
        return(as.numeric(.S2_lap(summary, n)))
      }
      else if (scenario == "S3") {
        if (length(summary) != 5) {
          stop("Incorrect length of summary statistics in scenario 3.")
        }
        return(as.numeric(.S3_lap(summary, n)))
      }
    }
    else if (dist == "logistic") {
      if (scenario == "S1") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 1.")
        }
        return(as.numeric(.S1_logit(summary, n)))
      }
      else if (scenario == "S2") {
        if (length(summary) != 3) {
          stop("Incorrect length of summary statistics in scenario 2.")
        }
        return(as.numeric(.S2_logit(summary, n)))
      }
      else if (scenario == "S3") {
        if (length(summary) != 5) {
          stop("Incorrect length of summary statistics in scenario 3.")
        }
        return(as.numeric(.S3_logit(summary, n)))
      }
    }
    
  }
  
  
}


#' @title Magic Wrapper for 'BLUE Estimators for Mean and Standard Deviation'
#' @description A magic wrapper for \code{metablue} to provide a colorful and detailed report for applicable distributions and methods in estimation, once and for all.
#'
#' @param summary a vector of ordered summary statistics in ascending order.
#' @param n the sample size.
#' @param scenario a character string indicating the scenario in which summary statistics is reported. The options for the \code{scenario} argument are as follows.
#' \tabular{ll}{
#' \code{"S1"} \tab for the median, minimum and maximum values. \cr
#' \code{"S2"} \tab for the median, first and third quartiles. \cr
#' \code{"S3"} \tab for the median, first and third quartiles, and minimum and maximum values. \cr
#' \code{"tertiles"} \tab for tertiles. \cr
#' \code{"quintiles"} \tab for quintiles. \cr
#' \code{"deciles"} \tab for deciles. \cr}
#' @param object an object of class \code{metabluex}.
#' @param ... additional arguments (not used).
#'
#' @returns A list with class "metabluex" containing the following components:
#' \tabular{ll}{
#' \code{summary} \tab summary statistics. \cr
#' \code{n} \tab the sample size. \cr
#' \code{scenario} \tab the scenario in which summary statistics is reported. \cr
#' \code{result} \tab a list containing results (\eqn{\hat\mu} and \eqn{\hat\sigma} as a vector) in all distributions and methods applicable for the selected scenario. \cr }
#'
#' @details
#' \code{metabluex} is a magic wrapper for \code{metablue}, which outputs all possible estimations in all distributions and methods applicable for the selected scenario.
#'
#' This function is highly useful in scenarios where you have determined summary statistics and want to explore different methods and distributional assumptions for estimating the mean and standard deviation.
#' You may acquire a tidy report to elegantly review and compare all possible BLUE methods applicable under the current scenario.
#'
#' Consistent with the applicability of \code{"dist"} and \code{"method"} in \code{metablue},
#' you only need to specify \code{"scenario"} to get the estimation results.
#' Generally, you don’t need to check about the details below:
#' \itemize{
#'    \item if \code{scenario} is \code{"S1"}, \code{"S2"} or \code{"S3"}, report contains estimates in \code{"normal"} for \code{"wls"}, \code{"bala"}, \code{"yang"}, in \code{"laplace"} for \code{"yang"}, \code{"approx"} and in \code{"logistic"} for \code{"yang"}, \code{"approx"}.
#'    \item if \code{scenario} is \code{"tertiles"}, \code{"quintiles"} or \code{"deciles"}, report contains estimates in \code{"normal"} for \code{"bala"}, \code{"yang"}, in \code{"laplace"} for \code{"yang"} and in \code{"logistic"} for \code{"yang"}.
#' }
#' \code{crayon} package is used to report in a colorful way.
#'
#' @examples
#' set.seed(1)
#' n = 100
#' # mu = 1, sigma = 1
#' r = rnorm(n, 1, 1)
#'
#' ### S3 ###
#' fivenumber = fivenum(r)
#' metabluex(fivenumber, n, "S3")
#'
#' ### tertiles ###
#' tertiles = quantile(r, c(1/3, 2/3), names = FALSE)
#' metabluex(tertiles, n, "tertiles")
#'
#' @author Sung Manhin \email{songwenxuan@ruc.edu.cn}
#'
#' @export
metabluex = function(summary = NULL,
                     n = NULL,
                     scenario = "S1") {
  if (is.null(summary)) {
    stop("No summary statistics are given.")
  }
  
  if (!is.numeric(summary)) {
    stop("Summary statistics are not of numeric type.")
  }
  
  if (is.null(n)) {
    stop("No sample size is given.")
  }
  
  if (!(is.numeric(n) && n > 0 && n %% 1 == 0)) {
    stop("Sample size is not a positive integer.")
  }
  
  if (!scenario %in% c("S1", "S2", "S3", "tertiles", "quintiles", "deciles")) {
    stop(
      "Paramter 'scenario' must be set to 'S1', 'S2', 'S3', 'tertiles', 'quintiles' or 'deciles'. For more information, use ?metabluex."
    )
  }
  
  if (length(scenario) != 1) {
    stop("Please handle one scenario at a time.")
  }
  
  if (scenario %in% c("tertiles", "quintiles", "deciles")) {
    param1 = list(
      summary = summary,
      n = n,
      scenario = scenario,
      dist = "normal",
      method = "bala"
    )
    param2 = list(
      summary = summary,
      n = n,
      scenario = scenario,
      dist = "normal",
      method = "yang"
    )
    param3 = list(
      summary = summary,
      n = n,
      scenario = scenario,
      dist = "laplace",
      method = "yang"
    )
    param4 = list(
      summary = summary,
      n = n,
      scenario = scenario,
      dist = "logistic",
      method = "yang"
    )
    results = list(
      summary = summary,
      n = n,
      scenario = scenario,
      result = list(
        normal_bala = do.call(metablue, param1),
        normal_yang = do.call(metablue, param2),
        laplace_yang = do.call(metablue, param3),
        logistic_yang = do.call(metablue, param4)
      )
    )
  }
  else{
    param5 = list(
      summary = summary,
      n = n,
      scenario = scenario,
      dist = "normal",
      method = "wls"
    )
    param6 = list(
      summary = summary,
      n = n,
      scenario = scenario,
      dist = "normal",
      method = "bala"
    )
    param7 = list(
      summary = summary,
      n = n,
      scenario = scenario,
      dist = "normal",
      method = "yang"
    )
    param8 = list(
      summary = summary,
      n = n,
      scenario = scenario,
      dist = "laplace",
      method = "yang"
    )
    param9 = list(
      summary = summary,
      n = n,
      scenario = scenario,
      dist = "laplace",
      method = "approx"
    )
    param10 = list(
      summary = summary,
      n = n,
      scenario = scenario,
      dist = "logistic",
      method = "yang"
    )
    param11 = list(
      summary = summary,
      n = n,
      scenario = scenario,
      dist = "logistic",
      method = "approx"
    )
    results = list(
      summary = summary,
      n = n,
      scenario = scenario,
      result = list(
        normal_wls = do.call(metablue, param5),
        normal_bala = do.call(metablue, param6),
        normal_yang = do.call(metablue, param7),
        laplace_yang = do.call(metablue, param8),
        laplace_approx = do.call(metablue, param9),
        logistic_yang = do.call(metablue, param10),
        logistic_approx = do.call(metablue, param11)
      )
    )
    
  }
  
  class(results) = "metabluex"
  return(results)
  
}

#' @rdname metabluex
#' @export
print.metabluex = function(object, ...) {
  summary = object$summary
  n = object$n
  scenario = object$scenario
  cat(crayon::bgBlue(
    paste0(
      strrep(" ", 8),
      "MetaBLUEX Summary2MeanSD Report",
      " (",
      scenario,
      ")",
      strrep(" ", 8)
    )
  ), "\n")
  cat(crayon::green(paste0(
    " summary = ", paste(round(summary, 3), collapse = ", ")
  )), "\n")
  cat(crayon::green(paste0("                   ", "sample size = ", n, "                 ")), "\n")
  for (name in names(object$result)) {
    name1 = strsplit(name, "_")[[1]][1]
    name2 = strsplit(name, "_")[[1]][2]
    convert_name = c(
      "normal" = " Normal ",
      "laplace" = "Laplace ",
      "logistic" = "Logistic",
      "wls" = " Wan, Luo & Shi ",
      "yang" = "      Yang      ",
      "bala" = "  Balakrishnan  ",
      "approx" = "  Our approx.   "
    )
    if (name2 == "wls") {
      cat(
        crayon::bgBlue(paste0(" ", convert_name[name1], " ")),
        crayon::bgGreen(paste0(" ", convert_name[name2], " ")),
        crayon::cyan(paste0(
          "mean = ",
          round(object$result[[name]][1], 3),
          " ",
          "sd = ",
          round(object$result[[name]][2], 3),
          " "
        )),
        "\n"
      )
    } else if (name2 == "yang") {
      cat(
        crayon::bgBlue(paste0(" ", convert_name[name1], " ")),
        crayon::bgYellow(paste0(" ", convert_name[name2], " ")),
        crayon::cyan(paste0(
          "mean = ",
          round(object$result[[name]][1], 3),
          " ",
          "sd = ",
          round(object$result[[name]][2], 3),
          " "
        )),
        "\n"
      )
    } else if (name2 == "bala") {
      cat(
        crayon::bgBlue(paste0(" ", convert_name[name1], " ")),
        crayon::bgCyan(paste0(" ", convert_name[name2], " ")),
        crayon::cyan(paste0(
          "mean = ",
          round(object$result[[name]][1], 3),
          " ",
          "sd = ",
          round(object$result[[name]][2], 3),
          " "
        )),
        "\n"
      )
    } else if (name2 == "approx") {
      cat(
        crayon::bgBlue(paste0(" ", convert_name[name1], " ")),
        crayon::bgMagenta(paste0(" ", convert_name[name2], " ")),
        crayon::cyan(paste0(
          "mean = ",
          round(object$result[[name]][1], 3),
          " ",
          "sd = ",
          round(object$result[[name]][2], 3),
          " "
        )),
        "\n"
      )
    }
  }
  invisible(object)
}

