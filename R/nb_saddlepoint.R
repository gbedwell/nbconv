#' Implements the saddlepoint approximation for the sum of arbitrary NB random variables
#'
#' Adapted from https://stats.stackexchange.com/questions/72479/generic-sum-of-gamma-random-variables/137318#137318
#' and https://www.martinmodrak.cz/2019/06/20/approximate-densities-for-sums-of-variables-negative-binomials-and-saddlepoint/
#'
#'@param mus Vector of individual mean values
#'@param phis Vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param ps Vector of individual probabilities.
#'@param counts.start The smallest number of counts at which the PMF is to be evaluated. Defaults to 0.
#'@param counts.end The largest number of counts at which the PMF is to be evaluated.
#'@param log.scale Boolean. If TRUE, saddlepoint approximation is done on the log scale.
#'@param normalize Boolean. If TRUE, the PMF is re-normalized to sum to 1.
#'
#'@export
#'
nb_saddlepoint <- function(mus, phis, ps, counts.start = 0, counts.end, log.scale = TRUE, normalize = TRUE){

  require(matrixStats)

  if (!missing(mus)){
    if (!missing(ps)){
      stop("'mus' and 'ps' both specified",
           call. = FALSE)
    }
  }

  if (missing(mus)){
    if (!missing(ps)){
      mus <- phis*(1 - ps)/ps
    }
    else{
      stop("One of 'mus' or 'ps' must be specified.",
           call. = FALSE)
    }
  }

  if (length(mus) != length(phis)){
    stop("'mus' and 'phis' must have the same length.", call. = FALSE)
  }

  if(isTRUE(log.scale)){
    K  <-  function(s) { sum(phis*(log(phis) - log(phis + mus*(1 - exp(s))))) }
    Kd <-  function(s) { logSumExp(log(phis) + log(mus) + s - log(phis + mus - mus * exp(s))) }
    Kdd <- function(s) { logSumExp(log(phis) + log(mus) + log(phis + mus) + s - 2 * log(phis + mus - mus * exp(s))) }
    pmf_eq <- function(s, x) { -0.5 * (log(2*pi) + Kdd(s)) + K(s) - s * x }
  }
  else{
    K  <-  function(s) { sum(phis*(log(phis) - log(phis + mus*(1 - exp(s))))) }
    Kd <-  function(s) { sum((phis*mus*exp(s))/(phis + mus*(1 - exp(s)))) }
    Kdd <- function(s) { sum((phis*mus*(phis + mus)*exp(s))/(phis + mus*(1 - exp(s)))^2) }
    pmf_eq <- function(s, x) { (1/sqrt(2 * pi * Kdd(s)))*exp(K(s) - s * x) }
  }

  if (counts.start == 0){
    x.start <- 1
    pmf0 <- prod(pnbinom(0, size = phis, mu = mus))
  }
  else{
    x.start <- counts.start
  }

  pmf <- sapply(X = x.start:counts.end,
                FUN = function(x) {
                  # suppressWarnings() used here to suppress "NaNs produced" warning.
                  # warnings() returns 'In log(phis + mus - mus * exp(s)) : NaNs produced'
                  # uniroot is trying to take the log of a negative number when mus * exp(s) > phis + mus
                  # this warning doesn't appear to affect uniroot's output
                  if(isTRUE(log.scale)){
                    s <- suppressWarnings(uniroot(function(s) Kd(s) - log(x),
                                                  lower = -1e2,
                                                  upper = 0,
                                                  extendInt = "yes",
                                                  tol = .Machine$double.eps)$root)
                  }
                  else{
                    s <- suppressWarnings(uniroot(function(s) Kd(s) - x,
                                                  lower = -1e2,
                                                  upper = 0,
                                                  extendInt = "yes",
                                                  tol = .Machine$double.eps)$root)
                  }
                  pmf <- pmf_eq(s, x)
                  return(pmf)
                })

  if (exists("pmf0")){
    if(isTRUE(log.scale)){
      saddlepoint.pmf <- c(pmf0, exp(pmf))
    }
    else {
      saddlepoint.pmf <- c(pmf0, pmf)
    }
  }
  else{
    if(isTRUE(log.scale)){
      saddlepoint.pmf <- exp(pmf)
    }
    else {
      saddlepoint.pmf <- pmf
    }
  }

  if ((saddlepoint.pmf[1] > 1e-5 & counts.start != 0) || saddlepoint.pmf[length(saddlepoint.pmf)] > 1e-5){
    near.zero <- FALSE
  }
  else{
    near.zero <- TRUE
  }

  if (!isTRUE(near.zero)){
    warning("The density values at one or both of the ends of the given range are > 1e-5. Consider increasing the evaluated range.",
            call. = FALSE)
  }

  if (isTRUE(normalize)){
    if (isTRUE(near.zero)){
      saddlepoint.pmf <- saddlepoint.pmf/sum(saddlepoint.pmf)
    }
    else{
      stop("The ends of the distribution are too far from zero to normalize. Increase the evaluated range to normalize.",
           call. = FALSE)
    }
  }

  return(saddlepoint.pmf)
}
