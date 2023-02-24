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
#'@param normalize Boolean. If TRUE, the PMF is re-normalized to sum to 1.
#'
#'@import matrixStats
#'@import parallel
#'
#'@export
#'
nb_sum_saddlepoint <- function(mus, phis, ps, counts.start = 0, counts.end, normalize = TRUE, n.cores = 1){

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

  saddlepoint_calc <- function(mus, phis, counts.start, counts.end){

    K  <-  function(s) { sum(phis*(log(phis) - log(phis + mus*(1 - exp(s))))) }
    Kd <-  function(s) { logSumExp(log(phis) + log(mus) + s - log(phis + mus - mus * exp(s))) }
    Kdd <- function(s) { logSumExp(log(phis) + log(mus) + log(phis + mus) + s - 2 * log(phis + mus - mus * exp(s))) }
    pmf_eq <- function(s, x) { -0.5 * (log(2 * pi) + Kdd(s)) + K(s) - s * x }

    if (counts.start == 0){
      x.start <- 1
      pmf0 <- prod(pnbinom(0, size = phis, mu = mus))
    }
    else{
      x.start <- counts.start
    }

    v <- x.start:counts.end

    pmf <- sapply(X = v,
                  FUN = function(x) {
                    # suppressWarnings() used here to suppress "NaNs produced" warning.
                    # warnings() returns 'In log(phis + mus - mus * exp(s)) : NaNs produced'
                    # uniroot is trying to take the log of a negative number when mus * exp(s) > phis + mus
                    # this warning doesn't appear to affect uniroot's output
                    s <- suppressWarnings(uniroot(function(s) Kd(s) - log(x),
                                                  lower = -1e2,
                                                  upper = 0,
                                                  extendInt = "yes",
                                                  tol = .Machine$double.eps)$root)
                    pmf <- pmf_eq(s, x)
                    return(pmf)
                  }
    )

    if (exists("pmf0")){
      pmf <- c(pmf0, exp(pmf))
    }
    else{
      pmf <- exp(pmf)
    }
    return(pmf)
  }

  if (n.cores == 1){
    saddlepoint.pmf <- saddlepoint_calc(mus = mus,
                                        phis = phis,
                                        counts.start = counts.start,
                                        counts.end = counts.end)
  }
  else{
    v <- counts.start:counts.end
    v.list <- split(v, ceiling((seq_along(v))/1000))

    pmf.list <- mclapply(X = v.list,
                         FUN = function(y) {
                           split.start <- min(y)
                           split.end <- max(y)
                           pmf <- saddlepoint_calc(mus = mus,
                                                   phis = phis,
                                                   counts.start = split.start,
                                                   counts.end = split.end)
                           return(pmf) },
                         mc.cores = n.cores)

    saddlepoint.pmf <- Reduce(c, pmf.list)
  }

  left <- max(saddlepoint.pmf)/saddlepoint.pmf[1]
  right <- max(saddlepoint.pmf)/saddlepoint.pmf[length(saddlepoint.pmf)]

  if ((left < 1e10 && counts.start != 0) || right < 1e10){
    small <- FALSE
  }
  else{
    small <- TRUE
  }

  if (isFALSE(small)){
    warning("The density values at one or both of the ends of the given range are insufficiently small. Consider expanding the evaluated range.",
            call. = FALSE)
  }

  if (isTRUE(normalize)){
    if (isTRUE(small)){
      saddlepoint.pmf <- saddlepoint.pmf/sum(saddlepoint.pmf)
    }
    else{
      stop("Refusing to normalize. See generated warning and expand the evaluated range.",
           call. = FALSE)
    }
  }

  return(saddlepoint.pmf)
}
