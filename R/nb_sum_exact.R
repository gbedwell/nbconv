#' Implements Furman's exact PMF for the evaluation of the sum of arbitrary NB random variables
#'
#'@param mus Vector of individual mean values
#'@param phis Vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param ps Vector of individual probabilities.
#'@param n.terms The number of terms to include in the series for evaluating the PMF at a given number of counts. Defaults to 1000.
#'@param counts.start The smallest number of counts at which the PMF is to be evaluated. Defaults to 0.
#'@param counts.end The largest number of counts at which the PMF is to be evaluated.
#'@param n.cores The number of CPU cores to use in the evaluation. Allows parallelization.
#'
#'@import parallel
#'
#'@export
#'
nb_sum_exact <- function(mus, phis, ps, n.terms = 1000, counts.start = 0, counts.end, n.cores = 1){
  # Adapted from https://github.com/slundberg/NBConvolution.jl
  # Implements the PMF described in https://ssrn.com/abstract=1650365

  if (!missing(ps)){
    if (!missing(mus)){
      stop("'mus' and 'ps' both specified",
           call. = FALSE)
    }
  }

  if (missing(ps)){
    if (!missing(mus)){
      ps <- phis/(phis + mus)
    }
    else{
      stop("One of 'mus' or 'ps' must be specified.",
           call. = FALSE)
    }
  }

  if (length(ps) != length(phis)){
    stop("'ps' and 'phis' must have the same length.", call. = FALSE)
  }

  pmax <- max(ps)
  qmax <- 1 - pmax
  phisum <- sum(phis)

  R <- 1

  for (j in 1:length(ps)){
    R <- R * (((1 - ps[j]) * pmax)/(qmax * ps[j]))^(-phis[j])
  }

  xi <- rep(x = 0, times = n.terms)
  xtmp <- rep(x = 0, times = length(ps))
  for (i in 1:n.terms){
    for (j in 1:length(ps)){
      xtmp[j] <- (phis[j] * (1 - qmax * ps[j]/((1 - ps[j]) * pmax))^i) / i
    }
    xi[i] <- sum(xtmp)
  }

  delta <- rep(x = 0, times = n.terms)
  dtmp <- rep(x = 0, times = n.terms)
  delta[1] <- 1

  for (k in 1:(n.terms - 1)){
    for (i in 1:k) {
      previndex <- k + 1 - i
      dtmp[i] <- i * xi[i] * delta[previndex]
      }
    delta[k + 1] <- sum(dtmp) / k
    }

  mass_calc <- function(x){
    total <- 0
    for (k in 0:(n.terms - 1)){
      v <- delta[k + 1]*exp(lgamma(phisum + x + k) - lgamma(phisum + k) - lfactorial(x) + (phisum + k) * log(pmax) + x * log(1 - pmax))
      total <- total + v
    }
    if (any(is.na(v) | any(is.infinite(v)))){
      stop("Values out of bounds. PMF cannot be evaluated.", call. = FALSE)
    }
    if ( (v[length(v)] > 1e-10) ||  (v[length(v) - 1] < v[length(v)] && v[length(v)] > 1e-50) ){
      stop("Insufficient sum expansion. Use more terms.", call. = FALSE)
    }
    masses <- total*R
    return(masses)
  }

  if (n.cores == 1){
    v <- counts.start:counts.end
    pmf <- mass_calc(x = v)
  }
  else {
    require(parallel)

    v <- counts.start:counts.end
    v.list <- split(v, ceiling((seq_along(v))/1000))

    pmf.list <- mclapply(X = v.list,
                         FUN = function(y) {
                           v.split <- min(y):max(y)
                           pmf <- mass_calc(x = v.split)
                           return(pmf) },
                         mc.cores = n.cores)

    pmf <- Reduce(c, pmf.list)
  }

  if (is.numeric(pmf)){
    if (sum(pmf) < 0.9999){
      warning("The sum of the evaluated distribution is less than 0.9999. Consider expanding the range.",
              call. = FALSE)
    }
    return(pmf)
  }
  else{
    error <- sub("Error : *", "", pmf[1])
    stop(paste0(error, "\n"), call. = FALSE)
  }
}

