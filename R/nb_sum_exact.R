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
  qmax <- 1-pmax
  phisum <- sum(phis)

  R <- 1

  for (j in 1:length(ps)){
    R <- R*(((1-ps[j])*pmax)/(qmax*ps[j]))^(-phis[j])
  }

  xi <- rep(0, n.terms)
  xtmp <- rep(0, length(ps))
  for (i in 1:n.terms){
    for (j in 1:length(ps)){
      xtmp[j] <- (phis[j]*(1-qmax*ps[j]/((1-ps[j])*pmax))^i)/i
    }
    xi[i] <- sum(xtmp)
  }

  delta <- rep(0, n.terms)
  dtmp <- rep(0, n.terms)
  delta[1] <- 1

  for (k in 1:(n.terms-1)){
    for (i in 1:k) {
      previndex <- k+1-i
      dtmp[i] <- i*xi[i]*delta[previndex]
    }
    delta[k+1] <- sum(dtmp)/k
  }

  mass_calc <- function(x = s){
    total <- 0
    lastv <- 0

    for (k in 0:(n.terms-1)){
      v <- delta[k+1]*exp(lgamma(phisum + x + k) - lgamma(phisum + k) - lfactorial(x) + (phisum + k)*log(pmax) + x*log(1-pmax))
      total <- total + v

      if (k == (n.terms-1)){
        if (any(is.na(v)) || any(is.infinite(v))){
          stop("Values out of bounds. PMF can not be computed.", call. = FALSE)
        }
        if (any(v > 1e-10) || any(lastv < v & v > 1e-50)){
          stop("Insufficient sum expansion. Use more terms.", call. = FALSE)
        }
      }
      lastv <- v
    }
    masses <- total*R
    return(masses)
  }

  if (n.cores == 1){
    s <- counts.start:counts.end
    pmf <- mass_calc(x = s)
  }
  else {
    require(parallel)

    count.vec <- counts.start:counts.end
    count.list <- split(count.vec, ceiling(seq_along(count.vec)/floor(sqrt(n.counts))))

    pmf.list <- mclapply(X = count.list,
                         FUN = function(y) {
                           s <- min(y):max(y)
                           pmf <- mass_calc(x = s)
                           return(pmf) },
                         mc.cores = n.cores)

    pmf <- Reduce(c, pmf.list)

  }

  if (pmf[1] > 1e-5 && counts.start != 0 || pmf[length(pmf)] > 1e-5){
    warning("The density values at one or both of the ends of the given range are > 1e-5. Consider increasing the evaluated range.",
            call. = FALSE)
  }

  return(pmf)
}

