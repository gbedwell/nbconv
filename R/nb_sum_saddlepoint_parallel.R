#' Parallelized implementation of the saddlepoint approximation for the sum of arbitrary NB random variables
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
#'@param n.cores The number of CPU cores to use for evaluation.
#'
#'@export
#'
nb_sum_saddlepoint_parallel <- function(mus, phis, ps, counts.start, counts.end,
                                        normalize = TRUE, log.scale = TRUE, n.cores){
  require(parallel)
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

  if (missing(n.cores)){
    n.cores <- parallel::detectCores()/2
  }

  count.vec <- counts.start:counts.end
  count.list <- split(count.vec, ceiling(seq_along(count.vec)/floor(sqrt(n.counts))))


  pmf.list <- mclapply(X = count.list,
                       FUN = function(x) {
                         nb_saddlepoint(mus = mus,
                                        phis = phis,
                                        counts.start = min(x),
                                        counts.end = max(x),
                                        normalize = FALSE,
                                        log.scale = log.scale) },
                       mc.cores = n.cores)

  saddlepoint.pmf <- Reduce(c, pmf.list)

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
