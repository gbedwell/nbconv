#' Implements the method of moments approximation for the sum of arbitrary NB random variables
#'
#'@param mus Vector of individual mean values
#'@param phis Vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param ps Vector of individual probabilities.
#'@param counts.start The smallest number of counts at which the PMF is to be evaluated. Defaults to 0.
#'@param counts.end The largest number of counts at which the PMF is to be evaluated.
#'
#'@export
#'
nb_sum_moments <- function(mus, phis, ps, counts.start = 0, counts.end){

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

  s <- counts.start:counts.end

  mu.moment <- sum(mus)
  phi.moment <- sum(mus)^2/sum(mus^2/phis)

  moments.pmf <- dnbinom(x = s,
                         size = phi.moment,
                         mu = mu.moment)

  if ((moments.pmf[1] > 1e-5 & counts.start != 0) || moments.pmf[length(moments.pmf)] > 1e-5){
    warning("The density values at one or both of the ends of the given range are > 1e-5. Consider increasing the evaluated range.",
            .call = FALSE)
  }

  return(moments.pmf)
}
