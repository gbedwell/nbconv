#' Calculate the empirical sum of arbitrary NB random variables
#'
#'@param mus Vector of individual mean values
#'@param phis Vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param ps Vector of individual probabilities.
#'@param n.samp The number of samples of each distribution
#'
#'@export
#'
nb_empirical_sum <- function(mus, phis, ps, n.samp){

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


  make_dists <- function(x, y, n.samp){
    mat <- cbind(x, y)
    dist <- rnbinom(n=n.samp, size=mat[2], mu=mat[1])
  }

  dists <- mapply(FUN=make_dists, x = mus, y = phis, MoreArgs = list(n = n.samp))
  sumdists <- rowSums(dists)

  return(sumdists)
}
