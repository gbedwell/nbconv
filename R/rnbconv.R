#' Generates random samples from the convolution of arbitrary negative binomial rv's
#'
#'@param mus Vector of individual mean values
#'@param phis Vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param ps Vector of individual probabilities.
#'@param n.samp The number of samples per distribution
#'
#'@export
#'
rnbconv <- function(mus, phis, ps, n.samp){

  if (!missing(ps) & !missing(mus)){
    stop("'mus' and 'ps' both specified", call. = FALSE)
  }

  if (missing(ps) & missing(mus)){
    stop("One of 'mus' and 'ps' must be specified", call. = FALSE)
  }

  if (missing(mus) & !missing(ps)){
    mus <- phis*(1 - ps)/ps
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
