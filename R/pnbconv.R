#' Calculates the CDF for the convolution of arbitrary negative binomial rv's
#'
#'@param counts The counts over which the convolution is evaluated. Should be a vector.
#'@param mus Vector of individual mean values
#'@param ps Vector of individual probabilities.
#'@param phis Vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param method The method by which to evaluate the PMF. One of "exact", "moments", or "saddlepoint".
#'@param n.terms The number of terms to include in the series for evaluating the PMF at a given number of counts. Defaults to 1000.
#'@param n.cores The number of CPU cores to use in the evaluation. Allows parallelization.
#'@param tolerance The acceptable difference between the sum of the K distribution and 1.
#'@param normalize Boolean. If TRUE, the PMF is re-normalized to sum to 1.
#'
#'@import parallel
#'@import matrixStats
#'
#'@export
#'
pnbconv <- function(targets, mus, ps, phis, method = c("exact", "moments", "saddlepoint"), 
                    n.terms = 1000, n.cores = 1, tolerance = 1e-3, normalize = TRUE){
  
  counts <- 0:max( counts )
  
  method <- match.arg(method, c("exact", "moments", "saddlepoint"))
  
  if( method != "exact" & method != "moments" & method != "saddlepoint"){
    stop("method must be one of 'exact', 'moments', or 'saddlepoint'.", call. = FALSE)
  }
  
  if (!missing(ps) & !missing(mus)){
    stop("'mus' and 'ps' both specified", call. = FALSE)
  }
  
  if (missing(ps) & missing(mus)){
    stop("One of 'mus' and 'ps' must be specified", call. = FALSE)
  }
  
  if( method == "exact" ){
    if (missing(ps) & !missing(mus)){
      ps <- phis / ( phis + mus )
    }
    if (length(ps) != length(phis)){
      stop("'ps' and 'phis' must have the same length.", call. = FALSE)
    }
  }
  
  if( method == "moments" | method == "saddlepoint"){
    if (missing(mus) & !missing(ps)){
      mus <- phis*(1 - ps)/ps
    }
    if (length(mus) != length(phis)){
      stop("'mus' and 'phis' must have the same length.", call. = FALSE)
    }
  }
  
  probs <- switch( method,
                   "exact" = nb_sum_exact( ps = ps, phis = phis, counts = counts, n.terms = n.terms, n.cores = n.cores, tolerance = tolerance ),
                   "moments" = nb_sum_moments( mus = mus, phis = phis, counts = counts ),
                   "saddlepoint" = nb_sum_saddlepoint( mus = mus, phis = phis, counts, normalize = normalize, n.cores = n.cores ) )
  
  cumprobs <- cumsum( probs )
  
  targetprobs <- cumprobs[targets + 1]
  
  return( targetprobs )
  
}
