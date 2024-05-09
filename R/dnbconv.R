#' Probability Mass Function
#'
#' Calculates the probability mass function for the convolution of arbitrary negative binomial random variables.
#'
#'@param counts Integer vector of counts over which the convolution is evaluated.
#'@param mus Numeric vector of individual mean values
#'@param ps Numeric vector of individual probabilities of success.
#'@param phis Numeric vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param method The method by which to evaluate the PMF. One of "exact", "moments", or "saddlepoint".
#'@param n.terms The number of terms to include in the series representation of the exact PMF. Defaults to 1000.
#'@param n.cores The number of CPU cores to use in the evaluation. Allows parallelization.
#'@param tolerance The acceptable difference between the sum of the K distribution and 1. Defaults to 1E-3.
#'@param normalize Boolean. If TRUE, the PMF is normalized to sum to 1.
#'
#'@returns A numeric vector of probability densities.
#'
#'@examples dnbconv(counts = 0:500, mus = c(100, 10), phis = c(5, 8), method = "exact")
#'
#'@import parallel
#'@import matrixStats
#'
#'@export
#'
dnbconv <- function(counts, mus, ps, phis, method = c("exact", "moments", "saddlepoint"),
                    n.terms = 1000, n.cores = 1, tolerance = 1e-3, normalize = TRUE){

  if( !is.integer( counts ) ){
    stop( "counts must be integer.",
          call. = FALSE )
  }

  method <- match.arg( method )

  if( method != "exact" & method != "moments" & method != "saddlepoint"){
    stop("method must be one of 'exact', 'moments', or 'saddlepoint'.", call. = FALSE)
  }

  if( !missing(ps) & !missing(mus) ){
    stop("mus and ps both specified", call. = FALSE)
  }

  if( missing(ps) & missing(mus) ){
    stop("One of mus and ps must be specified", call. = FALSE)
  }

  if( !missing( ps ) ){
    if( any( ps <= 0 ) | any( ps > 1 ) ){
      stop("ps must be 0 < ps <= 1", call. = FALSE)
    }
  }

  if( any( phis <= 0 ) ){
    stop("phis must be > 0.", call. = FALSE)
  }

  if( !missing( mus ) ){
    if( any( mus < 0 ) ){
      stop("mus must be > 0.", call. = FALSE)
    }
  }

  if( method == "exact" ){
    if(missing(ps) & !missing(mus)){
      ps <- phis / ( phis + mus )
    }
    if( length(ps) != length(phis) ){
      stop("ps and phis must have the same length.", call. = FALSE)
    }
  }

  if( method == "moments" | method == "saddlepoint"){
    if( missing(mus) & !missing(ps) ){
      mus <- phis * ( 1 - ps ) / ps
    }
    if( length(mus) != length(phis) ){
      stop("mus and phis must have the same length.", call. = FALSE)
    }
  }

  pmf <- switch( method,
                 "exact" = nb_sum_exact(
                   ps = ps,
                   phis = phis,
                   counts = counts,
                   n.terms = n.terms,
                   n.cores = n.cores,
                   tolerance = tolerance
                   ),
                 "moments" = nb_sum_moments(
                   mus = mus,
                   phis = phis,
                   counts = counts
                   ),
                 "saddlepoint" = nb_sum_saddlepoint(
                   mus = mus,
                   phis = phis,
                   counts,
                   normalize = normalize,
                   n.cores = n.cores
                   )
                 )

  return( pmf )
}

