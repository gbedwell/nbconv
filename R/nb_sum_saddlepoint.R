#' Implements the saddlepoint approximation for the sum of arbitrary NB random variables. Called by other functions. Not intended to be run alone.
#'
#' Inspired by https://www.martinmodrak.cz/2019/06/20/approximate-densities-for-sums-of-variables-negative-binomials-and-saddlepoint/
#'
#'@param mus Vector of individual mean values.
#'@param phis Vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param counts The vector of counts over which the PMF is evaluated.
#'@param normalize Boolean. If TRUE, the PMF is normalized to sum to 1.
#'@param n.cores The number of CPU cores to use in the evaluation. Allows parallelization.
#'
#'@returns A numeric vector of probability densities.
#'
#'@examples nb_sum_saddlepoint(mus = c(100, 10), phis = c(5, 8), counts = 0:500)
#'
#'@import matrixStats
#'@import parallel
#'@importFrom stats "pnbinom"
#'@importFrom stats "uniroot"
#'
#'@export
#'
nb_sum_saddlepoint <- function(mus, phis, counts, normalize = TRUE, n.cores = 1){

  # Inspired by https://www.martinmodrak.cz/2019/06/20/approximate-densities-for-sums-of-variables-negative-binomials-and-saddlepoint/

  saddlepoint_calc <- function(mus, phis, counts){

    K  <-  function(t) { sum( phis * ( log( phis ) - log( phis + mus * ( 1 - exp( t ) ) ) ) ) }
    ldK <-  function(t) { logSumExp( log( phis ) + log( mus ) + t - log( phis + mus - mus * exp( t ) ) ) }
    lddK <- function(t) { logSumExp( log( phis ) + log( mus ) + log( phis + mus ) + t - 2 * log( phis + mus - mus * exp( t ) ) ) }
    pmf_eq <- function(t, x) { -0.5 * ( log( 2 * pi ) + lddK( t ) ) + K( t ) - t * x }

    if ( min( counts ) == 0 ){
      pmf0 <- prod( pnbinom( 0, size = phis, mu = mus ) )
    }

    counts <- counts[counts != 0]

    pmf <- sapply(X = counts,
                  FUN = function(x) {
                    t <- uniroot(function(t) { ldK(t) - log(x) },
                                 lower = -1e2,
                                 upper = min( log( phis / mus + 1 ) ),
                                 f.upper = min( log( phis / mus + 1 ) ),
                                 extendInt = "yes",
                                 tol = sqrt( .Machine$double.eps ) )$root
                    pmf <- pmf_eq(t, x)
                    return(pmf)
                  }
    )

    if (exists("pmf0")){
      pmf <- c( pmf0, exp( pmf ) )
    }
    else{
      pmf <- exp( pmf )
    }
    return(pmf)
  }

  if (n.cores == 1){
    saddlepoint.pmf <- saddlepoint_calc(mus = mus,
                                        phis = phis,
                                        counts = counts)
  }
  else{
    counts.list <- split( counts, ceiling( ( seq_along( counts ) ) / 1000 ) )

    pmf.list <- mclapply(X = counts.list,
                         FUN = function(y) {
                           pmf <- saddlepoint_calc(mus = mus,
                                                   phis = phis,
                                                   counts = y )
                           return(pmf) },
                         mc.cores = n.cores)

    saddlepoint.pmf <- Reduce(c, pmf.list)
  }

  if (isTRUE( normalize )){
    saddlepoint.pmf <- saddlepoint.pmf / sum( saddlepoint.pmf )
  }

  return( saddlepoint.pmf )
}
