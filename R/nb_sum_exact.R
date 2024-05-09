#' Furman's PMF
#'
#' Implements Furman's exact PMF for the evaluation of the sum of arbitrary NB random variables. Called by other functions. Not intended to be run alone.
#'
#'@param phis Numeric vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param ps Numeric vector of individual probabilities of success.
#'@param n.terms The number of terms to include in the series representation of the exact PMF. Defaults to 1000.
#'@param counts Integer vector of counts over which the convolution is evaluated.
#'@param n.cores The number of CPU cores to use in the evaluation. Allows parallelization.
#'@param tolerance The acceptable difference between the sum of the K distribution and 1. Defaults to 1E-3.
#'
#'@returns A numeric vector of probability densities.
#'
#'@examples nb_sum_exact(ps = c(0.05, 0.44), phis = c(5, 8), counts = 0:500)
#'
#'@import parallel
#'@import matrixStats
#'
#'
nb_sum_exact <- function(phis, ps, n.terms = 1000, counts, n.cores = 1, tolerance = 1e-3){
  # Implements the PMF described in https://ssrn.com/abstract=1650365
  # Adapted from https://github.com/slundberg/NBConvolution.jl/blob/master/src/furman.jl

  qs <- 1 - ps
  pmax <- max(ps)
  qmax <- 1 - pmax
  phisum <- sum(phis)

  logR <- sum( -phis * ( log( qs * pmax ) - log ( ps * qmax ) ) )

  xi <- c()
  xtmp <- c()
  for ( i in seq_len( n.terms ) ){
    xtmp <- ( log( phis ) + i * log( 1 - qmax * ps / ( qs * pmax ) ) ) - log( i )
    xi[i] <- logSumExp( xtmp )
  }

  delta <- c()
  dtmp <- c()
  delta[1] <- 0

  for (k in 1:( n.terms - 1 )){
    for ( i in seq_len( k ) ) {
      previndex <- k + 1 - i
      dtmp[i] <- log( i ) + xi[i] + delta[previndex]
    }
    delta[k + 1] <- logSumExp( dtmp ) - log( k )
  }

  logKdist <- logR + delta
  Ktest <- all.equal( 1, sum(exp(logKdist)), tolerance = tolerance )

  if (!isTRUE( Ktest )){
    stop( "The sum of the K distribution is insufficiently close to 1. Use more terms.",
          "\n",
          Ktest,
          call. = FALSE )
  }

  mass_calc <- function(x){
    total <- 0
    for ( k in seq( 0, (n.terms - 1) ) ){
      probs <- delta[k + 1] + ( lgamma( phisum + x + k ) - lgamma( phisum + k ) - lfactorial( x ) + ( phisum + k ) * log( pmax ) + x * log( qmax ) )
      total <- total + exp( probs )
    }
    masses <- log( total ) + logR
    return( masses )
  }

  if (n.cores == 1){
    pmf <- mass_calc(x = counts)
  } else {
    nn <- floor( length( counts ) / n.cores )
    count.list <- split( counts, ceiling( ( seq_along( counts ) ) / nn ) )

    clust <- makeCluster( n.cores, type = "PSOCK" )
    clusterExport( cl = clust,
                   varlist = list( "count.list", "phisum", "n.terms", "pmax", "qmax", "logR", "mass_calc" ),
                   envir = environment() )

    pmf.list <- parLapply(
      cl = clust,
      X = count.list,
      fun = function(x){
        pmf <- mass_calc(x)
        return(pmf)
        }
      )

    stopCluster(clust)

    pmf <- do.call( c, pmf.list )
  }
  if (is.numeric( pmf )){
    return( exp ( pmf ) )
  } else{
    error <- sub( "Error : *", "", pmf[1] )
    stop(error,
         "\n",
         call. = FALSE)
  }
}
