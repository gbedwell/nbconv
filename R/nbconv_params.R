#' Summary statistics
#'
#' Calculates distribution parameters for the convolution of arbitrary negative binomial random variables.
#'
#'@param mus Numeric vector of individual mean values
#'@param phis Numeric vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param ps Numeric vector of individual probabilities of success.
#'
#'@returns A named numeric vector of distribution parameters.
#'
#'@examples nbconv_params(mus = c(100, 10), phis = c(5, 8))
#'
#'@import parallel
#'@import matrixStats
#'
#'@export
#'
nbconv_params <- function(mus, phis, ps){

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

  if( missing(mus) & !missing(ps) ){
    mus <- phis*(1 - ps)/ps
    }

  if( length(mus) != length(phis) ){
      stop("'mus' and 'phis' must have the same length.", call. = FALSE)
    }

  # First four cumulants expressed in terms of mu and phi
  k1 <- sum( mus )
  k2 <- sum( mus + mus^2 / phis )
  k3 <- sum( ( 2 * mus + phis ) * ( mus + phis ) * mus / phis^2 )
  k4 <- sum( ( 6 * mus^2 + 6 * mus * phis + phis^2)*( mus + phis ) * mus / phis^3 )

  # Central moments expressed in terms of cumulants
  # m1 <- k1
  # m2 <- k2
  # m3 <- k3
  # m4 <- k4 + 3 * k2^2

  mean <- k1
  sigma2 <- k2
  sigma <- sqrt( sigma2 )
  skewness <- k3 / k2^(3/2)
  ekurtosis <- k4 / k2^2

  pmax <- max( phis / ( phis + mus ) )
  qmax <- 1 - pmax

  K.mean <- ( mean * pmax / qmax ) - sum( phis )

  params <- c( mean = mean,
               sigma2 = sigma2,
               skewness = skewness,
               ex.kurtosis = ekurtosis,
               K.mean = K.mean )

  return( params )
}
