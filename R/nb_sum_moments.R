#' Method of Moments
#'
#' Implements the method of moments approximation for the sum of arbitrary NB random variables. Called by other functions. Not intended to be run alone.
#'
#'@param mus Numeric vector of individual mean values
#'@param phis Numeric vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param counts Integer vector of counts over which the convolution is evaluated.
#'
#'@returns A numeric vector of probability densities.
#'
#'@examples nb_sum_moments(mus = c(100, 10), phis = c(5, 8), counts = 0:500)
#'
#'@importFrom stats dnbinom
#'
#'
nb_sum_moments <- function(mus, phis, counts){

  mu.moment <- sum( mus )
  phi.moment <- sum( mus )^2 / sum( mus^2 / phis )

  moments.pmf <- dnbinom(x = counts,
                         size = phi.moment,
                         mu = mu.moment)

  return(moments.pmf)
}
