#' Implements the method of moments approximation for the sum of arbitrary NB random variables. Called by other functions. Not intended to be run alone.
#'
#'@param mus Vector of individual mean values.
#'@param phis Vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param counts The vector of counts over which the PMF is evaluated.
#'
#'@returns A numeric vector of probability densities.
#'
#'@examples nb_sum_moments(mus = c(100, 10), phis = c(5, 8), counts = 0:500)
#'
#'@importFrom stats "dnbinom"
#'
#'@export
#'
nb_sum_moments <- function(mus, phis, counts){

  mu.moment <- sum( mus )
  phi.moment <- sum( mus )^2 / sum( mus^2 / phis )

  moments.pmf <- dnbinom(x = counts,
                         size = phi.moment,
                         mu = mu.moment)

  return(moments.pmf)
}
