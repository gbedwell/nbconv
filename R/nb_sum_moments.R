#' Implements the method of moments approximation for the sum of arbitrary NB random variables
#'
#'@param mus Vector of individual mean values
#'@param phis Vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param counts The vector of counts over which the PMF is evaluated.
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
