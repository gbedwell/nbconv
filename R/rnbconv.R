#' Generates random samples from the convolution of arbitrary negative binomial random variables.
#'
#'@param mus Vector of individual mean values
#'@param phis Vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.
#'@param ps Vector of individual probabilities of success.
#'@param n.samp The number of samples per distribution
#'
#'@returns A numeric vector of random deviates.
#'
#'@examples rnbconv(mus = c(100, 10), phis = c(5, 8), n.samp = 10)
#'
#'@importFrom stats "rnbinom"
#'@import parallel
#'
#'@export
#'
rnbconv <- function(mus, phis, ps, n.samp, n.cores = 1){

  if (!missing(ps) & !missing(mus)){
    stop("mus and ps both specified", call. = FALSE)
  }

  if (missing(ps) & missing(mus)){
    stop("One of mus and ps must be specified", call. = FALSE)
  }

  if ( !missing( ps ) & ( any( ps <= 0 ) | any( ps > 1 ) ) ){
    stop("ps must be 0 < ps <= 1", call. = FALSE)
  }

  if ( any( phis <= 0 ) ){
    stop("phis must be > 0.", call. = FALSE)
  }

  if ( !missing( mus ) & any( mus < 0 ) ){
    stop("mus must be > 0.", call. = FALSE)
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

  if (n.cores == 1){
    dists <- mapply(FUN=make_dists, x = mus, y = phis, MoreArgs = list(n = n.samp))
    sumdists <- rowSums(dists)
  } else {
    ind <- split( 1:length( phis ), ceiling( ( seq_along( 1:length( phis ) ) ) / 1000 ) )
    splitparam <- lapply( X = ind,
                          FUN = function(x){
                            matrix( data = c( mus[ x ], phis[ x ] ),
                                    ncol = 2 )
                            })

    distlist <- mclapply( X = splitparam,
                          FUN = function(z) {
                            zmus <- z[ , 1 ]
                            zphis <- z[ , 2 ]
                            dists <- mapply(FUN=make_dists, x = zmus, y = zphis,
                                            MoreArgs = list( n = n.samp ) )
                            sumdists <- rowSums(dists) },
                          mc.cores = n.cores )

    sumdists <- rowSums( do.call( cbind, distlist ) )
  }

  return(sumdists)
}
