% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnbconv.R
\name{rnbconv}
\alias{rnbconv}
\title{Random Deviates}
\usage{
rnbconv(mus, phis, ps, n.samp, n.cores = 1)
}
\arguments{
\item{mus}{Numeric vector of individual mean values}

\item{phis}{Numeric vector of individual dispersion parameters. Equivalent to 'size' in dnbinom.}

\item{ps}{Numeric vector of individual probabilities of success.}

\item{n.samp}{The number of samples per distribution}

\item{n.cores}{The number of cores to use in the evaluation. Allows parallelization.}
}
\value{
A numeric vector of random deviates.
}
\description{
Generates random samples from the convolution of arbitrary negative binomial random variables.
}
\examples{
rnbconv(mus = c(100, 10), phis = c(5, 8), n.samp = 10)

}
